# Portfolio Rebalancer Engine — Theoretical Reference

A self-contained mathematical specification of the engine: the convex
optimization problem it sets up, the algorithm that solves it, the data it
consumes and produces, and the extensions that slot into the same
framework. The document is intentionally language- and
implementation-agnostic; any code that satisfies the specification is a
correct engine. Read it as a chapter in a textbook rather than as an API
reference.

---

## Table of contents

1. What the system does
2. The math
3. The five cash-flow modes
4. Data flow
5. The engine internals
6. Inputs
7. Outputs
8. Constants reference
9. Edge cases
10. Extensions
11. Reproduction checklist
12. Mathematical concepts index

---

## 1. What the system does

A portfolio carries many exposures at once: currency mix, sector weights,
average beta, average duration, region split. A mandate sets targets on
several of them simultaneously, and the targets routinely contradict each
other — no set of trades hits all of them exactly.

The engine takes a portfolio and a set of target exposures and returns the
single trade vector that minimizes the combined, priority-weighted miss
across every target at once, under a chosen cash-flow rule. It also reports
precisely which targets it could not reach and by how much.

It is a pure function: same input, same output, every time. No state, no
randomness — the engine is a mathematical object, and what follows is its
specification.

---

## 2. The math

### 2.1 Variables

The portfolio has `n` holdings. Holding `i` has current dollar value
`w0_i`; the total is `V0 = sum w0_i`.

Each holding is traded through two non-negative variables: a buy `u_i >= 0`
and a sell `v_i >= 0`. The net dollar change in holding `i` is
`x_i = u_i - v_i` — positive a net buy, negative a net sell — and after
trading the holding is worth `w0_i + x_i`. The decision vector is
`y = [u; v]`, length `2n`.

The split exists because the transaction-cost term (§2.3) charges the
*gross* dollars traded, `|x_i|`, which is not differentiable in a single
signed variable. With the split, gross traded dollars are `u_i + v_i`, a
linear quantity. At any optimum with a positive cost on holding `i`, `u_i`
and `v_i` are never both positive — buying and selling the same name only
burns cost — so `|x_i| = u_i + v_i` holds exactly. Always read the trade as
`x_i = u_i - v_i`; never read `u_i` or `v_i` alone.

### 2.2 Deviation, in dollars

A target dimension `d` has an exposure `a_i,d` per holding (1/0 for a category
like "USD", or the raw number for a metric like beta) and a target level
`t_d`. The natural thing to control is a weight,

```
weight_d(x) = sum_i a_i,d (w0_i + x_i)  /  sum_i (w0_i + x_i)
```

but `x` appears in the denominator — minimizing a function of that ratio is
non-convex. Instead the engine controls the **dollar deviation**:

```
dev_d(x) = sum_{i in S_d} (a_i,d - t_d)(w0_i + x_i)
```

`S_d` is the dimension's scope (all holdings unless restricted). This is
linear in `x`. A linear quantity squared is a convex quadratic; the whole
problem becomes a convex QP with exactly one minimum, solved exactly and
fast. `dev_d = 0` is exactly `weight_d = t_d` for a categorical dimension.

Split it into a constant part and a trade part:

```
b_i,d = a_i,d - t_d   for i in S_d, else 0
dev_d(x) = e_d + b_d . x      where  e_d = b_d . w0
```

`e_d` is the deviation before any trade; `b_d . x` is what a trade changes.

Because `dev_d` is linear in `x` and `x_i = u_i - v_i`, the deviation is
linear in `y` with coefficient vector

```
beta_d = [ b_d ; -b_d ]        length 2n
```

and the same constant `e_d = b_d . w0`. Everything downstream is assembled
from `beta_d` over `y` exactly as it was from `b_d` over `x`.

### 2.3 The objective

Each dimension contributes a scaled squared deviation `f_d` exactly as
before — `lambda_d` its priority, `g_d` the dual normalizer of §2.3's
original `g_d` formula (unchanged).

The objective also prices trading. Each holding has a commission rate
`k_i >= 0`, a fraction of traded notional (§6.1 — supplied per holding,
or a flat default when absent). Gross dollars traded on holding `i`
are `u_i + v_i`, so the commission cost is

```
cost(y) = sum_i k_i (u_i + v_i)
```

The full objective is

```
minimize  F(y) = sum_d f_d(y)  +  theta * cost(y) / V0
```

`theta >= 0` is the cost-aversion weight. Dividing `cost` by `V0` keeps the
term scale-free, so the engine behaves identically across book sizes
(§5.5). `theta` is not a user-facing knob — it is the variable swept to
trace the cost/completion frontier (§5.6); the user-facing control is the
completion fraction `rho` (§5.7).

This term replaces the former flat turnover penalty `TAU`. With a single
commission rate `k_i = k` for every holding, `theta * k * sum(u_i+v_i)/V0`
is exactly the old turnover penalty, so the previous behavior is the
uniform-rate special case.

In matrix form the deviation part is `(1/2) y' H y + c' y` with `H` and `c`
accumulated from `beta_d` (§2.2) over the `2n` variables; the commission
part adds `theta k_i / V0` to the `c` entry of every `u_i` and every `v_i`.
`H` is positive semidefinite; the §2.3 ridge still applies after rescaling
and makes the minimizer unique.

### 2.4 Constraints

Non-negativity on every traded variable: `u_i >= 0`, `v_i >= 0`. A holding
cannot be sold for more than it is worth, so `v_i <= w0_i`. Buys are
unbounded above. Shorting, if enabled, raises the sell bound to
`w0_i + 4 V0`.

One linear constraint on net cash:
`sum_i (u_i - v_i) = rhs`, with `rhs` pinned or bounded by the mode (§3).

Modes restrict trade direction by pinning a block to zero: mode A
(sells only) pins every `u_i = 0`; mode D (buys only) pins every `v_i = 0`;
modes B, C, E pin neither. The problem stays a box-bounded QP with a single
linear equality — the structure the active-set solver (§2.6, §5.3) is built
for — only now over `2n` variables.

### 2.5 Optimality — KKT

At the optimum `x*` there are multipliers `nu` (on the cash equality) and
`z_i` (on the box bounds) with:

- **Stationarity:** `H x* + c - nu g - z = 0` (`g` is all-ones).
- **Primal feasibility:** `x*` satisfies all constraints.
- **Complementary slackness:** an interior variable has `z_i = 0`; a variable
  at its lower bound has `z_i <= 0`; at its upper bound, `z_i >= 0`.

The sign rule is the stopping test: when every pinned variable's multiplier
has the correct sign, no single move improves the objective, and because the
problem is convex this is the global optimum.

The KKT system is assembled over the `2n` variables `y`; the stationarity,
feasibility and complementary-slackness conditions are otherwise unchanged.

### 2.6 The solver — active-set

If the active set (which variables end at a bound) were known, the problem
would collapse to one small linear system. It is not known, so: guess, solve,
check the KKT signs, correct, repeat. Each correction improves the objective;
the number of possible active sets is finite, so it terminates — in practice
in a handful of steps. Full step-by-step procedure in §5.3.

The active-set routine is unchanged — it already handles any box-bounded,
single-equality QP; the cost model only enlarges the system from `n` to `2n`
variables.

### 2.7 A worked example

Three holdings: A = 60,000 (USD), B = 30,000 (EUR), C = 10,000 (USD).
`V0 = 100,000`. One target: USD share 0.50. Mode B (self-finance,
`sum x = 0`).

Currently USD is 70,000 / 100,000 = 70%, over the 50% target.
`b = [0.5, -0.5, 0.5]`, `e = b.w0 = 20,000` over, `s = 70,000`,
`g = 0.5(1/100000^2 + 1/70000^2) = 1.520e-10`. After rescaling the linear
term is `c = [40000, -40000, 40000]`.

The solver returns `x = [-10000, +20000, -10000]`: sell 10,000 each of the
two USD names, buy 20,000 of the EUR name. Resulting holdings
`[50000, 50000, 0]` — USD is now 50%, exactly on target. `cash_flow = 0`,
`turnover = 40,000`, `rms_deviation = 0`.

With one target the engine had freedom in how to split the USD selling
between A and C. The ridge tie-break picks the least-`x'x` split, which sells
C entirely. Add a second target and that freedom is spent satisfying it too
— the engine trades off every dimension at once.

In the `u/v` formulation the same solution reads `u = [0, 20000, 0]`,
`v = [10000, 0, 10000]`: holding A and C are pure sells, B a pure buy, no
holding both bought and sold. The net trade `x = u - v` is the
`[-10000, +20000, -10000]` above. The numbers are unchanged; the split is a
reformulation, not a different answer.

---

## 3. The five cash-flow modes

The mode is the only behavioral switch. It fixes the cash constraint and the
allowed trade signs; the objective is identical across all five.

| Mode | Cash constraint | Trade signs | Parameter |
|---|---|---|---|
| A withdraw | `sum x = -W` | `x_i <= 0` (sells only) | withdrawal `W > 0` |
| B self-finance | `sum x = 0` | free within box | none |
| C min inject | `sum x >= 0` | free | penalty `mu` |
| D deposit | `sum x = +D` | `x_i >= 0` (buys only) | deposit `D > 0` |
| E perfect | `sum x >= 0` | free | none |

- **A — Withdraw.** Raise exactly `W` in cash by selling; the engine chooses
  the liquidation that leaves the rest closest to target.
- **B — Self-finance.** No cash moves; sells fund buys. The routine
  rebalance.
- **C — Min injection.** Cash may be added but each dollar costs `mu`, so the
  engine adds the least it can. When reallocation alone suffices, it adds
  nothing and equals B.
- **D — Deposit.** Fresh money in, buys only, existing positions untouched.
- **E — Perfect.** Best fit at any cash cost (C with `mu = 0`). A miss that
  survives E is structural — unreachable from the current universe.

A residual that B leaves but E removes was a cash constraint. A residual that
survives E is a limit of the holdings themselves.

---

## 4. Data flow

The engine is a deterministic pipeline from inputs to outputs.

```
inputs:  portfolio,  target dimensions,  mode + parameters
   |
   v
validation
   |
   v
per-dimension construction
   exposure vector  a_d
   scope mask       m_d
   deviation coef.  b_d  (= a_d - t_d on scope, 0 off)
   constant         e_d  (= b_d . w0)
   normalizer       g_d  (dual scaling, §2.3)
   |
   v
QP assembly
   Hessian          H  (sum of dimension contributions + impact)
   linear term      c  (sum of dimension contributions + commission)
   box bounds       (lo, hi)
   cash constraint  sum_i (u_i - v_i) = rhs   (from the mode)
   |
   v
QP solve (active-set, §5.3)  ->  trade vector x
   |
   v
diagnose:  per-dimension residuals,  completion fraction rho,
           cost figures,  knee location,  ...
   |
   v
outputs:  trade vector,  diagnostic block
```

The pipeline runs once per `theta` to trace the cost/completion frontier
(§5.6), and once more at the final chosen `theta` — the knee (§5.7) or a
user override — to produce the delivered solution. The engine is a pure
function of its inputs: same inputs, same outputs.

---

## 5. The engine internals

### 5.1 Forming the QP

For each dimension build the constants of §2.2 — `b_d`, `e_d`, `s_d`, `g_d`
— and the `2n` deviation coefficient `beta_d = [b_d; -b_d]`. Accumulate

```
H += 2 g_d (beta_d outer beta_d)        (deviation Hessian)
c += 2 g_d e_d beta_d                   (deviation linear term)
```

over all dimensions. Then add the cost contributions: for every holding `i`,

```
c[u_i] += theta k_i / V0        (commission, buy block)
c[v_i] += theta k_i / V0        (commission, sell block)
```

and, when tier-2 impact data is present (§10.1), the four-block contribution
of `q_i (u_i - v_i)^2` to `H`. Rescale the assembled system by `max|H|` to
bring entries to order 1 (§5.5) and add the `1e-7` ridge to the diagonal of
`H` to make it strictly positive definite (§2.3).

Bounds: `u_i in [0, +inf)`, `v_i in [0, w0_i]` (the sell bound becomes
`[0, w0_i + 4 V0]` with shorts enabled). Mode A pins every `u_i = 0`; mode
D pins every `v_i = 0`. The cash constraint is `sum_i (u_i - v_i) = rhs`
per the mode table (§3).

### 5.2 Modes C and E

Their cash constraint is `sum x >= 0`. Solve box-only first (no equality); if
the result has `sum x >= -1e-6`, keep it; else re-solve once with `sum x = 0`
pinned.

### 5.3 The active-set routine

Inputs `H, c, lo, hi, has_eq, rhs`.

**Feasible start.** `x_i = clamp(0, lo_i, hi_i)`. If `has_eq`, project onto
the equality by a greedy sweep: with `resid = rhs - sum x`, walk holdings
adding room toward `hi` (if `resid > 0`) or `lo` (if `resid < 0`) until
`|resid| < 1e-12`.

**Working set** `W_i`: 0 free, -1 pinned low, +1 pinned high. All start 0.

**Iterate** up to `40n + 100` times:

1. Free set `F = {i : W_i = 0}`.
2. Gradient `grad = H x + c`.
3. Build the KKT system, size `d = |F| + (1 if has_eq)`: top-left block is
   `H` on `F`; if `has_eq`, the last row/column on the free part is 1 and the
   corner is 0; right-hand side is `-grad` on `F` (and 0 for the equality
   row). Solve it. The step `p` is `0` off `F` and the solution on `F`; `nu`
   is the equality multiplier (0 if none).
4. If `max|p| < 1e-9`: the free variables are EQP-optimal. For each pinned
   `i`, the bound multiplier is `z = grad_i - nu`; the violation is `-z` if
   pinned low, `z` if pinned high. Release the pinned variable with the
   largest violation if it exceeds `1e-9`; otherwise the KKT signs all hold —
   return `x`, the optimum.
5. Else ratio-test along `p`: largest `alpha` in `[0,1]` keeping every free
   variable inside its box. Take `x += alpha p`; if a variable blocked, pin
   it exactly at the bound and set `W`.

### 5.4 Post-processing

After the solver returns: if long-only, clamp `x_i = max(x_i, -w0_i)` to kill
sub-dollar drift. Then snap: `|x_i| < 40` becomes `x_i = 0`.

When the engine traces the cost/completion frontier (§5.6), the `$40` snap
is disabled for every sweep solve — snap is post-processing, not
optimization, and leaving it on puts spurious steps in the frontier. Snap
is applied once, to the final chosen solution only (§5.7). The long-only
clamp is a feasibility correction and stays on throughout.

### 5.5 Numerical scaling — why

Raw `H` carries a factor `~1/V0^2`, around `1e-11` for a million-dollar book
— working with numbers that small loses precision and makes the solver slow.
Rescaling by `max|H|` brings `H` to order 1 without moving the minimizer
(scaling an objective by a positive constant does not change its argmin).
This is why the engine behaves identically on a 100-thousand and a
100-million dollar portfolio.

### 5.6 The cost/completion frontier

A single solve of `F(y)` answers one `theta`. To expose the cost/quality
tradeoff the engine solves a ladder of `theta` values and records, for
each, the pair (commission cost, deviation achieved). `theta = 0` is the
cheapest-quality-ignored solve and reaches the best deviation the mode
allows; large `theta` drives trading toward zero.

Deviation is reported as a completion fraction. Let `J_dev` be the
deviation objective (the `f_d` sum, cost excluded). With `R_0 = J_dev` at
no trade and `R_perfect = J_dev` at the `theta = 0` solve,

```
rho = (R_0 - J_dev) / (R_0 - R_perfect)
```

`rho = 0` is no progress, `rho = 1` a perfect rebalance. `rho` is unit-free
and book-independent. It measures completion in *residual*, not turnover —
the early, cheap dollars of trading close most of the miss, so a high `rho`
can cost far less than a proportional share of full-rebalance turnover.
That gap is the diminishing-returns shape of the frontier and the reason a
knee exists (§5.7).

`R_0`, `R_perfect` and the frontier are all mode-specific. If
`R_0 - R_perfect < 1e-12` the portfolio is already on target — the engine
skips the sweep, returns the no-trade solution, and reports `rho = 1`.
Otherwise the `theta` ladder is anchored to the portfolio
(`anchor = (R_0 - R_perfect) / cost_max`, where `cost_max` is `cost(x)` at
the `theta = 0` solution) and refined where the frontier is sparsely
sampled; full procedure in `COST_SPEC.md` §S7.

### 5.7 The knee — the default operating point

The frontier is concave: each additional unit of completion costs more than
the last. Its knee is the point where the curve turns from steep to flat —
before it, cost buys substantial completion; after it, cost buys almost
none. The knee is the engine's default operating point: it is where the
marginal cost of further completion has risen to meet its marginal benefit,
so trading past it spends a dollar to save less than a dollar.

The knee is located by maximum distance to the chord. Normalize both
frontier axes to `[0,1]` (`rho` already is; `c_norm = cost / cost_max`,
where `cost_max` is the full-rebalance cost at `theta = 0`). Draw the chord
from the no-trade endpoint `(0,0)` to the full-rebalance endpoint `(1,1)`
— the chord is the rebalance's average cost per unit of completion. For
each frontier point the engine measures the perpendicular distance to that
chord,

```
d_j = ( rho_j - c_norm_j ) / sqrt(2)
```

and the knee is the point of greatest distance, `argmax_j d_j`.
Equivalently the knee is where the curve's local slope equals the chord's
slope — below the knee the curve completes faster than its average rate,
above it slower. The two views describe the same point.

The longest distance `d_max` doubles as the knee's sharpness. A pronounced
knee bows far from the chord; a gently bowed frontier with no real corner
hugs it. If `d_max` falls below `KNEE_MIN_LEG` (§8) the engine reports
**gradual tradeoff — no sharp optimum**, suppresses the knee marker,
defaults the recommendation to `rho = 0.90`, and prompts the user to
choose. An honest "no clear optimum" beats a confident marker on a flat
curve.

The knee is computed on the true (un-snapped) frontier. The final solution
— the knee, or a user override along `rho` — is solved once at the
corresponding `theta`, then snapped (§5.4); the engine reports the small
shift snapping causes so the delivered figures, not the idealized ones, are
what the user sees.

Because the frontier is the *cost-inclusive* one, the knee location depends
on the commission rates: a costlier book moves the knee toward less
trading, a cheap book toward more. This coupling is automatic.

---

## 6. Inputs

The engine consumes three pieces of structured information: a portfolio, a
list of target dimensions, and the choice of cash-flow mode with any
associated parameter.

### 6.1 The portfolio

A set of `n` holdings. Each holding `i` carries:

- a unique identifier,
- a current dollar value `w0_i >= 0` (with total `V0 = sum_i w0_i > 0`),
- zero or more *attribute values* — strings ("USD", "Tech") or numbers
  (beta = 1.2, duration = 7.4) — keyed by attribute names shared across
  holdings.

Attribute names are the link between holdings and target dimensions: a
target dimension on "currency" only makes sense if "currency" is one of the
attributes the holdings expose.

Optionally, a holding may carry cost data:

- a per-asset *commission rate* `k_i >= 0` (a fraction of traded notional)
  enables tier 1 of the cost model (§10.1);
- a per-asset *impact coefficient* `q_i >= 0`, or the pair (average daily
  volume, volatility) from which `q_i ~ c sigma_i / ADV_i` can be derived,
  enables tier 2.

In the absence of all cost data, the engine uses the flat `COST_DEFAULT`
rate (§8) for every holding — tier 0.

### 6.2 The target dimensions

A list of target dimensions. Each dimension carries:

- an *attribute name*, naming which of the portfolio's attributes the
  dimension is about;
- a *target level* `t_d` — a fraction in `[0,1]` for a categorical
  dimension, a number for a numeric one;
- a *match value* for categorical dimensions (the category whose weight is
  being targeted); absent for numeric dimensions;
- a *priority* `lambda_d >= 0` (default 1);
- optionally a *scope* — an (attribute, value) pair restricting the
  dimension to a sub-universe of holdings; absent means all holdings.

From these the engine constructs, for each dimension, the per-holding
exposure value `a_i,d` (1 if the holding's value on the named attribute
equals the match value, 0 otherwise, for categorical; the attribute's
number for numeric) and the scope mask `m_i,d`. The deviation coefficient
is then `b_i,d = a_i,d - t_d` on the scope (`m_i,d = 1`), and 0 off it,
exactly as in §2.2.

### 6.3 The mode and its parameter

The cash-flow mode (one of A–E, §3) and any associated parameter — a
withdrawal amount `W`, a deposit amount `D`, or an injection penalty `mu`
— selects the cash constraint and the trade-direction restrictions used in
the QP (§5.1).

---

## 7. Outputs

A solve produces a trade vector and a block of diagnostic quantities.

### 7.1 The trade vector

The signed trade vector `x = (x_1, ..., x_n)`, one entry per holding, in
dollars: positive a net buy, negative a net sell. Together with the
pre-trade values `w0` it determines the post-trade portfolio `w0 + x`.

Equivalently, the engine yields the split `(u, v)` of length `2n` for
which `x = u - v` (§2.1). At any solve with positive cost on a holding the
split is exact (`u_i v_i = 0`); the engine reports the signed `x`, never
`u` or `v` in isolation.

### 7.2 Diagnostic quantities

The diagnostic block summarises both the QP's solution and the math
problem's structure.

**Solve summary.**

- *Status* — optimal, or trivial when there are no target dimensions.
- *Cash flow* `sum_i x_i`, *turnover* `sum_i |x_i|`, *post-trade value*
  `V0 + sum_i x_i`.
- *Objective value*, split into the two halves of §2.3: the deviation
  half `sum_d f_d` and the cost half `theta * cost / V0`.

**Per-dimension residuals.** For each dimension `d`:

- the target `t_d`,
- the current exposure (pre-trade) and achieved exposure (post-trade),
- a percentage residual,
- the dollar residual `dev_d` (§2.2),
- the weighted residual `lambda_d g_d dev_d^2` (the dimension's
  contribution to the objective).

The percentage residual uses one of three formulas, applied in order:
`(achieved - t_d) / t_d` when the target is nonzero; else, for an unscoped
dimension, `achieved / (V0 + sum_i x_i)`; else `achieved / S_v` where
`S_v = sum_{i in S_d} (w0_i + x_i)`. Indeterminate forms collapse to 0.

**Aggregate residual quantities.**

- *RMS deviation* — root-mean-square of the percentage residuals across
  dimensions with `lambda_d > 0`;
- *Worst dimension* — the dimension carrying the largest weighted residual.

**Cost figures.**

- *Total cost* in dollars, and the same in basis points
  (`total_cost / V0 * 1e4`);
- *Per-trade cost* `k_i |x_i| + q_i x_i^2`;
- *Cost tier* (0 / 1 / 2 — §10.1) and a flag recording whether any holding
  used a defaulted rate;
- *Cost effect* — `selection` in the pinned-cash modes (A, D), where total
  turnover is fixed by the cash target and cost steers only *which*
  holdings trade; `selection+volume` in B, C, E, where cost also reduces
  total turnover.

**Frontier and knee.**

- The swept *frontier* — the set of `(rho, cost, J_dev)` points produced
  by the `theta` sweep (§5.6);
- The delivered solution's *achieved completion* `rho`;
- The *knee* — its `rho` and its perpendicular leg `d_max` (the sharpness
  measure, §5.7); the knee `rho` is null when the frontier has no sharp
  knee;
- The *snap delta* `(d_rho, d_cost)` — the shifts caused by post-processing
  snapping on the final solution (§5.4).

**Shorts.** When shorting is enabled, the identifiers of any holdings
whose post-trade value falls below zero.

---

## 8. Constants reference

| Constant | Value | Role | Where |
|---|---|---|---|
| ridge | `1e-7` | diagonal regularizer after rescale | §5.1 |
| snap | `40` | trades below this (abs $) zeroed | §5.4 |
| short bound | `4 * V0` | extra room below `-w0_i` with shorts | §5.1 |
| active-set tol | `1e-9` | step / multiplier zero-threshold | §5.3 |
| active-set cap | `40n + 100` | iteration limit | §5.3 |
| short flag | `1.0` | post-trade holding below this (negative) = short | §7 |
| `mu` presets | `1e-8` / `1.5e-7` / `2e-6` | mode C cash penalty (low/med/high) | §3 |
| `COST_DEFAULT` | `8e-4` (8 bps) | tier-0 flat commission rate (no cost data) | §10.1 |
| blank-cell percentile | 90th | default for missing cells of a present cost field | §6.1 |
| `theta` ladder | `10^-3..10^3`, 25 pts + 0 | initial frontier sweep grid | §5.6 |
| `theta` refine gap | `0.05` in `rho` | adaptive resampling trigger | §5.6 |
| sweep solve cap | 60 | max solves per frontier | §5.6 |
| `KNEE_MIN_LEG` | calibrated | min perpendicular distance for a sharp knee | §5.7 |

These are part of the specification, not tunable knobs — changing them
changes the output. `COST_DEFAULT`, the percentile, the ladder spacing and
the solve cap are specification. `KNEE_MIN_LEG` is a detector-sensitivity
threshold (a distance on the normalized unit square, bounded by
`1/sqrt(2)`), calibrated on real frontiers rather than derived.

---

## 9. Edge cases

- **No target dimensions.** The QP is trivial (`H = 0`, `c = 0`); the
  solver returns the zero trade vector and the engine reports a trivial
  status.
- **Empty scope.** A dimension whose scope mask is identically zero has
  `s_d = 0`; the normalizer falls back to whole-book scaling
  (`g_d = lambda_d / V0^2`) so the dimension still contributes, only
  unscaled by its (nonexistent) sleeve.
- **Zero-value holding.** With `w0_i = 0` the sell bound `v_i in [0,0]`
  pins `v_i = 0` automatically; the holding can be bought into in any mode
  that allows buys (B, C, D, E), cannot be sold from, and can only go
  negative if shorts are enabled.
- **Categorical match hitting nothing.** A target whose match value
  appears on no holding has the all-zero exposure vector; the dimension
  contributes a constant term `-t_d * V_d` to the deviation that the
  engine cannot move, so it is reported as a residual the trade cannot
  reach.
- **Already on target.** When `R_0 - R_perfect < 1e-12` the engine skips
  the frontier sweep and returns the no-trade solution with `rho = 1`
  (§5.6).
- **Degenerate buy/sell split.** At the sweep endpoint `theta = 0` the
  cost on a variable vanishes and `u_i, v_i` can both be positive — the
  ridge tie-break minimizes `y'y` and so suppresses the degeneracy
  numerically. The signed `x_i = u_i - v_i` is correct regardless; the
  engine reports the signed trade, never the split.
- **`V0 = 0`.** The whole problem is ill-posed (no portfolio to
  rebalance) and is rejected.
- **Invariant.** The post-trade total value equals the pre-trade total
  value plus the net cash flow exactly: `V0 + sum_i x_i`.

---

## 10. Extensions

The math is genericised so that several common extensions slot in without
changing the core. The most-used extension — the cost model — is set out
fully (§10.1); the framework for further extensions follows
(§10.2 – §10.5).

### 10.1 The cost model and its tiers

Transaction cost is part of the core engine (§2.3) and runs in three tiers,
each gated on what cost data is supplied per holding (§6.1). The QP stays
convex at every tier, so the active-set solver (§5.3) is unaffected.

- **Tier 0 — no cost data.** Every holding uses the flat `COST_DEFAULT`
  rate (§8). The cost term is in the objective, but a uniform `k` cancels
  out of the normalized chord (§5.7) and the knee location is identical to
  the no-cost baseline.
- **Tier 1 — per-asset commission.** A per-asset rate `k_i` enters `c`
  over the buy/sell-split variables: `c[u_i] += theta k_i / V0` and
  `c[v_i] += theta k_i / V0`. Linear, convex, no solver change.
- **Tier 2 — commission + impact.** A per-asset quadratic impact term
  `sum_i q_i x_i^2` enters `H`: the Hessian of `q_i (u_i - v_i)^2`
  populates the four `n x n` blocks linking `u_i` and `v_i`. Still convex.
  `q_i` may be supplied directly or derived as `q_i ~ c sigma_i / ADV_i`
  from volatility and average daily volume; absent both, the engine stays
  at tier 1. Full assembly in `COST_SPEC.md` §S2, §S5.

### 10.2 Adding a target dimension

A new target on any holding-level attribute enters as another
`(b_d, e_d, lambda_d, g_d)` tuple in the QP assembly (§5.1). The engine
is generic over the number and kind of dimensions; adding one does not
change the math, only its inputs.

### 10.3 Adding a cash-flow mode

A mode is defined by (a) the cash constraint — a linear equality or a
one-sided inequality on `sum_i (u_i - v_i)` — and (b) any sign restrictions
on the trade variables (`u_i = 0`, `v_i = 0`, or unrestricted). Any such
pair yields a box-bounded QP with a single linear equality, the structure
the active-set solver (§5.3) is built for; no change to the solver is
required.

### 10.4 Adding an objective term

A new soft objective term enters as additional contributions to `H` and
`c` in the QP assembly (§5.1). As long as the term is convex — its
Hessian contribution to `H` is PSD — the QP stays convex and the
active-set routine is unaffected. Linear terms enter only `c`; convex
quadratic terms enter only `H`. The cost model (§10.1) is exactly this —
its commission half is linear-in-`y`, its impact half is convex
quadratic-in-`y`.

### 10.5 What is out of scope

Risk-based objectives — return variance, tracking error against a
benchmark, drawdown — are quadratic-or-worse in the holdings and become
non-convex once shorting is admitted, so they would break the active-set
solver. Fixed per-ticket fees (`k_fixed * 1[x_i != 0]`) introduce a 0/1
indicator and are non-convex (cardinality); concave economies-of-scale
pricing is non-convex. None of these are admissible at the convex-QP
level; the post-processing snap (§5.4) is the pragmatic stand-in for
fixed-fee avoidance.

---

## 11. Reproduction checklist

A faithful rebuild satisfies all of:

1. Inputs as defined in §6 — a portfolio of holdings with attributes and
   optional cost data; a list of target dimensions with attribute, target
   level, optional match value, priority, optional scope; a mode with its
   parameter.
2. Per-dimension exposure `a_d` and scope mask `m_d` per §6.2; deviation
   coefficient `b_i,d = a_i,d - t_d` on `m_i,d = 1`, else 0; constant
   `e_d = b_d . w0`.
3. Dual-scaled normalizer `g_d = lambda_d * 0.5 * (1/V0^2 + 1/s_d^2)` with
   `s_d = sum_{i in S_d} w0_i` (or `V0` when the scope is empty), per §2.3.
4. Variables `y = [u; v]` of length `2n` per §2.1; deviation coefficient
   `beta_d = [b_d; -b_d]` per §2.2.
5. Deviation objective `H, c` over `y` per §5.1: accumulate
   `H += 2 g_d (beta_d outer beta_d)` and `c += 2 g_d e_d beta_d` across
   dimensions.
6. Cost terms by tier (§10.1): commission (tier 1) contributes
   `c[u_i] += theta k_i / V0` and `c[v_i] += theta k_i / V0`; impact
   (tier 2) contributes the `q_i (u_i - v_i)^2` Hessian to `H`. Then
   rescale by `max|H|` and add the `1e-7` ridge.
7. Bounds (`u_i >= 0`, `v_i in [0, w0_i]` or `[0, w0_i + 4 V0]` with
   shorts) and cash constraint `sum_i (u_i - v_i) = rhs` per mode (§3,
   §5.1); mode pinning of `u_i = 0` (A) or `v_i = 0` (D).
8. Active-set routine exactly as §5.3 — feasible-start projection, KKT
   assembly, multiplier sign test, ratio test — operating on the `2n`
   system.
9. C and E handled by solve-free-then-pin per §5.2.
10. Cost/completion frontier per §5.6: `theta` ladder anchored to the
    portfolio, adaptive refinement on `rho` gaps,
    `rho = (R_0 - J_dev)/(R_0 - R_perfect)`, sweep solve cap.
11. Knee detection per §5.7: chord-distance argmax with `d_max` as
    sharpness, `KNEE_MIN_LEG` gate for "no sharp knee".
12. Post-processing: long-only clamp throughout; snap (§5.4) applied only
    to the final chosen solution, not to sweep solves.
13. Diagnostic block per §7 — trade vector and the listed solve-summary,
    per-dimension residual, aggregate residual, cost, frontier-and-knee,
    and shorts quantities — including the three-way percentage-residual
    formula.
14. All edge cases of §9.
15. The worked example of §2.7 reproduced exactly: trade
    `x = [-10000, +20000, -10000]`, post-trade holdings `[50000, 50000, 0]`,
    USD share `0.50`.

### Scope and limits

The engine optimizes deviation from explicit exposure targets. It does not
optimize risk-based objectives — variance, tracking error, drawdown —
which are quadratic-or-worse and non-convex in the holdings and would
break the active-set solver. Targets are soft goals balanced by priority,
not hard constraints; the cash-flow rule of the chosen mode is the only
hard constraint, plus the no-short rule unless shorting is enabled.

Transaction cost is priced and optimized in tiers: a linear commission
term (tier 1) and a convex quadratic market-impact term (tier 2), each
gated on whether the corresponding per-holding data is supplied, with a
flat fallback at tier 0. Cost is a priced soft term, not a hard constraint
— the cash-flow rule of the mode remains the only hard constraint.

---

## 12. Mathematical concepts index

A reference for the moderate-to-higher-difficulty concepts the WIKI invokes.
One or two sentences each, organized by area so the index doubles as a
study path. Section references point back into this document.

### 12.1 Convex optimization and QP theory

**Convex optimization.** Minimizing a convex function over a convex feasible
set. The crucial consequence: any local minimum is the global minimum, so
the active-set routine can stop at the first KKT-satisfying point without
worrying about being stuck in a basin.

**Quadratic programming (QP).** A convex optimization problem whose
objective is quadratic in the variables, `(1/2) y' H y + c' y`, under linear
equality and inequality constraints. The rebalancer's per-`theta` solve is
exactly one QP (§2.3).

**Quadratic form.** A scalar expression `y' M y` for a symmetric matrix
`M`. The sign of `y' M y` for nonzero `y` classifies `M` as positive
(semi)definite, negative (semi)definite, or indefinite — which in turn
fixes whether the parent QP is strictly convex, convex, or non-convex.

**Positive semidefinite (PSD) and positive definite (PD).** A symmetric
matrix `H` is PSD if `y' H y >= 0` for every `y`, PD if the inequality is
strict for every nonzero `y`. PSD makes a QP convex; the `1e-7 * I` ridge
(§2.3) lifts it to PD so the minimizer is unique.

**Hessian matrix.** The matrix of second partial derivatives of a scalar
function. For a QP objective the Hessian is the matrix `H` itself, and
convexity of the objective is exactly the PSD property of `H` (§2.3,
§10.1).

**Box constraints.** Component-wise bounds `lo_i <= y_i <= hi_i`. The
active-set method in §5.3 is built specifically for problems that are box
plus a single linear equality — anything more general would change the
solver.

### 12.2 KKT conditions and the active set

**KKT (Karush–Kuhn–Tucker) conditions.** First-order conditions that any
local optimum of a smooth constrained problem must satisfy: stationarity,
primal feasibility, dual feasibility, and complementary slackness. For a
convex problem they are also *sufficient* — satisfying them means you've
found the global optimum (§2.5).

**Stationarity.** The condition that the gradient of the Lagrangian
vanishes at the optimum: in the rebalancer, `H y + c - nu * 1 - z = 0`.
Geometrically, no feasible direction strictly improves the objective.

**Primal and dual feasibility.** Primal feasibility: the solution `y*`
satisfies all original constraints. Dual feasibility: every inequality
multiplier has the correct sign — non-positive on a lower-bound multiplier,
non-negative on an upper-bound multiplier.

**Complementary slackness.** At the optimum, for every inequality
constraint either the constraint is binding or its multiplier is zero —
never both nonzero. The active-set method is engineered around this fact:
it guesses which constraints are binding and tests the multiplier signs.

**Lagrange multiplier (dual variable).** One scalar per constraint,
measuring how much the optimal objective would change if that constraint
were relaxed by one unit. Here `nu` rides the cash equality and `z_i` rides
each box bound.

**Active-set method.** A QP algorithm that maintains a guess of which
inequality constraints are binding (the "active set"), solves the
equality-constrained subproblem on that set, then releases or pins one
constraint at a time: release if a multiplier sign is wrong, pin if a
search step would violate a box. Each step strictly improves the objective
and termination is finite for convex QPs (§5.3).

**Equality-constrained QP (EQP).** A QP whose only constraints are linear
equalities. It collapses to a single dense linear system — the KKT system
— and is the subproblem the active-set method solves at every iteration.

**Ratio test.** Given a search direction `p`, find the largest scalar
`alpha in [0,1]` such that every free variable stays within its box after
`y + alpha p`. The blocking variable (first to hit a bound) is added to
the active set at that bound.

### 12.3 Numerical linear algebra

**Outer product.** Given column vectors `u`, `v`, the outer product `u v'`
is the rank-1 matrix with entries `(u v')_ij = u_i v_j`. The deviation
Hessian `H = sum_d 2 g_d (beta_d outer beta_d)` is a non-negative-weighted
sum of such rank-1 PSD matrices, which is why it is automatically PSD.

**Gradient.** The vector of partial derivatives of a scalar function; for
the QP objective `grad = H y + c`. The active-set routine uses it both to
generate search directions and to estimate the bound multipliers at pinned
variables (`z = grad_i - nu`).

**Gaussian elimination with partial pivoting.** A direct method for solving
a dense linear system `A x = b`: at each pivot step, swap the row with the
largest pivot entry into place, then eliminate below. Partial pivoting
prevents division by tiny numbers and gives backward-stable behaviour; it
is the linear-system step at the heart of each active-set iteration
(§5.3).

**Tikhonov (ridge) regularization.** Adding `epsilon * I` to a PSD Hessian
to make it strictly PD. The `(H + epsilon I) y = -c` system becomes
well-conditioned and, where the original problem had ties, the regularizer
picks the minimum-`y'y` representative — here `epsilon = 1e-7` is small
enough not to bias the answer at any digit a user reads (§2.3).

**Condition number.** A scalar that bounds how much rounding error in `b`
or `H` can amplify into error in the solution of `H y = b`. A condition
number near 1 is a stable solve; near `1 / eps_machine` the answer is
noise. The §5.5 rescale keeps the QP's condition number in a workable
range.

**Numerical rescaling.** Multiplying both sides of the QP by a positive
constant (here `1 / max|H|`) so the matrix entries are of order 1.
Rescaling does not move the minimizer — multiplying an objective by a
positive constant leaves its argmin unchanged — but it dramatically
improves floating-point conditioning.

### 12.4 Modelling techniques

**Variable splitting (buy/sell decomposition).** Replacing a signed
variable `x_i` with the difference of two non-negatives,
`x_i = u_i - v_i`, `u_i, v_i >= 0`. This linearizes the non-differentiable
`|x_i| = u_i + v_i` (exact at any optimum with positive cost on the
variable) and lets a turnover-like term enter the QP's linear part cleanly
(§2.1, §S4).

**Piecewise-linear convex function.** A convex function expressible as a
maximum of finitely many linear pieces, e.g. `|x| = max(x, -x)`. It is
convex but not differentiable at the kinks; to enter a smooth QP it has to
be linearized (e.g. variable splitting) or epigraphed first.

**Linearization of non-convex ratios.** Rewriting a ratio like
`weight_d(x) = sum a x / sum x` (non-convex in `x`) into a quantity that
is linear in `x`, here `dev_d(x) = sum (a - t) x`. Squaring a linear
quantity is a convex quadratic, which is admissible in a QP; the
rebalancer's convex form is built on this trick (§2.2).

**Scale-invariant normalization.** Dividing a problem term by a constant
proportional to the natural scale (here `V0`) so the term's numerical
value is unit-free. This keeps tuning constants like `theta` portable
across portfolios of different sizes and makes the engine's behaviour
book-independent (§2.3, §5.5).

**Convex relaxation.** Replacing a non-convex term with the closest convex
surrogate so the solver still applies. Fixed per-ticket fees
(`k_fixed * 1[x != 0]`) and concave economies-of-scale are non-convex and
out of scope here; the `$40` snap (§5.4) is a heuristic post-processing
stand-in, not a true relaxation.

### 12.5 Frontier analysis

**Pareto / efficient frontier.** The set of solutions to a multi-objective
problem that are not dominated by any other in every objective
simultaneously. The rebalancer's cost-vs-completion curve, swept by
`theta`, is one slice of such a frontier — each point is the best
deviation reachable at a given cost (§5.6).

**Concavity and diminishing returns.** A frontier `(rho, cost)` is concave
in `rho` when each additional unit of completion costs strictly more than
the last. This shape — the early dollars closing most of the miss, later
dollars closing only a sliver — is the engine's diminishing-returns curve
and the geometric reason a knee exists (§5.7).

**Chord.** The straight line connecting two distinguished points on a
curve. In the knee detector the chord runs from the no-trade endpoint
`(0,0)` to the full-rebalance endpoint `(1,1)` of the normalized frontier;
its slope is the rebalance's average cost-per-completion.

**Perpendicular distance to a chord (Kneedle algorithm).** For each curve
point, drop a perpendicular onto the chord and measure its length; the
knee is the point of maximum distance. With the chord running corner to
corner of the unit square the signed distance reduces to
`(rho - c_norm)/sqrt(2)`, and its maximum doubles as the knee's sharpness
`d_max` (§5.7).

**Argmin / argmax.** The *argument* — the input value — at which a
function attains its minimum (resp. maximum), as distinct from the
optimum's value. The rebalancer solves `argmin_y F(y)` for the trade
vector and locates the knee as `argmax_j d_j` on the swept frontier.
