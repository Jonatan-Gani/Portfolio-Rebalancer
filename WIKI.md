# Portfolio Rebalancer Engine — Wiki

The single reference for understanding the system and developing it further.
This file is self-contained: it absorbs the conceptual overview and the
reproduction methodology, and adds the architecture, data flow, file-by-file
internals, and a guide to extending the system.

Visual styling is described where it touches behavior, but pixel-level design
is not the focus — the focus is the engine, its contract, and how to build on
it.

---

## Table of contents

1. What the system does
2. The math
3. The five cash-flow modes
4. Architecture and data flow
5. File-by-file reference
6. The engine internals
7. The input contract — workbook schema
8. The output contract — diagnostics
9. Constants reference
10. Edge cases
11. Extending the system
12. Reference test cases
13. Reproduction checklist
14. Mathematical concepts index

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
network, no randomness.

Two ways to run it:

- **As a library.** `rebalancer.py` is the standalone engine. Import it, call
  `run(portfolio_df, targets_df, mode, params)`, get trades and diagnostics.
- **As an app.** `app.py` is a Dash web app. Upload an Excel workbook, pick a
  mode, drag a slider; results update live. The solver in the app is a
  JavaScript port of the same engine, running in the browser.

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
`k_i >= 0`, a fraction of traded notional (§7 — read from the portfolio
sheet, or a flat default when absent). Gross dollars traded on holding `i`
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
(§6.5). `theta` is not a user-facing knob — it is the variable swept to
trace the cost/completion frontier (§6.6); the user-facing control is the
completion fraction `rho` (§6.7).

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
linear equality — the structure the active-set solver (§2.6, §6.3) is built
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
in a handful of steps. Full step-by-step procedure in §6.3.

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

## 4. Architecture and data flow

The system has two execution paths that share one engine, expressed twice —
once in Python, once in JavaScript.

### 4.1 Library path

```
Excel workbook
   |
   |  read_sheets()           parse the two sheets
   v
portfolio_df, targets_df
   |
   |  build_dimensions()      validate, build Dimension objects
   v
pf, tk, val, dims
   |
   |  rebalance()             form QP, solve, diagnose
   v
{ trades, diagnostics }
```

`run()` chains these; `run_from_workbook()` adds the file read. All in
`rebalancer.py`, all pure Python.

### 4.2 App path

```
Browser: user uploads .xlsx
   |
   v
Dash server (app.py)
   |  read_sheets() + build_dimensions()      validate
   |  widget.serialize()                      portfolio+targets -> JSON
   |  widget.build_srcdoc()                   JSON + engine.js + widget_ui.js
   |                                          + CSS -> one HTML string
   v
<iframe srcDoc="...">      sent to the browser, once
   |
   v
Browser, inside the iframe:
   engine.js     the solver, ported from rebalancer.py
   widget_ui.js  controls, renders results
   |
   user drags slider / switches mode
   |  -> solve happens here, in the browser, no server round-trip
   v
results re-render instantly
```

The key architectural decision: **the server parses, the browser solves.**
The Dash server touches the workbook exactly once, to validate and serialize
it. Everything interactive — every slider tick, every mode switch — runs the
JavaScript engine locally in the iframe. This is why the slider is smooth: no
network is involved in a solve.

### 4.3 Why two copies of the engine

`rebalancer.py` is the reference engine and the library API. `engine.js` is a
line-for-line port for the browser. They implement the identical math and are
verified to agree numerically. If you change the math, you must change both —
see §11.

---

## 5. File-by-file reference

```
rebalancer.py     the engine + library API + CLI self-test   (Python)
app.py            Dash web server: upload, validate, serve    (Python)
widget.py         builds the iframe HTML (CSS + body + data)  (Python)
engine.js         the solver, ported for the browser          (JavaScript)
widget_ui.js      iframe controls and result rendering        (JavaScript)
make_example.py   generates the example workbook              (Python)
portfolio_example.xlsx   sample input, 30 assets
requirements.txt  Python dependencies
WIKI.md           this file
```

### rebalancer.py

The whole engine. Key functions:

- `_norm(s)` — normalize a string: lowercase, strip non-alphanumeric.
- `_pick(cols, *cands)` — match a column name against candidates (exact, then
  substring fallback).
- `read_sheets(src)` — load the workbook, return the two DataFrames.
- `Dimension` — dataclass: one target dimension with its precomputed exposure
  and scope vectors.
- `_portfolio_cols(pf)` — find the ticker and value columns.
- `build_dimensions(pf, tg)` — validate the targets sheet against the
  portfolio, build the `Dimension` list. Raises on any bad row.
- `_solve_as(H, c, lo, hi, has_eq, rhs)` — the active-set QP solver.
- `_diagnose(...)` — compute the diagnostics block from a trade vector.
- `rebalance(pf, tk, val, dims, mode, params)` — form the QP for a mode,
  solve, diagnose.
- `run(portfolio_df, targets_df, mode, params)` — top-level entry; returns
  `{ok, trades, diagnostics, ...}` or `{ok: False, error}`.
- `run_from_workbook(src, mode, params)` — `run` plus the file read.
- `make_example_workbook(path)` — writes a 10-asset CLI test workbook.
  NOTE: this overwrites `portfolio_example.xlsx`; the 30-asset example comes
  from `make_example.py` instead.
- `selftest()` — runs all five modes on the example and asserts validation;
  invoked by `python rebalancer.py`.

### app.py

The Dash server. Builds the page, handles the upload callback: parse the
workbook, validate, on success call `widget.build_srcdoc` and put the
returned HTML into the iframe; on failure show the error. The server does no
solving.

### widget.py

Turns a validated workbook into the iframe. `serialize(pf, tk, val, dims)`
flattens the portfolio and dimensions into a JSON payload (`holdings`,
`dims`, `cols`). `build_srcdoc(...)` assembles the full self-contained HTML:
CSS, the body markup, the serialized data, `engine.js`, and `widget_ui.js`.
`EMPTY_SRCDOC` is the placeholder shown before any upload.

### engine.js

The browser engine. `dimExposure` / `dimScope` build a dimension's vectors
from the raw data. `solveLin` is a dense linear solver (Gaussian elimination,
partial pivoting). `solveAS` is the active-set QP. `rebalance` forms the QP
and solves; `diagnose` builds the diagnostics. Mirrors `rebalancer.py`.

### widget_ui.js

The iframe UI. `buildTableHead` builds the merged portfolio+trades table
header from the data's column list. `run` solves and renders everything —
table, exposure bars, diagnostics. `solveCached` memoizes solves by
mode+amount+mu+shorts. `scheduleRun` coalesces rapid slider input to one
solve per animation frame. `warmup` pre-solves slider stops in idle time.
`setMode`, `buildParams`, `bindSeg` wire the controls.

### make_example.py

Generates `portfolio_example.xlsx`: 30 assets (equities and bonds), profiling
variables `asset_class, sector, currency, region, beta, duration`, and a
targets sheet with 10 dimensions.

---

## 6. The engine internals

### 6.1 Forming the QP

For each dimension build `b_d`, `e_d`, `s_d`, `g_d` exactly as before, then
the `2n` deviation coefficient `beta_d = [b_d; -b_d]`. Accumulate
`H += 2 g_d (beta_d outer beta_d)` and `c += 2 g_d e_d beta_d`.

Add commission: for every holding `i`, `c[u_i] += theta k_i / V0` and
`c[v_i] += theta k_i / V0`.

Rescale by `max|H|` and add the `1e-7` ridge as in the original §6.1, now
on the `2n` system.

Bounds: `u_i in [0, +inf)`, `v_i in [0, w0_i]` (sell bound `w0_i + 4 V0`
with shorts). Mode A pins all `u_i = 0`; mode D pins all `v_i = 0`. Cash
constraint `sum_i (u_i - v_i) = rhs` per the §3 mode table.

### 6.2 Modes C and E

Their cash constraint is `sum x >= 0`. Solve box-only first (no equality); if
the result has `sum x >= -1e-6`, keep it; else re-solve once with `sum x = 0`
pinned.

### 6.3 The active-set routine

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

### 6.4 Post-processing

After the solver returns: if long-only, clamp `x_i = max(x_i, -w0_i)` to kill
sub-dollar drift. Then snap: `|x_i| < 40` becomes `x_i = 0`.

When the engine traces the cost/completion frontier (§6.6), the `$40` snap
is disabled for every sweep solve — snap is post-processing, not
optimization, and leaving it on puts spurious steps in the frontier. Snap
is applied once, to the final chosen solution only (§6.7). The long-only
clamp is a feasibility correction and stays on throughout.

### 6.5 Numerical scaling — why

Raw `H` carries a factor `~1/V0^2`, around `1e-11` for a million-dollar book
— working with numbers that small loses precision and makes the solver slow.
Rescaling by `max|H|` brings `H` to order 1 without moving the minimizer
(scaling an objective by a positive constant does not change its argmin).
This is why the engine behaves identically on a 100-thousand and a
100-million dollar portfolio.

### 6.6 The cost/completion frontier

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
knee exists (§6.7).

`R_0`, `R_perfect` and the frontier are all mode-specific. If
`R_0 - R_perfect < 1e-12` the portfolio is already on target — the engine
skips the sweep, returns the no-trade solution, and reports `rho = 1`.
Otherwise the `theta` ladder is anchored to the portfolio
(`anchor = (R_0 - R_perfect) / cost_max`, where `cost_max` is `cost(x)` at
the `theta = 0` solution) and refined where the frontier is sparsely
sampled; full procedure in `COST_SPEC.md` §S7.

### 6.7 The knee — the default operating point

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
hugs it. If `d_max` falls below `KNEE_MIN_LEG` (§9) the engine reports
**gradual tradeoff — no sharp optimum**, suppresses the knee marker,
defaults the recommendation to `rho = 0.90`, and prompts the user to
choose. An honest "no clear optimum" beats a confident marker on a flat
curve.

The knee is computed on the true (un-snapped) frontier. The final solution
— the knee, or a user override along `rho` — is solved once at the
corresponding `theta`, then snapped (§6.4); the engine reports the small
shift snapping causes so the delivered figures, not the idealized ones, are
what the user sees.

Because the frontier is the *cost-inclusive* one, the knee location depends
on the commission rates: a costlier book moves the knee toward less
trading, a cheap book toward more. This coupling is automatic.

---

## 7. The input contract — workbook schema

One spreadsheet, two tables.

### 7.1 Locating the tables

Portfolio sheet: first whose normalized name contains `port`. Targets sheet:
first whose normalized name contains `target` or `tgt`. If not found by name,
fall back to order (first, second). Fewer than two sheets raises `workbook
needs two sheets: portfolio, targets`.

### 7.2 Column matching — the `pick` rule

Normalize (lowercase, strip non-alphanumeric), then: exact match against the
candidate list first; substring fallback second; absent if neither.

Portfolio candidates:
- ticker: `ticker, symbol, security, name, id, asset`
- value: `value, marketvalue, mv, amount, dollars, notional, weight`
- commission rate (optional, tier 1): `cost, costbps, commission, spread, tc, ticketcost`
- impact coefficient (optional, tier 2): `impact, impactcoef, marketimpact, qimpact`
- ADV (optional, used to derive impact): `adv, avgdailyvolume, advusd`
- volatility (optional, used to derive impact): `vol, volatility, sigma`

Targets candidates:
- label: `label, name, dimension, target, id`
- attribute: `attribute, column, field, factor, variable, var, exposure`
- match value: `matchvalue, value, category, bucket, level, member, class`
- target: `target, targetvalue, goal, weighttarget`
- lambda: `lambda, priority, lam, importance`
- scope attribute: `scopeattribute, scopecolumn, scopefield`
- scope value: `scopevalue, scopecategory, scopebucket`

Collision overrides: if the label column resolved to the same column as
target, label becomes absent; same for match value. Missing ticker/value, or
missing attribute/target, raises an error.

### 7.3 Row processing

Portfolio: drop rows with empty ticker; coerce value to numeric (unparseable
→ 0); every other column is a free-form profiling variable.

Cost columns are the exception to the "unparseable → 0" rule. A blank cost
cell takes the **90th percentile** of that column's known values, not 0 — a
0 would read as "this asset is free to trade" and route all flow through it.
When a cost column is absent entirely, every `k_i` defaults to
`COST_DEFAULT` (§9). A rate column given in basis points (`costbps`, or any
column whose values exceed 1) is divided by `1e4` on read so `k_i` is a
plain fraction. Full schema in `COST_SPEC.md` §S3.

Targets: each row is one dimension, unless its attribute cell is empty
(skipped). The attribute must name a portfolio column, else the row records
an error. Match value present → categorical (`unit = "share"`); absent →
numeric (`unit = "num"`). Target must parse to a number. Lambda defaults to
1.0. Label defaults to the match value or attribute name. Scope is optional;
a given scope attribute must be a portfolio column. One bad row rejects the
whole workbook — all errors are joined and raised together.

### 7.4 Exposure and scope vectors

Categorical exposure: 1.0 where the attribute column's stripped string equals
the match value (case-sensitive), else 0.0. Numeric exposure: the attribute
column coerced to number, unparseable → 0.0. Scope mask: all true if no
scope; else true where the scope column's stripped string equals the scope
value.

The example workbook's `targets` sheet uses columns `label, attribute,
match_value, target, lambda, scope_attribute, scope_value` — but any names
matching the candidate lists work.

---

## 8. The output contract — diagnostics

`run()` returns `{ok: True, trades, diagnostics, v0, n_holdings}` on success,
or `{ok: False, error}` on a validation failure.

`trades` — one entry per holding, portfolio order: `ticker, now (= w0_i),
dx (= x_i), after (= w0_i + x_i), side` (`buy` if `dx > 1e-6`, `sell` if
`< -1e-6`, else `hold`).

`diagnostics`:

- `solver_status` — `optimal`, or `trivial` when there are no dimensions.
- `mode` — the mode letter.
- `cash_flow = sum x_i`, `turnover = sum |x_i|`, `new_value = V0 + cash_flow`.
- `per_dim_residuals` — per dimension: `label, attribute, unit, target,
  current` (exposure before), `achieved` (exposure after), `residual_pct`,
  `residual_dollars` (the dollar deviation `dev`), `residual_weighted` (the
  dimension's weighted objective term), `lambda`.
- `rms_deviation` — RMS of `residual_pct` over dimensions with `lambda > 0`.
- `worst_dim` — label of the dimension with the largest
  `|residual_weighted|`.
- `shorts_opened` — tickers whose resulting holding is below `-1.0`.
- `objective_value` — the objective at the solution.

`residual_pct` uses one of three formulas: `(achieved - target)/target` if
the target is nonzero; else `achieved/new_value` if the dimension is
unscoped; else `achieved/Sv` (scoped, zero target). NaN/inf → 0.

The exact field names are part of the contract — `widget_ui.js` and any
integration bind to them.

The cost model adds to `diagnostics`: `total_cost` (cost paid, $),
`cost_bps` (`total_cost / V0 * 1e4`), `objective_deviation` and
`objective_cost` (the two terms of `objective_value`, reported
separately), `rho_achieved` (completion fraction of the delivered
solution), `knee_rho` (knee location, or null when no sharp knee),
`knee_leg` (`d_max`, the longest perpendicular distance — the sharpness
measure), `frontier` (the swept `(rho, cost, J_dev)` points),
`snap_delta` (`{ d_rho, d_cost }`, pre-snap vs post-snap shifts on the
final solution), `cost_tier` (`0` flat fallback, `1` commission only, `2`
commission + impact), and `cost_default_used` (bool — true if any asset
used a defaulted rate). Each `trades` entry gains `trade_cost`, the
per-asset cost `k_i|x_i| + q_i x_i^2`.

In pinned-cash modes (A, D) total traded volume is largely fixed by the
cash target, so commission steers *which* holdings trade, not how much; in
B, C, E it also reduces total turnover. The field `cost_affects`
(`"selection"` or `"selection+volume"`) records which.

---

## 9. Constants reference

| Constant | Value | Role | Where |
|---|---|---|---|
| ridge | `1e-7` | diagonal regularizer after rescale | QP build |
| snap | `40` | trades below this (abs $) zeroed | post-process |
| short bound | `4 * V0` | extra room below `-w0_i` with shorts | bounds |
| active-set tol | `1e-9` | step / multiplier zero-threshold | solver |
| active-set cap | `40n + 100` | iteration limit | solver |
| short flag | `1.0` | resulting holding below this (negative) = short | diagnostics |
| `mu` presets | `1e-8` / `1.5e-7` / `2e-6` | mode C cash penalty (low/med/high) | mode C |
| `COST_DEFAULT` | `8e-4` (8 bps) | flat commission rate when no cost column is present | QP build |
| blank-cell percentile | 90th | default for blank cells of a present cost column | input parsing |
| `theta` ladder | `10^-3..10^3`, 25 pts + 0 | initial frontier sweep grid | §6.6 |
| `theta` refine gap | `0.05` in `rho` | adaptive resampling trigger | §6.6 |
| sweep solve cap | 60 | max solves per frontier | §6.6 |
| `KNEE_MIN_LEG` | calibrated | min perpendicular distance for a sharp knee | §6.7 |

These are part of the specification, not tunable knobs — changing them
changes the output. `COST_DEFAULT`, the percentile, the ladder and the
solve cap are specification — changing them changes output. `KNEE_MIN_LEG`
is a detector-sensitivity threshold (a distance on the normalized unit
square, bounded by `1/sqrt(2)`), calibrated on real frontiers, not derived.
`TAU` (formerly the turnover-penalty weight) is superseded by the
commission term in §2.3 and no longer exists as a constant.

---

## 10. Edge cases

- **Zero dimensions** — valid workbook, no target rows: returns all-`hold`,
  `solver_status = "trivial"`.
- **Empty scope** — `s = 0`, so `s_eff = V0`; the dimension still enters with
  whole-book scaling only.
- **Zero-value holding** — in sells-only / self-finance it is pinned at 0 and
  cannot be funded; in deposit mode it can be bought into; with shorts it can
  go positive from 0.
- **Non-numeric value cells** — coerced to 0.
- **Non-numeric target** — rejects that row; **non-numeric lambda** — falls
  back to 1.0.
- **Categorical match hitting nothing** — exposure all-zero; the engine
  trades toward the (unreachable) target as far as the universe allows.
- **`V0 <= 0`** — raises `portfolio total value must be positive`.
- Invariant: `new_value` always equals `V0 + cash_flow` exactly.

---

## 11. Extending the system

### 11.1 Add a target dimension or a new metric

No code change. Add a row to the targets sheet; if it is a new metric, also
add the column to the portfolio sheet. The only rule: the targets-row
`attribute` must name a portfolio column. `build_dimensions` builds one more
`Dimension` and the solver includes it. This is the system working as
designed — it is generic over dimensions.

### 11.2 Change the math

If you change the objective, the constraints, or the solver, you must change
**both** `rebalancer.py` and `engine.js` — they are two copies of one engine.
After any such change, re-run the cross-check: solve the same workbook with
both and confirm the trade vectors and diagnostics agree (see §12). A change
to one without the other silently desynchronizes the library and the app.

### 11.3 Add a cash-flow mode

A mode is defined entirely by (a) the cash constraint and (b) the trade-sign
bounds. To add one: in `rebalance` (both files) add a branch that sets `lo`,
`hi`, and the cash constraint, then add the mode letter to the UI mode picker
in `widget_ui.js` and its description in `MDESC`. No solver change is needed
— the active-set routine already handles any box + single-equality QP.

### 11.4 Add a diagnostic field

Compute it in `_diagnose` (Python) and `diagnose` (JS), add it to the
returned `diagnostics` object in both, and render it in `widget_ui.js`'s
`run` (the diagnostics grid). Keep the field name identical across all three.

### 11.5 The cost model and its tiers

Transaction cost is part of the core engine (§2.3) and runs in three tiers,
each gated on what the portfolio sheet supplies. The QP stays convex at
every tier, so the active-set solver is unaffected.

- **Tier 0 — no cost column.** Every asset uses the flat `COST_DEFAULT`
  rate (§9). The cost term is in the objective but a uniform `k` cancels
  out of the normalized chord (§6.7), so the knee location is identical
  to the no-cost baseline.
- **Tier 1 — commission column present.** Per-asset `k_i` enters `c` over
  the buy/sell-split variables: `c[u_i] += theta k_i / V0` and
  `c[v_i] += theta k_i / V0`. Linear, convex, no solver change.
- **Tier 2 — commission + impact.** A per-asset quadratic impact term
  `sum_i q_i x_i^2` drops straight onto `H` (the Hessian of
  `q_i (u_i - v_i)^2` lives in the four `n x n` blocks linking `u_i` and
  `v_i`). Still convex, still no solver change. `q_i` can be read from an
  `impact` column directly or derived from ADV and volatility columns
  (`q_i ~ c sigma_i / ADV_i`); absent both, the engine stays at tier 1.
  Full assembly in `COST_SPEC.md` §S2, §S5.

A non-convex cost — fixed per-ticket fees (a 0/1 indicator on whether a
holding trades) or concave economies-of-scale pricing — would break the
active-set solver and is out of scope; the `$40` snap (§6.4) is the
pragmatic stand-in for fixed-fee avoidance.

### 11.6 Change the input format

The engine never sees Excel directly — `read_sheets` returns two DataFrames
and everything downstream is format-agnostic. To accept CSV pairs, a
database, or an API payload, replace `read_sheets` with something that
returns the same `(portfolio_df, targets_df)` shape. Nothing else changes.

### 11.7 Performance notes

The active-set solver is exact and fast for tens of holdings (single-digit
milliseconds). For hundreds of holdings the dense KKT solve (O(n^3) per
iteration) would dominate; a sparse factorization or a warm-started
factorization update would be the next step. The browser path additionally
memoizes solves (`solveCached`) and pre-warms slider stops (`warmup`), so
revisiting a slider position is instant.

---

## 12. Reference test cases

### 12.1 The 10-asset CLI workbook

`make_example_workbook()` writes this; `python rebalancer.py` solves it.
Portfolio columns `ticker, value, sector, fx, beta`; 10 holdings totalling
`V0 = 1,000,000`. Targets: USD 50%, EUR 50%, Tech 40%, Financials 40%, beta
1.00; all `lambda = 1`, no scope.

Expected (A withdraws 150,000; D deposits 200,000; C `mu = 1.5e-7`;
long-only):

```
A:  optimal,  cash -150,000,  rms ~20.9%
B:  optimal,  cash        0,  rms ~0.00%
C:  optimal,  cash        0,  rms ~0.00%
D:  optimal,  cash +200,000,  rms ~12.1%
E:  optimal,  cash        0,  rms ~0.00%
```

B, C, E coincide and inject nothing — the targets are reachable by
reallocation alone. A and D are one-directional and leave real residuals.

### 12.2 Validation

Point a targets-row attribute at a non-existent portfolio column;
`build_dimensions` must raise `target row K: attribute 'NAME' not in
portfolio columns [...]`.

### 12.3 The three-asset worked example

§2.7: holdings 60k/30k/10k, one USD-share-0.50 target, mode B. Solution
`x = [-10000, +20000, -10000]`, holdings `[50000, 50000, 0]`, USD share
exactly 0.50, `cash_flow = 0`, `turnover = 40,000`, `rms_deviation = 0`.

### 12.4 Python / JavaScript agreement

Solve any workbook with `rebalancer.py` and with `engine.js` and confirm the
trade vectors agree to within a dollar per holding and the diagnostics to 4
significant figures. This is the regression test after any math change.

---

## 13. Reproduction checklist

A faithful rebuild satisfies all of:

1. Two-table input located by the name/order rule (§7.1).
2. Columns identified by the exact `pick` rule and candidate lists (§7.2),
   including the two `target`-collision overrides.
3. Row processing with the exact skip/coerce/default rules (§7.3); one bad
   row rejects the whole workbook.
4. Exposure and scope vectors per §7.4 (categorical = case-sensitive stripped
   equality; numeric = coerced number).
5. Variables `y = [u; v]` of length `2n` per §2.1; objective `H, c` over `y`
   per §6.1 with the dual `g` scaling and the `beta_d = [b_d; -b_d]` stack.
6. Cost terms by tier (§11.5): commission (tier 1) contributes
   `c[u_i] += theta k_i / V0` and `c[v_i] += theta k_i / V0`; impact
   (tier 2) contributes the `q_i (u_i - v_i)^2` Hessian onto `H` per
   `COST_SPEC.md` §S5. Then rescale by `max|H|` and add the `1e-7` ridge.
7. Bounds (`u_i >= 0`, `v_i in [0, w0_i]` or `[0, w0_i + 4 V0]` with shorts)
   and cash constraint `sum_i (u_i - v_i) = rhs` per mode (§3, §6.1); mode
   pinning of `u_i = 0` (A) or `v_i = 0` (D).
8. Active-set solver exactly as §6.3 — feasible-start projection, KKT
   assembly, multiplier sign test, ratio test — operating on the `2n`
   system.
9. C and E handled by solve-free-then-pin (§6.2).
10. Cost/completion frontier per §6.6: `theta` ladder anchored to the
    portfolio, adaptive refinement on `rho` gaps,
    `rho = (R_0 - J_dev)/(R_0 - R_perfect)`, sweep solve cap.
11. Knee detection per §6.7: chord-distance argmax with `d_max` as
    sharpness, `KNEE_MIN_LEG` gate for "no sharp knee".
12. Post-processing: long-only clamp throughout; snap (§6.4) applied only to
    the final chosen solution, not to sweep solves.
13. Diagnostics per §8 with the exact field names — including `total_cost`,
    `cost_bps`, `objective_deviation`, `objective_cost`, `rho_achieved`,
    `knee_rho`, `knee_leg`, `frontier`, `snap_delta`, `cost_tier`,
    `cost_default_used`, `cost_affects`, and per-trade `trade_cost` — and
    the three-way `residual_pct` rule.
14. All edge cases of §10.
15. The §12 test cases reproduced to the stated precision.

### Scope and limits

The engine optimizes deviation from explicit exposure targets. It does not
optimize risk-based objectives — variance, tracking error, drawdown — which
are quadratic-or-worse and non-convex in the holdings and would break the
active-set solver. Targets are soft goals balanced by priority, not hard
constraints; the cash-flow rule of the chosen mode is the only hard
constraint, plus the no-short rule unless shorting is enabled.

Transaction cost is priced and optimized in tiers: a linear commission
term (tier 1) and a convex quadratic market-impact term (tier 2), each
gated on the corresponding portfolio column, with a flat fallback at
tier 0. Cost is a priced soft term, not a hard constraint — the cash-flow
rule of the mode remains the only hard constraint.

---

## 14. Mathematical concepts index

A reference for the moderate-to-higher-difficulty concepts the WIKI invokes.
One or two sentences each, organized by area so the index doubles as a
study path. Section references point back into this document.

### 14.1 Convex optimization and QP theory

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
§11.5).

**Box constraints.** Component-wise bounds `lo_i <= y_i <= hi_i`. The
active-set method in §6.3 is built specifically for problems that are box
plus a single linear equality — anything more general would change the
solver.

### 14.2 KKT conditions and the active set

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
and termination is finite for convex QPs (§6.3).

**Equality-constrained QP (EQP).** A QP whose only constraints are linear
equalities. It collapses to a single dense linear system — the KKT system
— and is the subproblem the active-set method solves at every iteration.

**Ratio test.** Given a search direction `p`, find the largest scalar
`alpha in [0,1]` such that every free variable stays within its box after
`y + alpha p`. The blocking variable (first to hit a bound) is added to
the active set at that bound.

### 14.3 Numerical linear algebra

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
is what `solveLin` in `engine.js` implements.

**Tikhonov (ridge) regularization.** Adding `epsilon * I` to a PSD Hessian
to make it strictly PD. The `(H + epsilon I) y = -c` system becomes
well-conditioned and, where the original problem had ties, the regularizer
picks the minimum-`y'y` representative — here `epsilon = 1e-7` is small
enough not to bias the answer at any digit a user reads (§2.3).

**Condition number.** A scalar that bounds how much rounding error in `b`
or `H` can amplify into error in the solution of `H y = b`. A condition
number near 1 is a stable solve; near `1 / eps_machine` the answer is
noise. The §6.5 rescale keeps the QP's condition number in a workable
range.

**Numerical rescaling.** Multiplying both sides of the QP by a positive
constant (here `1 / max|H|`) so the matrix entries are of order 1.
Rescaling does not move the minimizer — multiplying an objective by a
positive constant leaves its argmin unchanged — but it dramatically
improves floating-point conditioning.

### 14.4 Modelling techniques

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
book-independent (§2.3, §6.5).

**Convex relaxation.** Replacing a non-convex term with the closest convex
surrogate so the solver still applies. Fixed per-ticket fees
(`k_fixed * 1[x != 0]`) and concave economies-of-scale are non-convex and
out of scope here; the `$40` snap (§6.4) is a heuristic post-processing
stand-in, not a true relaxation.

### 14.5 Frontier analysis

**Pareto / efficient frontier.** The set of solutions to a multi-objective
problem that are not dominated by any other in every objective
simultaneously. The rebalancer's cost-vs-completion curve, swept by
`theta`, is one slice of such a frontier — each point is the best
deviation reachable at a given cost (§6.6).

**Concavity and diminishing returns.** A frontier `(rho, cost)` is concave
in `rho` when each additional unit of completion costs strictly more than
the last. This shape — the early dollars closing most of the miss, later
dollars closing only a sliver — is the engine's diminishing-returns curve
and the geometric reason a knee exists (§6.7).

**Chord.** The straight line connecting two distinguished points on a
curve. In the knee detector the chord runs from the no-trade endpoint
`(0,0)` to the full-rebalance endpoint `(1,1)` of the normalized frontier;
its slope is the rebalance's average cost-per-completion.

**Perpendicular distance to a chord (Kneedle algorithm).** For each curve
point, drop a perpendicular onto the chord and measure its length; the
knee is the point of maximum distance. With the chord running corner to
corner of the unit square the signed distance reduces to
`(rho - c_norm)/sqrt(2)`, and its maximum doubles as the knee's sharpness
`d_max` (§6.7).

**Argmin / argmax.** The *argument* — the input value — at which a
function attains its minimum (resp. maximum), as distinct from the
optimum's value. The rebalancer solves `argmin_y F(y)` for the trade
vector and locates the knee as `argmax_j d_j` on the swept frontier.
