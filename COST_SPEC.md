# Cost-Aware Rebalancing — Design Spec

Companion to `WIKI.md`. Adds a transaction-cost model and a completion-fraction
knob with an auto-located economic optimum (the "knee"). Section references
(§n) point into `WIKI.md` unless prefixed `S` (this file).

Notation follows `WIKI.md`: `w0_i`, `V0`, `x_i`, `dev_d`, `e_d`, `b_d`, `H`,
`c`, `g_d`, `lambda_d`, modes A–E. New symbols introduced where used.

---

## S1. Summary

Two changes, one new control:

1. The flat `TAU` turnover term (§8) is replaced by a real per-asset cost
   model — linear commission plus optional quadratic impact. Each term is
   gated on an optional portfolio column; absent column → term off → current
   behavior. Commission is tier 1, impact tier 2.
2. The objective gains a cost-aversion weight `theta`. `theta` is not a
   user-facing opinion knob. It is the internal axis of a frontier sweep.
3. The user-facing control is a completion fraction `rho in [0,1]`: 0 = no
   trade, 1 = full rebalance. The engine traces the cost-vs-`rho` frontier,
   auto-locates the knee (the economically optimal point), solves there by
   default, and exposes a knob for the PM to override.

The active-set solver (§5.3) is structurally unchanged — still box bounds
plus one linear equality. The QP *assembly* in `rebalance` is rewritten
(buy/sell variable split). Both `rebalancer.py` and `engine.js` change; they
must stay in sync.

---

## S2. Cost model

Total cost of a trade vector `x`:

```
cost(x) = sum_i k_i |x_i|        (commission / spread, linear)
        + sum_i q_i x_i^2        (market impact, quadratic)
```

- `k_i >= 0` — per-asset commission rate, a fraction of traded notional.
  `k_i |x_i|` is dollars.
- `q_i >= 0` — per-asset impact coefficient, units chosen so `q_i x_i^2` is
  dollars.

Both terms are convex, so the QP stays convex (§10.1). The linear term is
piecewise-linear, handled by the buy/sell split (S4). The quadratic term
drops into `H` directly (S5).

### S2.1 Data degradation ladder

Each term is gated on its source column. Behavior by tier:

| Tier | Columns present | Active terms |
|---|---|---|
| 0 | none | flat fallback rate only — equals current `TAU` behavior |
| 1 | commission | linear commission |
| 2 | commission + impact | linear + quadratic |

"Always start from commission": linear is the baseline cost model; impact is
opt-in. Impact is a derived quantity (`q_i ~ c * sigma_i / ADV_i`); few
clients supply it ready-made. If a client supplies ADV and volatility columns
but no `impact_coef`, derive `q_i`; if neither, tier stays 1, nothing breaks.

---

## S3. Column schema additions

Extends §6.1. New portfolio-sheet columns, matched by the existing `pick`
rule (normalize, exact-then-substring).

- commission rate — candidates: `cost, costbps, commission, spread, tc,
  ticketcost`
- impact coefficient — candidates: `impact, impactcoef, marketimpact, qimpact`
- (optional, for deriving impact) ADV — `adv, avgdailyvolume, advusd`;
  volatility — `vol, volatility, sigma`

If a rate column is given in basis points (`costbps`), divide by `1e4` on
read so `k_i` is a plain fraction. Detect by column name suffix or by
magnitude (values > 1 are almost certainly bps).

### S3.1 Missing-value defaults — the partial-column trap

A blank cost cell must NOT default to 0. Zero reads as "this asset is free to
trade," and the engine routes all flow through it — actively wrong, worse
than no cost model. Defaults:

- Column entirely absent → flat fallback rate for every asset. Fallback =
  `COST_DEFAULT` (S13). Knee-safe (S8.3).
- Column present, some cells blank → blank cells take the **90th percentile**
  of the column's known values. Conservative-high: cheap-to-trade must be
  proven by data, never assumed.

---

## S4. Buy/sell variable split

The linear term `k_i |x_i|` is non-differentiable, cannot enter `H`. Split
each holding into a non-negative buy and a non-negative sell.

```
x_i = u_i - v_i ,    u_i >= 0 ,    v_i >= 0
```

At any optimum with `k_i > 0`, `u_i` and `v_i` are never both positive
(buying and selling the same name only burns cost for zero net effect), so
`|x_i| = u_i + v_i` exactly. The decision vector becomes `y = [u; v]`, length
`2n`.

Exception: at the sweep endpoint `theta = 0` the cost term vanishes and the
split can be degenerate (`u_i, v_i` both nonzero). `x_i = u_i - v_i` is still
correct — always read `x` from the difference, never read `u, v`
individually. The ridge tie-break (§2.3) minimizes `y'y`, which already
discourages the degenerate split.

### S4.1 Mapping the deviation objective onto `y`

`dev_d` is linear in `x`, hence linear in `y`. The deviation coefficient
vector over `y` is

```
beta_d = [ b_d ; -b_d ]        length 2n
```

with `e_d` unchanged. Assemble `H_y` (`2n x 2n`) and `c_y` (length `2n`) by
the **existing** §5.1 accumulation, substituting `beta_d` for `b_d`. No new
formula — the same accumulation over the doubled coefficient vector.

### S4.2 Bounds and cash constraint over `y`

```
u_i in [0, +inf)
v_i in [0,  w0_i]           ( [0, w0_i + 4*V0] with shorts enabled )
cash:  sum_i (u_i - v_i) = rhs        (rhs per mode, §3)
```

Modes via pinning, replacing the §5.1 sign bounds:

- Mode A (sells only): pin `u_i = 0` for all `i`.
- Mode D (buys only): pin `v_i = 0` for all `i`.
- Modes B/C/E: no extra pinning.

Still box + one equality — the active-set routine (§5.3) is untouched, only
larger (`2n`).

---

## S5. Objective with cost

```
F(y) = J_dev(y) + theta * cost(y) / V0
```

- `J_dev(y)` — the deviation objective, `H_y`/`c_y` from S4.1.
- `theta >= 0` — cost-aversion weight, the sweep variable (S7). Not
  user-facing.
- `/ V0` — mandatory normalization. Keeps the cost term scale-free so the
  engine behaves identically across book sizes (§5.5). `k_i` and `q_i x_i^2 /
  V0` are already fractions; dividing by `V0` makes the whole term
  dimensionless and `theta` portable.

Cost contributions into the QP:

- Commission (linear → `c_y`):
  `c_y[i]   += theta * k_i / V0`   (buy block)
  `c_y[n+i] += theta * k_i / V0`   (sell block)
- Impact (quadratic → `H_y`), with `q'_i = theta * q_i / V0`:
  `H_y[i,i]     += 2 q'_i`
  `H_y[n+i,n+i] += 2 q'_i`
  `H_y[i,n+i]   += -2 q'_i`
  `H_y[n+i,i]   += -2 q'_i`
  (this is the Hessian of `q_i (u_i - v_i)^2`).

A single `theta` scales total cost. The commission/impact mix is fixed by the
data (`k_i`, `q_i` are measured rates), not by `theta`.

The old `TAU` term is removed. At tier 0, `k_i = COST_DEFAULT` flat, which
reproduces the `TAU` tiebreaker role.

Rescale and ridge (§5.1) proceed unchanged on the `2n` system.

---

## S6. The completion axis `rho`

Scalar residual measure: `J_dev` (the deviation objective value, cost term
excluded). Per mode:

- `R_0` — `J_dev` at no trade (`x = 0`).
- `R_perfect` — `J_dev` at the full rebalance for that mode (the `theta = 0`
  solution).

Completion fraction of any solution:

```
rho = (R_0 - J_dev) / (R_0 - R_perfect)
```

`rho = 0` no progress, `rho = 1` perfect. `rho` is unit-free and
book-independent — `rho = 0.85` means the same thing on a $10M and a $400M
book.

`rho` measures completion in **residual, not turnover**. 85% of the miss
closed may cost far less than 85% of full-rebalance turnover — the cheap
early dollars do most of the work. That gap is the diminishing-returns curve
and the reason the knee exists. Defining `rho` on turnover would straighten
the curve and destroy the knee.

`R_perfect` is mode-specific; `rho` and the knob are always relative to the
currently selected mode. If `R_0 - R_perfect < 1e-12` (already on target),
skip the sweep, return the no-trade solution, report `rho = 1`.

---

## S7. The frontier sweep

For each `theta`, solve `F(y)` (S5) and record `(cost, J_dev) -> (cost,
rho)`. `theta = 0` gives `rho = 1`; `theta -> inf` gives `rho -> 0`.

**Penalty route, not budget route.** A direct `rho` target would be a
quadratic inequality constraint — outside the solver's box+one-equality
structure. Sweeping `theta` instead means every solve is the standard QP
(S5). The PM-facing axis is still `rho`; `theta` is purely internal.

### S7.1 Auto-scaled `theta` ladder

`theta` magnitude is not intuitive (cost term ~bps, deviation term ~O(0.01-1)).
Anchor it to the portfolio:

```
anchor = (R_0 - R_perfect) / cost_max
```

`cost_max` = `cost(x)` of the `theta = 0` solution. Ladder:

```
theta_k = anchor * 10^t ,   t in linspace(-3, +3, 25)
```

plus `theta = 0`. Then **refine adaptively**: where consecutive solves are
more than `0.05` apart in `rho`, bisect `theta` in that gap and re-solve.
Stop when no gap exceeds `0.05` or a solve cap (`60`) is hit. Uneven `theta`
sampling bunching points near the knee is desirable — that is where the
detector needs resolution to place the longest perpendicular leg accurately.

### S7.2 Snap disabled during the sweep

The `$40` snap (§5.4) is **off** for every sweep solve. Snap is
post-processing, not optimization; leaving it on makes trades cross the $40
threshold as `theta` moves, putting stair-steps in the curve that distort the
perpendicular distances and produce fake knees. With snap off, the swept
frontier is the true QP frontier — smooth, because the QP solution varies
continuously in `theta`. No curve smoothing or spline fitting is needed.

The long-only clamp (§5.4) stays on during the sweep — it is a feasibility
correction, not cosmetic.

---

## S8. Knee detection

The frontier is concave: each unit of completion costs more than the last.
The knee is the point where the curve stops being steep and goes nearly
flat — before it, cost buys a lot of completion; after it, cost buys almost
none. The knee is the economically optimal default: it is the point where the
marginal cost of one more unit of completion has risen to meet its marginal
benefit.

### S8.1 Normalize both axes

The knee's location depends on how the two axes are scaled, so both are
normalized to `[0,1]` before detection:

```
rho       already in [0,1]
c_norm    = cost / cost_max
```

`cost_max` = `cost(x)` of the `theta = 0` (full-rebalance) solution.
Normalizing also makes a flat tier-0 commission default knee-invariant
(S8.4): a flat `k` is a global factor that cancels in `cost / cost_max`.

### S8.2 Distance-to-chord detector

Detection is the maximum-distance-to-chord method (Kneedle).

Draw the straight **chord** connecting the two endpoints of the normalized
frontier: the no-trade point `(rho=0, c_norm=0)` and the full-rebalance point
`(rho=1, c_norm=1)`. This chord is the curve's average cost-per-completion
over the whole rebalance.

For each swept point `(rho_j, c_norm_j)`, compute its **perpendicular
distance** to that chord — the length of the perpendicular leg dropped from
the point onto the chord. With the chord running corner to corner of the unit
square, the signed perpendicular distance reduces to:

```
d_j = ( rho_j - c_norm_j ) / sqrt(2)
```

The frontier bows below the chord (the curve completes faster than its
average rate early on), so `d_j > 0` along it.

The knee is the point of **maximum perpendicular distance**:

```
knee = argmax_j  d_j
```

Equivalently — and this is the same point — the knee is where the curve's
local slope equals the chord's slope: before the knee the curve completes
faster than its all-rebalance average, after it slower. The first-derivative
and longest-perpendicular views agree by construction; the agreement is the
method's internal consistency check.

This detector needs no second derivative, so it does not amplify the small
sweep-sampling wiggles a curvature detector would. It is global — it weighs
the whole curve at once — not a three-point local estimate.

### S8.3 Sharpness — read directly from the leg length

The maximum perpendicular distance `d_max` **is itself** the knee's sharpness
measure. No separate threshold constant is invented for it. A sharp knee bows
far from the chord (`d_max` large); a gently bowed curve with no real corner
hugs the chord (`d_max` small). On the normalized unit square `d_max` is
bounded by `1/sqrt(2) ≈ 0.707` (a perfect right-angle corner).

If `d_max` is below a sharpness cutoff the frontier has no real corner: the
engine reports **gradual tradeoff — no sharp optimum**, suppresses the knee
marker, defaults to the stated convention `rho = 0.90`, and prompts the PM to
choose. An honest "no clear optimum" beats a confident marker on a flat curve.

The cutoff is the one detector parameter (`KNEE_MIN_LEG`, S13). Unlike a raw
curvature threshold it is an interpretable geometric quantity — "how far the
curve must bulge from the straight line to count as a knee" — but its exact
value still depends on what real client frontiers look like, so it is
calibrated on real sweeps, not derived. State it as calibratable; do not
treat it as a fixed specification constant.

### S8.4 Coupling to cost data

The knee location depends on `k_i`. Higher commissions → the cost axis
weights turnover more → the knee moves left → the engine rebalances less.
This is automatic and correct: the knee is computed on the cost-inclusive
frontier. Cheap-to-trade book → knee moves right. Per-asset `k_i` (tier 1+,
partial or full column) genuinely moves the knee; a flat tier-0 default does
not (S8.1).

---

## S9. Final solution and snap delta

1. Determine target `rho`: the knee (S8.2), or the PM's knob override.
2. Solve `F(y)` at the corresponding `theta` (interpolate `theta` from the
   ladder, then one exact solve).
3. Apply post-processing (§5.4) — long-only clamp, then `$40` snap — to this
   **final** solution only.
4. Recompute `J_dev`, `rho`, and `cost` of the snapped solution. Report the
   delta from the pre-snap marker (`snap_delta`, S10) so the PM sees the real
   delivered numbers, not the idealized frontier point.

The knee is found on the true (un-snapped) frontier; snap is a bounded
final rounding (sub-$40 per trade) whose effect is measured and shown.

---

## S10. Diagnostics additions

Extends §7. New fields in `diagnostics`:

- `total_cost` — `cost(x)` in dollars of the delivered solution.
- `cost_bps` — `total_cost / V0 * 1e4`.
- `objective_deviation` — `J_dev` of the delivered solution.
- `objective_cost` — `theta * cost(x) / V0`.
- `objective_value` — their sum (existing field, now split-reported).
- `rho_achieved` — completion fraction of the delivered solution.
- `knee_rho` — knee location, or `null` if no sharp knee.
- `knee_leg` — `d_max`, the longest perpendicular leg length from S8.2;
  the sharpness measure.
- `frontier` — the swept `(rho, cost, J_dev)` points, for plotting.
- `snap_delta` — `{ d_rho, d_cost }`, pre-snap vs post-snap (S9).
- `cost_tier` — 0 / 1 / 2 (S2.1).
- `cost_default_used` — bool: any asset used a defaulted rate.

Per-trade additions to each `trades` entry: `trade_cost` (`k_i|x_i| +
q_i x_i^2`).

Mode note: in pinned-cash modes (A withdraw, D deposit) `sum|x|` is largely
fixed by the cash target, so cost steers *which* assets trade, not the total
volume. In B/C/E cost also shrinks total turnover. Surface this so a PM is
not surprised cost "does nothing" to turnover in mode A — emit
`cost_affects = "selection"` (A, D) or `"selection+volume"` (B, C, E).

Field names are contract (§7) — identical across `rebalancer.py`,
`engine.js`, `widget_ui.js`.

---

## S11. UI changes (`widget_ui.js`)

1. **Frontier plot** — cost (commission dollars) vs `rho`, the concave curve,
   drawn from `diagnostics.frontier`.
2. **Knee marker** — dropped at `knee_rho`, labelled with the `rho` value.
   When `knee_rho` is null, show the "gradual tradeoff — no sharp optimum"
   state instead.
3. **`rho` knob** — `[0,1]`, defaulting to the knee. PM drags to override;
   table, exposure bars, diagnostics re-render. Reuses `solveCached` /
   `warmup`. The knob is an audit/override control, not the deciding act —
   the engine has already asserted the knee.
4. **Cost axis label** — always commission dollars. At tier 0 label it
   "assumes N bps blanket rate — upload a cost column for per-asset
   accuracy." A defaulted number must never look like measured data.

---

## S12. Reference test cases

Companion to the worked example in `WIKI.md` §2.7 — these are the
cost-model-specific cases.

- **Tier 0 equivalence.** Workbook with no cost column → flat
  `k_i = COST_DEFAULT` → trade vector matches the pre-change `TAU` engine to
  within snap tolerance.
- **Buy/sell split exactness.** For any tier-1 solve with `theta > 0`, assert
  `u_i * v_i = 0` for every `i` (no simultaneous buy and sell).
- **Knee monotonicity.** Raise every `k_i` uniformly → `knee_rho` does not
  increase (knee moves left or stays, S8.4).
- **Flat-default knee invariance.** Tier 0 with two different `COST_DEFAULT`
  values → identical `knee_rho` and `knee_leg` (S8.1, S8.3).
- **Py/JS agreement.** Sweep the same workbook with both engines → frontier
  points agree to 4 significant figures, `knee_rho` agrees to `0.01`.

---

## S13. Constants

| Constant | Value | Role |
|---|---|---|
| `COST_DEFAULT` | `8e-4` (8 bps) | tier-0 flat commission rate; tier-0 fallback |
| blank-cell percentile | `90th` | default for blank cells in a present column |
| `theta` ladder span | `10^-3 .. 10^3`, 25 pts | initial sweep grid |
| `theta` refine gap | `0.05` in `rho` | adaptive bisection trigger |
| sweep solve cap | `60` | max solves per sweep |
| `KNEE_MIN_LEG` | calibrate on real frontiers | min `d_max` (perpendicular leg) for a sharp knee (S8.3) |

`COST_DEFAULT`, percentile, ladder params and solve cap are specification
(changing them changes output). `KNEE_MIN_LEG` is the one tunable — it sets
detector sensitivity, a UX threshold, not a result-changing value. It is a
perpendicular distance on the normalized unit square (bounded by
`1/sqrt(2) ≈ 0.707`); no default is fixed here — calibrate on real sweeps by
finding the gap between frontiers that visibly have a knee and those that do
not.

The old `TAU` (`1e-4`) is removed from the objective.

---

## S14. Out of scope

- Fixed per-ticket fees `k_fixed * 1[x_i != 0]` — a 0/1 indicator,
  non-convex (cardinality), breaks the active-set solver. The `$40` snap
  (§5.4) is the pragmatic stand-in. Heuristic upgrade if needed: solve
  convex, snap, re-solve with snapped assets pinned at 0.
- Concave cost (economies of scale) — non-convex, out.
- Hard cost budget (`cost <= B`) — a second linear inequality, a real solver
  extension. If ever needed without the extension: binary-search `theta`
  until `cost <= B` (cost is monotone decreasing in `theta`), reusing the
  sweep machinery.
- Dollarizing the deviation objective via a rebalance horizon — an
  alternative route to the same optimum, not taken here. The knee method
  reaches the economically optimal point from curve geometry and needs no
  horizon input; `lambda_d` stays a dimensionless priority (§2.3).

---

## S15. Implementation order

1. `rebalancer.py`: column matching for cost columns (S3), defaults (S3.1).
2. `rebalancer.py`: buy/sell split — rewrite `rebalance` QP assembly to the
   `2n` system (S4, S5). Solver (`_solve_as`) untouched.
3. `rebalancer.py`: sweep function — `theta` ladder, adaptive refine, snap
   off during sweep (S7).
4. `rebalancer.py`: knee detector (S8), final-solve + snap-delta (S9),
   diagnostics fields (S10).
5. Port 2–4 to `engine.js` line-for-line.
6. Cross-check both engines (S12).
7. `widget_ui.js`: frontier plot, knee marker, `rho` knob (S11).
8. Update `WIKI.md` §8 (constants), §10.1 (cost term now implemented).

Steps 2–5 are the math change — `rebalancer.py` and `engine.js` must land
together.
