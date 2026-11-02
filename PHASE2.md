# CryoSand — Phase 2 Plan

**Period:** November–December 2026
**Starting point:** branch `first-report-baseline` · `python -m pytest` → 57 passed, 1 skipped
**Phase 2 goal:** answer both project questions with a fair comparison and sourced
inputs, and build an interactive app for experimenting with every part of the trade.

---

## 0. How to use this file

**Start of every session**
1. Read PROJECT.md §3 (conventions) and this file.
2. Run `python -m pytest`. Do not start work on a red baseline.
3. Pick the first unchecked box in the build order (§5).

**While working**
- Where a test is listed, write it before the code.
- When a task needs a number that can't be sourced, stop and flag it (§3). Don't guess.

**End of every session**
- Tick finished boxes and add one line to the status log (§11).
- Do not commit, push or open PRs. The owner does that.

**Which document wins**
- This file: what to build in Phase 2, and in what order.
- PROJECT.md: conventions and physics (after WP0 brings it up to date).
- CURRENT_STATE.md: a 1 Nov snapshot. Planning input only.

---

## 1. Goals

| # | Goal | Delivered by |
|---|---|---|
| G1 | Objective R: closed-tank crossover under every closure, with a fair comparison | WP2, WP4 |
| G2 | Objective E: first crossover map (propellant × radius × environment × fill × vent pressure) | WP7 |
| G3 | Interactive app with live visualizers | WP3, WP8 |
| G4 | Sourced parameter register with ranges; t* reported as a band | WP5 |
| G5 | Validation: bracket vs measured self-pressurization; MHTB | WP6 |
| G6 | Report 2 figures and numbers, generated, never typed | WP9 |

---

## 2. Decisions (settled — reopen only with the owner)

### D1 — One cooler-mass model: Strobridge

Every code path that produces a reported number prices the cooler with
`cooler_mass_strobridge` (Plachta & Kittel 2003):

```
m_cool = 0.2 * Q_lift**0.7 * ((T_h - T_c) / T_c)**1.45      [kg, Q_lift in W]
```

A constant kg/W survives only as an **analytic limit**, inside the closed-form
functions in `l5_closure.py` and the test that checks them.

Correct wording for docs and docstrings:
- End-point optimality (report 1, App. B.2) holds **with** Strobridge, because
  linear + concave is still concave.
- Only the closed form `t* = h_fg * mu_act` needs a constant kg/W. With
  Strobridge it is approximate.

### D2 — Three designs, and a fair closed-tank trade

Report 1 credited the closed tank's stored energy `E_store` to the passive design
only. The ZBO design still lifted the whole leak. But a cooler on a closed tank
only needs to lift enough heat that the tank reaches `P_vent` right at mission end.

| Design | Cooler lift | Meaning |
|---|---|---|
| `PASSIVE` | 0 | no cooler; vents after the hold time |
| `PARTIAL` | `x_k = clip(Q_leak - E_store/t, 0, Q_leak)` | "just enough" cooler; tank reaches P_vent at mission end |
| `ZBO` | `Q_leak` | holds storage pressure indefinitely |

Why only these three. Closed tank, fixed insulation, lift `x` in `[0, Q_leak]`:

```
m(x) = m_ins + m_cool(x) + m_pow(x) + m_rad(x) + max((Q_leak - x)*t - E_store, 0) / h_fg
```

- For `x < x_k` the boil-off term is linear and the cooler term is concave, so the
  sum is concave → minimum at `x = 0` or `x = x_k`.
- For `x >= x_k` nothing boils off and mass rises with `x` → minimum at `x = x_k`.
- So the optimum is PASSIVE or PARTIAL. ZBO is optimal only when `E_store = 0`
  (vented tank), where PARTIAL = ZBO.
- Ties go to PASSIVE (no hardware).

With constant kg/W at fixed insulation:

```
t*_closed = max(t*_vented, t_hold)     # fair definition
t*_closed = t*_vented + t_hold         # report-1 definition
```

**Both definitions are kept and reported:**
- `active=(Design.ZBO,)` — report-1 definition. Kept for continuity and for
  store-forever architectures (depots with no end date).
- `active=(Design.PARTIAL, Design.ZBO)` — fair definition. For missions with an end date.

In a vented tank the two must agree exactly.

### D3 — The conduction anchor is computed, never hard-coded

`f ≈ 0.10` (LH2) came from one case. It depends on fluid, heat leak, interface
area and time. Compute it per case with `m_surf_conduction` + `brentq`. No literal
0.10 or 0.07 in any code path.

### D4 — Report figures wait for Gate B

Sweeps and the app can run on provisional numbers. Figures and tables for
report 2 cannot (§5).

### D5 — No silent defaults

Trade functions take `closure` and the `active` design set as explicit arguments
with no default. At startup the app shows all closures side by side, so none is
privileged.

---

## 3. Rules

Carried from PROJECT.md §3 and §12:
- SI base units in all physics. Convert only for display.
- CoolProp is the only fluid-property source.
- Pure physics functions; state in frozen dataclasses.
- Test first wherever an analytic limit exists.
- Parameter values live only in `data/parameters.yaml`. Report numbers are generated macros.
- If a term is negligible, report the number. Don't drop it silently.
- `python -m pytest` must pass before work is handed back.

New for Phase 2:
- **Never invent a value or a data point.** An unsourced parameter goes in the
  register with a range, `status: provisional` and the note `SOURCE NEEDED`.
  Measured data enters only through files the owner adds to `data/validation/`.
  No values or digitized points from memory.
- **The app holds no physics.** It builds inputs, calls `model` / `trade`, and
  converts units for display.
- **No tool credits anywhere.** No generated-by notes, co-author trailers or tool
  names in code, comments, docs, commit messages or the report.
- **Git is owner-only.** No commits, pushes, PRs or branch changes.

---

## 4. Reference numbers (regression anchors)

Baseline: LH2, LEO 400 km, r = 2 m, L = 1.5 r, 90 % fill, 1 atm storage, vent 3 bar.

| Quantity | Value |
|---|---|
| Outer surface | T_s = 259.0 K (q_abs = 219.5 W/m², F_E = 0.331) |
| One sub-blanket | R = 3.158 K/W (10 layers, 8.33 mm, 0.2 kg/m²) |
| Passive optimum at 57 d | n = 7, Q_leak = 11.99 W (struts 1.19 W), 254.7 kg |
| ZBO optimum | n = 6, Q_leak = 13.79 W, W_in = 2,524 W, A_rad = 6.14 m², 255.4 kg |
| Specific power φ | 183 (LH2), 11.6 (LO2), 8.4 (LCH4) |
| Radiator flux | 413.4 W/m² at 300 K, ε = 0.9 |
| Vented t* | 57.3 d (LH2), 1.6 d (LO2), 2.7 d (LCH4) |
| Homogeneous E_store, hold | 2.11×10⁸ J, 204 d at 11.99 W |
| Report-1 closed t* (LH2) | 60.4 d (f = 1e-3), 84.6 d (anchor), 377 d (homogeneous) |
| With controller mass 1.4× | 91.2 d (insulation fixed, report 1 App. A.5) |

Any change that moves these must be intended and explained in the status log.

---

## 5. Build order and gates

| Order | Work package | Gate |
|---|---|---|
| 1 | WP0 Housekeeping | |
| 2 | WP1 One cooler-mass model | |
| 3 | WP2 Fair closed-tank trade | |
| 4 | WP3 App v1: single-case explorer | **Gate A** — owner starts experimenting |
| 5 | WP4 Closed-tank completeness (+ app controls) | |
| 6 | WP5 Register, ranges, model options (+ app controls) | |
| 7 | WP6 Validation | **Gate B** — report figures allowed |
| 8 | WP7 Sweeps and crossover map | |
| 9 | WP8 App v2: map, sensitivity, validation, register pages | |
| 10 | WP9 Report 2 outputs | **Gate C** — clean, committed tree |
| 11 | WP10 2-D wall port | optional before report 2 |

- **Gate A:** WP0–WP3 ticked, all tests green, app runs at baseline and shows vented t* = 57.3 d.
- **Gate B:** every parameter that drives t* is sourced, or is provisional with a
  range and included in the t* band. Bracket check and MHTB checks run and recorded.
- **Gate C:** report outputs written by `make_figures.py` from a committed tree (`git_dirty = false`).

---

## 6. Work packages

### WP0 — Housekeeping

- [x] `requirements.txt` pinned to the versions now installed: numpy, scipy,
      CoolProp, pyyaml, pandas, streamlit, plotly, pytest. Add venv setup steps to the README.
      (matplotlib kept too — used by make_figures; plotly not yet installed, left unpinned until it is.)
- [x] Bring PROJECT.md up to date (L3 two-zone form was already current):
  - §4 L3: replace the c_p surface bound with the two-zone form below. Note that
    the c_p form failed the ordering test and survives only as a regression test.
  - §4 L0: the sphere-to-Earth view factor F_E is in. Attitude dependence is still spring.
  - §4 L5 and §6: add the `Design` enum and the three designs (D2).
  - §7: `test_bounds_ordered` checks stored energy at every f, plus pressurization rate.
  - §11: current status.
- [x] `scripts/run_case.py`: replace the `f = 0.01` default with the per-case
      conduction anchor (D3). Keep `--f` as an override.
      (shared helper `model.conduction_anchor_fraction`; override flag is `--surf-fraction`.)

Two-zone form for PROJECT.md §4 L3:

```
m_w     = m_surf + m_v
V_w     = V - (m_l - m_surf) / rho_l
du_w/dt = Q_net / m_w
P       = PropsSI('P', 'D', m_w / V_w, 'U', u_w, fluid)
```

**Done when:** 57 passed, 1 skipped; `run_case.py --days 57` still prints vented t* = 57.3 d.

### WP1 — One cooler-mass model (D1)

- [x] Test first: `test_run_steady_matches_trade` — for the baseline ZBO design,
      `model.run_steady` cooler, power, radiator and total masses equal
      `trade.zbo_breakdown` to 1e-9 relative.
- [x] `model.run_steady`: price the cooler with `cooler_mass_strobridge`
      (`l5_closure.mass_breakdown` now takes a pre-priced `m_cooler`).
- [x] `Cryocooler.specific_mass` renamed to `s_cool_const`, docstring "analytic
      limit only". Read only by `mu_active`, `t_crossover_vented` and
      `test_crossover_matches_linearity_result` (rewritten to the constant-kg/W limit).
- [x] Fix docstrings in `l5_closure.py` to the D1 wording.

**Done when:** one code path prices the cooler; all reference numbers unchanged.

### WP2 — Fair closed-tank trade (D2) — research-critical

Tasks:
- [ ] `core/types.py`: `class Design(Enum): PASSIVE, PARTIAL, ZBO`.
- [ ] `layers/l5_closure.py` (pure functions):
  - `lift_partial(Q_leak, E_store, t) -> float` [W]
  - `mass_at_lift(...)` — mass breakdown for any lift, so tests can brute-force it.
- [ ] `trade.py`:
  - `design_mass(c, n, t, closure, design, E_store)`
  - `best_design(c, t, closure, design, E_store) -> (mass, n, Q_lift)`
  - `t_crossover(c, closure, E_store, active)` — `active` is a tuple from
    {PARTIAL, ZBO}; no default.
  - `partial_breakdown(...)`: insulation, cooler, power, radiator, boil-off, Q_lift, hold time.
- [ ] `run_case.py`: print t* under both definitions for every closure.
- [ ] Add the D2 derivation to METHODOLOGY.md as App. B.5.
- [ ] Stretch: `cross_check(c, n, t, sized_under, run_under, f)` — size the PARTIAL
      cooler under one closure, run it under another; return time of first vent and
      propellant lost. Shows what the mixing assumption costs a real design.

Tests (write first):

| Test | Assertion |
|---|---|
| `test_partial_equals_zbo_when_vented` | E_store = 0 ⇒ PARTIAL mass = ZBO mass (1e-12) |
| `test_partial_equals_passive_within_hold` | t ≤ E_store/Q_leak ⇒ lift = 0 and PARTIAL mass = PASSIVE mass |
| `test_closed_optimum_is_zero_or_hinge` | brute force 10⁴ lifts in [0, Q_leak] with Strobridge at several (t, E_store): minimum within one grid step of 0 or x_k |
| `test_closed_crossover_linear_limit` | constant s_cool, fixed n: t* = max(t*_vented, t_hold) to 1e-6 relative |
| `test_partial_never_heavier_than_zbo` | for all t and all E_store ≥ 0 |
| `test_vented_definitions_agree` | VENTED: both definitions give 57.3 d at baseline |
| `test_crossover_sign_changes_once` | on a 200-point log grid of t, (best passive − best active) changes sign at most once |

Sanity expectations (not tests — if the code disagrees by a lot, look for a bug first):
- LH2, conduction anchor: hold ≈ 25 d < 57 d, so the fair t* should land near 57 d, not 84.6 d.
- LH2, homogeneous: the fair t* should fall below 377 d. The fixed-insulation
  estimate is max(57, 204) = 204 d; re-optimizing insulation can push it higher.

**Done when:** all tests pass; `run_case.py` prints both definitions.

### WP3 — App v1: experiment with one case

Run with `streamlit run cryosand/app/streamlit_app.py`.

```
cryosand/app/
  streamlit_app.py            # home page: Case Explorer
  pages/
    1_Inside_the_Tank.py
    2_Insulation_and_Cooler.py
  components/
    sidebar.py                # all inputs → one frozen AppCase
    compute.py                # cached calls into model/trade; no physics
    figures.py                # plotly builders: results in, figure out
    units.py                  # SI ↔ display (days, bar, kg/day, US units)
```

**Sidebar** (shared by all pages, held in `st.session_state`)

| Input | Control | Range |
|---|---|---|
| Propellant | select | LH2, LO2, LCH4 |
| Environment | select | presets in the register |
| Tank radius | slider | 0.5–5 m |
| Fill fraction | slider | 10–95 % |
| Mission duration | log select-slider | 1 d – 10 yr |
| Closures | multiselect | vented, homogeneous, two-zone (all three at startup) |
| Layer fraction f | log select-slider + "use conduction anchor" checkbox | 1e-4 – 1 |
| Designs | multiselect | passive, partial, ZBO (all three at startup) |
| Crossover definition | radio | report-1 (vs ZBO) / fair (vs best active) |
| Register overrides | expander of number inputs, bounded by register ranges | "modified" badge when changed |
| Units | toggle | SI / US |

**Page — Case Explorer**
- Number cards: outer-surface T; heat leak split into blanket and strut paths;
  boil-off (kg/day); cooler power; radiator area; hold time per closure; t* per
  closure; winning design at the chosen duration.
- Mass vs duration: best mass of each design vs duration (log x), one panel per
  closure, crossover markers, vertical line at the chosen duration.
- Mass breakdown: stacked bars per design at the chosen duration (insulation,
  cooler, power, radiator, boil-off, and tank-wall penalty once WP4 adds it).
- Heat-flow diagram (plotly Sankey): leak → boil-off or cooler lift → cooler
  electrical input → radiator, every link labeled in W. Show the environment's
  absorbed power as a separate card; it is ~1,600× the leak and would swamp the diagram.
- Fidelity banner: closures shown, number of provisional parameters in use, model
  limits in force (constant k_eff, full sun, 4 K radiator sink, 1-D wall).

**Page — Inside the Tank**
- Pressure history to P_vent for homogeneous and two-zone at the chosen f, with a
  shaded band for f from 1e-4 to 1, the vent line and hold-time markers.
- Tank cross-section (schematic, plotly shapes, not draggable): liquid level from
  fill, warm layer drawn to scale from m_surf, ullage colored by pressure at a time slider.
- E_store vs f (log x): homogeneous as a horizontal line, conduction anchor marked,
  dry-out region shaded.
- If the liquid-full guard trips, show its message instead of the plots.
- Stretch: cross-check panel (WP2).

**Page — Insulation and Cooler**
- Heat leak vs number of sub-blankets, stacked blanket and strut paths.
- Mass vs sub-blankets for each design at the chosen duration, optimum marked;
  warning when the optimum sits at n = 1 or n_max.
- Specific power vs cold-end temperature: Carnot and the chosen η, propellant markers.
- Radiator-temperature explorer: active mass vs T_rad (with T_h = T_rad), showing
  the interior optimum. The trade still holds T_rad fixed; this page only shows the
  optimum exists.

**Behaviour**
- One-case update under 1 s on a laptop. Cache CoolProp-backed results with
  `st.cache_data` keyed on the frozen case. Vectorize `Q_leak` over n with numpy if needed.
- One color per design and one line style per closure, defined once in
  `figures.py` and used everywhere. Propellant colors match report 1 Fig. 5.1.
- Every plot downloadable (plotly toolbar). A "Download case" button saves inputs,
  outputs, overrides and provenance as JSON.

Tests:
- `test_app_pages_run`: `streamlit.testing.v1.AppTest` runs each page at startup
  inputs without an exception.
- `test_app_baseline_numbers`: the home page at baseline shows vented t* = 57.3 d.
- `test_figures_build`: each builder returns a figure with the expected traces for
  the baseline results.

**Gate A** after this package.

### WP4 — Closed-tank completeness

Each item gets its control in the app when it lands.

- [ ] **Vent pressure** as a TradeCase field and sweep variable (slider from just
      above storage pressure to a register maximum). E_store recomputed per P_vent.
- [ ] **Boil-off after first vent** uses h_fg at P_vent for closed closures.
      Test `test_vent_hfg_reduces_to_vented`: P_vent = P_store and E_store = 0 ⇒
      identical to vented.
- [ ] **Storage margin** `storage_margin` in [0, 1): usable energy =
      (1 − margin)·E_store, applied to PASSIVE and PARTIAL alike.
      Tests: margin 0 reproduces WP2; margin → 1 reproduces vented.
- [ ] **Tank-wall pressure penalty**, pure function, thin wall with a minimum gauge:
      cylinder t = SF·P·r/σ_allow, hemispheres half that.
      Penalty = m_wall(P_max of the design) − m_wall(P_base), floored at 0, where
      P_base is the pressure the tank is built for anyway (register). In orbit the
      wall carries the absolute internal pressure. Applies to closed PASSIVE and
      PARTIAL; vented and ZBO stay at storage pressure. App toggle, off until its
      register entries are sourced; the fidelity banner says when it is off.
      Tests: zero when P_vent ≤ P_base; thickness matches a hand calculation.
- [ ] **Liquid-full in sweeps** is a flagged result (`status = "liquid_full"`),
      never a crash and never dropped. Test at 95 % fill.

### WP5 — Register, ranges and model options

- [ ] Register schema per entry: `value, low, high, unit, source, status (sourced |
      provisional | assumption), note`. Test `test_register_schema`: every entry has
      every field; every non-sourced entry has low ≤ value ≤ high.
- [ ] New entries, all `provisional` + `SOURCE NEEDED` until the owner fills them:
  - cooler controller mass factor (Plachta & Kittel use 1.4 × cooler mass)
  - cooler reject ΔT above T_rad
  - radiator α, ε and view factor to Earth
  - strut scaling exponent
  - foam (SOFI) conductivity and thickness (needed for MHTB)
  - tank wall σ_allow, density, safety factor, minimum gauge, P_base
- [ ] **Controller mass** m_ctrl = factor · m_cool. Tests: factor 0 reproduces the
      baseline; factor 1.4 with insulation fixed (n = 7 passive, n = 6 ZBO) gives 91.2 d.
- [ ] **Radiator environment load**: q_net = ε_r·σ·(T_rad⁴ − T_sink⁴) − q_env,rad,
      with q_env,rad from the L0 functions using the radiator's α, ε and view factor.
      Option toggle. Tests: q_env,rad = 0 reproduces 413.4 W/m²; q_net ≤ 0 raises a
      clear error.
- [ ] **Strut scaling**: G = G_ref·(r / r_ref)^k. Test: r = r_ref gives 0.005 W/K.
- [ ] **Optional foam layer** in L1, in series with the blankets. Test: zero
      thickness reproduces the baseline.
- [ ] **Sensitivity**: `sensitivity(case, closure, active)` returns t* at each
      non-sourced parameter's low and high (one at a time), plus a seeded Monte Carlo
      band (uniform over ranges, N = 500). Saved to `runs/`.
- [ ] **Plachta & Kittel preset** (`pk2003`): fixed 243 K environment temperature,
      273 K reject, controller factor 1.4, k_eff matched to their MLI heat flux.
      Lockheed coefficients are not entered until sourced. Output: our t* for LH2,
      LO2 and LCH4 under their assumptions beside their values (about two months,
      one week and two weeks).

### WP6 — Validation

- [ ] **Data format** (owner supplies): `data/validation/<name>.yaml` — fluid, tank
      volume and shape, fill, heat input (W) or flux, heating configuration, initial
      pressure, citation with figure or table number — plus `<name>.csv` with
      columns `time_s, P_Pa`. Code never invents points.
- [ ] **Bracket check** `scripts/bracket_check.py`: for each dataset, run homogeneous
      and two-zone (f band) pressure histories at the experiment's own heat input;
      plot the data over the bounds; output a table — inside, above the fast bound,
      or below the slow bound, with times. This is a finding, not a pass/fail test.
- [ ] Unit test `test_bracket_machinery`: a synthetic dataset built from the
      homogeneous bound lies exactly on the slow bound.
- [ ] Fast-bound note: Aydelott's top-heating runs exceeded the surface-evaporation
      limit. If data rises faster than the two-zone bound at its smallest f, report
      it. Never tune f to hide it.
- [ ] **MHTB**, two separate checks, each named in its test docstring:
  - (a) L1: predicted heat flux vs measured 0.085 W/m² at a 164 K warm boundary and
    0.22 W/m² at 305 K (Martin & Hastings 2001). Expect constant k_eff to miss the
    trend: measured flux rises 2.6× while ΔT rises 2.0×. That gap is the evidence
    for the k_eff(T) upgrade.
  - (b) L3: boil-off from the measured heat input vs measured boil-off, within
    10 % → un-skip `test_boiloff_published`.

### WP7 — Sweeps and the crossover map

- [ ] `trade.sweep(grid, closures, active_sets)` → tidy pandas DataFrame, one row
      per (propellant, environment, radius, fill, P_vent, closure, f, active set):
      t*, optimum n per design, Q_leak, Q_lift, hold time, and flags
      (`liquid_full`, `optimizer_at_edge`, `radiator_cannot_reject`). Saved through
      `run_case.py --save` with provenance.
- [ ] Cache E_store per (fluid, fill, P_vent, closure, f).
- [ ] `optimizer_at_edge` flag when the optimum n is 1 or n_max. Test.
- [ ] Map script: t* over radius × environment for each propellant and closure,
      plus a difference map (homogeneous minus two-zone at the anchor) showing where
      mixing matters.
- [ ] Target: the default grid finishes in under 5 minutes on a laptop. Add
      multiprocessing only if needed.

### WP8 — App v2

- [ ] `3_Crossover_Map.py`: heatmap of t* over two chosen axes (radius,
      environment, fill, P_vent); closure and definition selectors; difference map;
      hover shows designs and flags. Loads the latest saved sweep, or runs a small
      grid live with a progress bar.
- [ ] `4_Sensitivity.py`: tornado chart of t*; Monte Carlo histogram with the band.
- [ ] `5_Validation.py`: bracket plot per dataset; MHTB comparison; inside/outside table.
- [ ] `6_Register.py`: every parameter with value, range, source and status; filter
      to provisional; session overrides highlighted.
- [ ] Stretch: compare tray — pin up to 3 cases and overlay their mass-vs-duration curves.
- [ ] Tests: `AppTest` for each new page.

### WP9 — Report 2 outputs (requires Gate B)

- [ ] `make_figures.py` regenerates: crossover map (fig3), hold-time band (fig4),
      sensitivity table, bracket plot, fair-trade table (t* by propellant × closure ×
      definition), and `results_macros.tex`.
- [ ] Every figure: lettering at least 3 mm tall at printed size; SI with US units
      (secondary axis or labels).
- [ ] Gate C check: refuse to write report outputs when `git_dirty` is true or a
      t*-driving parameter is provisional without a range, unless `--draft` (adds a
      draft watermark). Test.
- [ ] Rebuild: `latexmk -pdf`, then `python scripts/build_docx.py`.

### WP10 — 2-D wall port (optional before report 2)

- [ ] Port the ThermoSand solver into `l2_wall.py`: axisymmetric, harmonic-mean
      faces, radiative boundary via h_r = εσ(T² + T∞²)(T + T∞), iterated.
- [ ] Tests: `test_robin_limits`; `test_2d_matches_1d_slab` (thin wall, 1e-6
      relative); existing `test_harmonic_mean`.
- [ ] App page: wall temperature field. Report: the wall resistance number, stated
      as negligible.

---

## 7. Out of scope (Spring 2027)

- Two-phase CFD dataset and surrogate; PINN surrogate.
- Lockheed MLI correlation coefficients (until sourced).
- Attitude dependence, eclipses, orbit-transient heating.
- Thermodynamic vent system.
- Optimizing T_rad inside the trade (WP3 only visualizes the optimum).
- Draggable or free-form UI elements.

## 8. Owner-only tasks

- Source register values and ranges (WP5 list).
- Digitize Aydelott, Lin and Van Dresar pressure histories and MHTB boil-off into
  `data/validation/`, with citations.
- Choose P_base (the pressure the tank is built for anyway).
- All git operations.
- Report: own-words rewrite of the 25-Sep passages, advisor block in the LaTeX,
  reference checks (Zumbo and Hastings initials, Daigle pages, ZBOT authors).

## 9. Suggested schedule

| Week of | Work |
|---|---|
| 2 Nov | WP0, WP1, WP2 |
| 9 Nov | WP3 (Gate A), WP4 |
| 16 Nov | WP5; owner sourcing and digitizing in parallel |
| 23 Nov | WP6 (Gate B), WP7 |
| 30 Nov | WP8 |
| 7 Dec | WP9 (Gate C); WP10 if time |

## 10. Phase 2 is done when

- [ ] All tests green, including the WP2 fair-trade tests and `test_boiloff_published`.
- [ ] t* for LH2, LO2 and LCH4 under every closure and both definitions, with a sensitivity band.
- [ ] Crossover map and difference map generated from a committed tree.
- [ ] Bracket check run on every supplied dataset, with the result recorded either way.
- [ ] The app runs all pages and reproduces the reference numbers.
- [ ] PROJECT.md matches the code.

## 11. Status log

| Date | Session | Done | Notes |
|---|---|---|---|
| 2026-11-02 | planning | PHASE2.md written | Baseline 57 passed / 1 skipped |
| 2026-11-02 | WP0 | requirements pinned + README venv; PROJECT.md L0/§4 L5/§6 Design/§7/§11; run_case default f → conduction anchor (shared helper in model.py) | 57 passed / 1 skipped; run_case --days 57 vented t*=57.3 d, surface anchor f=0.0997 → 84.6 d (matches §4); plotly still to install before WP3 |
| 2026-11-02 | WP1 | run_steady prices cooler via Strobridge (mass_breakdown takes pre-priced m_cooler); Cryocooler.specific_mass → s_cool_const (analytic-limit only); l5_closure D1 docstrings; added test_run_steady_matches_trade; linearity test rewritten to constant-kg/W limit | 58 passed / 1 skipped; reference numbers unchanged (254.7 / 255.4 kg, t*=57.3 / 84.6 / 377.3 d) |
