# CryoSand — PROJECT.md

Single source of truth for this repository. Read this first in any new session.

**Status:** MS project, accepted. SJSU Aerospace, Fall 2026 – Spring 2027.

---

## 1. What this is

An open, interactive tool for the in-space cryogenic storage trade: **insulate,
or refrigerate?** The user picks a mission and propellant, adjusts the insulation
layup and cryocooler capacity, and watches boil-off rate, ullage pressure
history, radiator area and total system mass respond live.

Two results come out of it.

**The engineering result.** A published map of where the passive/active
crossover sits, across propellant, tank scale, thermal environment and mission
duration. No such map exists in the open literature; the tools that could
produce one (GFSSP, Thermal Desktop, BoilFAST, CryoTank, Simcenter Amesim) all
sit behind licences.

**The research result.** Whether the universal well-mixed-ullage assumption
actually displaces that crossover, or merely alters the vent schedule. If the
crossover is insensitive, the assumption is demonstrably safe for architecture
decisions. If it moves, the fidelity required has been quantified.

---

## 2. The research method — read this before touching L3

The second result is delivered by a **bounding argument**, not by CFD. This is
the single most important design decision in the project and it governs how L3
is implemented.

Self-pressurisation is analytically bracketed:

- **Homogeneous limit** — the entire fluid mass warms together, internal energy
  rises at the heat rate, pressure follows from the equation of state. This is
  the *slowest* pressure rise physically available.
- **Surface limit** — all incoming energy is absorbed by a thin layer at the
  liquid–vapour interface, and saturation pressure tracks that layer alone. This
  is the *fastest*.

Reality is provably between them. So: propagate **both** bounds through the full
trade and observe whether the crossover moves. That answers the research
question with no simulation at all.

The spring CFD campaign then locates reality *within* bounds already
established — refining the answer rather than producing it. If two-phase
cryogenic CFD proves intractable, the project forfeits refinement, not its
result.

**Implication for the code:** both closures are first-class, both are always
available, and every trade function must be runnable under either. Never let the
homogeneous closure become a default that the surface closure is bolted onto.

---

## 3. Conventions — non-negotiable

- **SI base units everywhere, internally.** Metres, kilograms, seconds, kelvin,
  watts, pascals, joules. No kJ. No bar. No centimetres. Convert only at the
  display layer, never inside physics code.
- **CoolProp is the single source of fluid properties.** Never hand-roll an
  equation of state or hard-code a latent heat. `PropsSI` returns SI already,
  which is why the rule above costs nothing.
- **Every physics function is pure.** Inputs in, number out, no hidden state.
  State lives in frozen dataclasses passed explicitly.
- **Naming:** `T_*` temperatures (K), `Q_*` heat rates (W), `E_*` energies (J),
  `m_*` masses (kg), `mdot_*` mass flows (kg/s), `A_*` areas (m²), `R_*` thermal
  resistances (K/W). A variable named `Q` holding joules is the bug that costs a
  week.
- **Anything checkable against an analytic limit gets a test that checks it.**
  Write the test first where practical.

**Environment:** Python 3.11+, numpy, scipy, CoolProp, streamlit, pytest.
Run tests with `pytest -q` before every commit.

---

## 4. Physics specification

### L0 — Orbital environment

```
q_abs = alpha * (G_sun + G_albedo) + eps * G_IR      [W/m^2]
```

`alpha` = solar absorptivity, `eps` = IR emissivity. Radiates to a ~4 K sink, so
the outbound term is `eps * sigma * T_surf**4`.

**Current scope:** one absorbed-flux number per mission preset. View factors and
attitude dependence are a spring upgrade — do not build them now.

### L1 — Insulation stack

Series thermal resistances from outer surface to tank wall.

```
R_blanket = t / (k_eff * A)                          [K/W]
R_strut   = L / (k_bar * A_strut)
            k_bar = (1/(T_h - T_c)) * integral(k dT)
Q_total   = (T_hot - T_cold) / sum(R_series)  +  parallel strut path
```

**Struts and penetrations are a parallel path, not a series element.** They
frequently dominate a well-insulated tank. Model them as a conductance added to
the blanket conductance.

MLI `k_eff` is strongly temperature-dependent (roughly 1e-4 to 3e-5 W/m·K in
hard vacuum). Start from published effective values per blanket. Residual gas
conduction is absorbed into `k_eff` at this fidelity — state that explicitly in
docstrings rather than silently omitting it. The modified Lockheed correlation
is a spring upgrade; do **not** hard-code its coefficients until sourced.

### L2 — Tank wall

Fourier conduction, `div(k grad T) = 0`, ported from the existing ThermoSand
solver. Harmonic-mean face conductivity at material interfaces:

```
k_face = 2 * k_P * k_N / (k_P + k_N)
```

Two half-cell resistances in series. The arithmetic mean can overstate flux by
orders of magnitude, and mesh refinement does not repair it.

**Expect wall resistance to be negligible** — aluminium at 2 mm is nearly a short
circuit. Include it, verify it, then report that it is small. Knowing a term is
negligible because you computed it is not the same as assuming it.

### L3 — Ullage and liquid

Three closures. All three are required.

**Vented (constant pressure):**
```
mdot_boil = Q_net / h_fg
```

**Closed, homogeneous bound** — constant volume, constant mass:
```
du/dt = Q_net / m_total
rho   = m_total / V                       (constant)
P     = PropsSI('P', 'D', rho, 'U', u, fluid)
```
Exact for a well-mixed tank. No correlation required.

**Closed, surface bound** — all energy into a *warm zone*: an interface layer
of liquid mass `m_surf` plus the ullage vapour. The remaining liquid is inert at
its initial state and occupies a fixed volume.
```
m_w   = m_surf + m_v                     V_w = V - (m_liquid - m_surf) / rho_l
du_w/dt = Q_net / m_w
P     = PropsSI('P', 'D', m_w / V_w, 'U', u_w, fluid)   # = P_sat(T_int) while two-phase
```
`m_surf / m_liquid` is a swept parameter, and the sweep is a result, not a
nuisance. As `m_surf -> m_liquid` the warm zone is the whole tank, so this
collapses onto the homogeneous bound **exactly, by construction**.

*History (Sept 2026):* the textbook form `dT_int/dt = Q_net / (m_surf * cp)`,
`P = P_sat(T_int)` was implemented first. `test_bounds_ordered` caught it
storing more energy than the homogeneous bound at 90 % fill (it ignores the
vapour, and the condensation of a compressed ullage). It is kept in `l3_ullage`
only as a regression record. A thin layer may evaporate completely before
`P_max` (superheated warm zone); re-condensation onto the cold liquid is
neglected there, which keeps the result a valid fast bound.

Closed-tank energies raise `LiquidFullError` if the absorbing zone becomes
liquid-full before `P_max` (e.g. 95 % fill at 3 bar for all three propellants).

### L4 — Pressure control and active cooling

```
COP_carnot = T_cold / (T_hot - T_cold)
W_in       = Q_lift / (eta_carnot * COP_carnot)
Q_reject   = Q_lift + W_in
A_rad      = Q_reject / (eps * sigma * (T_rad**4 - T_sink**4))
```

| T_cold | ideal W/W | % Carnot | real W/W |
|---|---|---|---|
| 20 K (LH2) | 14 | 5–10 % | 140–280 |
| 90 K (LO2) | 2.3 | 15–25 % | 9–15 |

**This asymmetry is the headline engineering result.** Zero boil-off is nearly
free for oxygen and punishing for hydrogen.

`T_rad` has an interior optimum: raising it shrinks radiator area but is bounded
by the reject temperature the cooler tolerates. Hold it fixed for now; note the
optimum exists.

**Current scope:** thermodynamic vent system deferred. Passive vent plus
cryocooler only.

### L5 — System closure

```
m_insulation = areal_density * A_tank * n_blankets
m_cooler     = Q_lift * specific_mass_cooler      [kg/W lift, T-dependent]
m_power      = W_in   * specific_mass_power       [kg/W]
m_radiator   = A_rad  * areal_density_radiator    [kg/m^2]
m_lost       = mdot_boil * t_mission

m_penalty = m_insulation + m_cooler + m_power + m_radiator + m_lost
```

**The crossover is the minimum of `m_penalty` over the design space.** Passive
wins where `m_lost` is cheap; active wins where it is not.

---

## 5. Repository layout

```
cryosand/
  core/
    constants.py        # sigma, G_sun, etc. Named, sourced, SI.
    properties.py       # CoolProp wrapper. Cached. Single entry point.
    network.py          # series/parallel resistance primitives
  layers/
    l0_environment.py
    l1_insulation.py
    l2_wall.py          # ported from ThermoSand
    l3_ullage.py        # vented + homogeneous bound + surface bound
    l4_cooling.py
    l5_closure.py
  model.py              # assembles layers; steady + transient drivers
  trade.py              # parameter sweeps, crossover map, both closures
  app/
    streamlit_app.py
  tests/
    test_limits.py      # analytic limits — write these first
    test_layers.py
    test_validation.py  # against published data
  data/
    materials.yaml      # k(T), areal densities, with sources
    missions.yaml       # presets
```

---

## 6. Core data structures

```python
@dataclass(frozen=True)
class Propellant:
    name: str                 # CoolProp name: 'Hydrogen', 'Oxygen', 'Methane'

@dataclass(frozen=True)
class TankGeometry:
    radius: float             # m
    length: float             # m, cylindrical section
    wall_thickness: float     # m
    wall_material: str
    fill_fraction: float      # liquid volume / total volume
    # derived: volume, surface_area, wetted_area

@dataclass(frozen=True)
class InsulationStack:
    blankets: list            # each: (thickness_m, k_eff, areal_density)
    strut_conductance: float  # W/K, parallel path

@dataclass(frozen=True)
class Environment:
    q_absorbed: float         # W/m^2
    T_sink: float             # K

@dataclass(frozen=True)
class Cryocooler:
    Q_lift: float             # W; 0 means fully passive
    T_cold: float             # K
    T_reject: float           # K
    eta_carnot: float
    specific_mass: float      # kg per W of lift

@dataclass(frozen=True)
class Mission:
    name: str
    duration: float           # s
    environment: Environment

class UllageClosure(Enum):
    VENTED = auto()
    HOMOGENEOUS = auto()      # slow bound
    SURFACE = auto()          # fast bound
```

Every trade function takes `UllageClosure` as an explicit argument. There is no
default.

---

## 7. Validation targets — write as tests

| test | assertion |
|---|---|
| `test_zero_leak` | `Q_net == 0` ⇒ `mdot_boil == 0` exactly |
| `test_harmonic_mean` | `k_face(2, 6) == 3.0` |
| `test_series_chain` | network matches analytic 1-D resistance sum to 1e-10 |
| `test_robin_limits` | `h -> inf` gives Dirichlet at `T_inf`; `h -> 0` gives insulated |
| `test_zbo_balance` | `Q_lift == Q_leak` ⇒ `mdot_boil == 0` to 1e-6 |
| `test_bounds_ordered` | surface-bound `dP/dt` ≥ homogeneous `dP/dt`, at every condition |
| `test_bounds_collapse` | `m_surf -> m_liquid` ⇒ surface bound = homogeneous bound (to 1e-10) |
| `test_energy_balance` | transient energy closes to better than 0.5 % |
| `test_boiloff_published` | matches a published demonstration-tank result within 10 % |

**`test_bounds_ordered` protects the entire research claim.** If it ever fails,
the bounding argument is broken and every downstream conclusion is void. Write
it before any L3 physics code.

---

## 8. Build order

Each step ends with tests green before the next begins.

1. `constants.py` + `properties.py`. Test: CoolProp returns sane `h_fg` for LH2.
2. `network.py`. Test: series and parallel resistance against hand calculations.
3. `l1_insulation.py`. Test: analytic chain.
4. `l5_closure.py` + `l4_cooling.py`. Test: ZBO balance, Carnot limits.
5. `l3_ullage.py`, vented mode only. Test: zero leak, published boil-off.
6. `model.py` — steady assembly. **First end-to-end run.**
7. `l3_ullage.py`, both closed-tank bounds. Test: ordering and collapse.
8. `l2_wall.py` port from ThermoSand. Test: harmonic mean, Robin limits.
9. `trade.py` — sweeps under both closures. **First crossover map.**
10. `app/streamlit_app.py`. Sliders and plots only; nothing draggable.

Steps 1–6 are the minimum viable tool. **Step 7 and step 9 together carry the
research claim.** Step 10 makes it usable.

---

## 9. Reference values

| quantity | value |
|---|---|
| sigma (Stefan–Boltzmann) | 5.670374419e-8 W/m²K⁴ |
| Solar constant | 1361 W/m² |
| Earth IR | ~237 W/m² |
| Earth albedo | ~0.3 |
| Deep-space sink | ~4 K |
| LH2: T_sat, rho_l, h_fg | 20.3 K, 70.8 kg/m³, 446 kJ/kg |
| LO2: T_sat, rho_l, h_fg | 90.2 K, 1141 kg/m³, 213 kJ/kg |
| LCH4: T_sat, rho_l, h_fg | 111.7 K, 422 kg/m³, 511 kJ/kg |
| Good MLI k_eff | ~3e-5 W/m·K |
| 1 W into LH2 | 0.19 kg/day boil-off |
| Published depot ZBO need | ~80–100 W @ 80 K (LO2); ~100–120 W @ 20 K (LH2) |

For orientation and sanity checks only. **Production code queries CoolProp.**

---

## 10. Known traps

- **Mixing J and kJ.** Most likely source of a factor-1000 error.
- **Struts in series instead of parallel.** They are a parallel conductance.
- **Treating MLI `k_eff` as constant.** Varies strongly with boundary
  temperatures. Acceptable now; must be flagged in the report.
- **Forgetting the radiator rejects `Q_lift + W_in`, not `Q_lift`.** For LH2,
  `W_in` exceeds `Q_lift` by two orders of magnitude.
- **Letting one ullage closure become the default.** Both bounds are first-class.
- **Choosing `m_surf` arbitrarily.** It is swept; the sweep is a result.
- **`T_rad` treated as fixed without noting the optimum exists.**

---

## 11. Status

**Complete (Summer 2026):** 2-D FDM thermal solver — Jacobi, red-black SOR at
optimal omega, explicit transient under CFL, Dirichlet/Neumann/Robin boundary
conditions, harmonic-mean variable-k stencils, Streamlit interface. Becomes
Layer L2. PINN module in progress.

**Complete (Sept 2026, first written report):** build steps 1–7. 54 tests pass
under CoolProp: analytic limits, bound ordering at every `m_surf` for three
propellants and four fills, exact collapse, energy closure, liquid-full guard.
Vented crossover map (step 9, vented closure) and closed-tank crossover at the
LEO baseline. All report numbers regenerated from CoolProp.

**Next:**
- Closed-tank crossover sweeps across radius, environment and fill (step 9 under
  all three closures).
- One cooler-mass model: `model.run_steady` uses linear `specific_mass`, while
  `trade.py` and the figures use the Strobridge correlation.
- Load `data/parameters.yaml` in `make_figures.py` instead of duplicating values;
  source every `provisional` entry.
- Step 8 (2-D L2 port), `test_boiloff_published` (MHTB), step 10 (Streamlit).

**Deferred to spring:** two-phase CFD dataset and trained surrogate; modified
Lockheed MLI correlation; view factors and attitude dependence; thermodynamic
vent system.

**December target:** complete tool released, first crossover map produced, both
ullage bounds implemented and their effect on the crossover measured.

---

## 12. Working agreement

When helping on this repository:

- Write the test before the physics wherever an analytic limit exists.
- Never introduce a unit conversion inside a physics function.
- Never hard-code a fluid property — call CoolProp.
- Cite the source in a comment for any empirical constant or correlation.
- If a term turns out negligible, say so with the number, don't silently drop it.
- Flag rather than guess when a value is unsourced; unsourced constants become
  report liabilities later.
