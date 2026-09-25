# Running CryoSand yourself

This guide takes you from a fresh laptop to running your own experiments,
recording them, and turning them into figures and text for the report. It
starts with copy-and-paste steps and gets more technical as it goes.

---

## 0. The whole pipeline in one picture

```
data/parameters.yaml          <- the ONLY place parameter values live
        |
        v
cryosand/scenario.py          <- turns the register into model inputs
        |
        +--> scripts/run_case.py     you, exploring: one case or a sweep -> runs/*.json, *.csv
        |
        +--> scripts/make_figures.py the report: figures + results_macros.tex
                     |
                     v
             report/CryoSand_FirstReport.tex  (every number is a macro, never typed by hand)
```

The physics lives in `cryosand/layers/` (L0–L5), gets assembled in
`cryosand/model.py`, and gets optimised in `cryosand/trade.py`. Section 7
covers that code.

---

## 1. Set up on your laptop (once)

You need Python 3.11 or newer. Check with `python3 --version`.

**macOS / Linux**
```bash
git clone https://github.com/VSLEGION/CryoSand.git
cd CryoSand
git checkout first-report-baseline         # until this branch is merged
python3 -m venv .venv                      # a private Python just for this project
source .venv/bin/activate                  # do this in every new terminal
pip install -r requirements.txt
```

**Windows (PowerShell)**
```powershell
git clone https://github.com/VSLEGION/CryoSand.git
cd CryoSand
git checkout first-report-baseline
py -3.11 -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

*Why a venv?* It keeps this project's library versions separate from
everything else on your machine, so results don't change when some other
project upgrades numpy.

To rebuild the PDF you also need a LaTeX distribution: MacTeX on macOS,
MiKTeX on Windows, or `texlive-latex-extra` plus `latexmk` on Linux.

---

## 2. Level 1: prove the install works

```bash
python -m pytest
```

You should see `57 passed, 1 skipped`. The skip is `test_boiloff_published`,
which is waiting for MHTB data.

**What the tests are telling you.** Each test is a claim about the physics
that has to stay true. These are the ones to know by name, because the report
leans on them:

| Test | The claim it defends |
|---|---|
| `test_bounds_ordered` | The surface bound always pressurises at least as fast as the fully mixed bound. **If this fails, the research argument is void.** |
| `test_bounds_collapse` | Growing the surface layer to all the liquid gives exactly the fully mixed answer. |
| `test_cp_form_violates_ordering_at_high_fill` | Records why the textbook surface model was dropped. |
| `test_zbo_balance` | A cooler that lifts exactly the heat leak gives zero boil-off. |
| `test_crossover_matches_linearity_result` | Brute-force optimisation reproduces t* = h_fg · μ_active. |
| `test_run_case_reproduces_report_baseline` | Your run tool and the report agree. |

To watch one test run in detail: `python -m pytest -v -s tests/test_coolprop.py -k collapse`.

---

## 3. Level 2: your first case

```bash
python scripts/run_case.py --days 57
```

This runs the baseline: liquid hydrogen, LEO at 400 km, a 2 m radius tank, 90 %
full, venting at 3 bar. Every value not given on the command line comes from
`data/parameters.yaml`. Read the output top to bottom:

1. **Propellant / Environment / Tank.** These are inputs, resolved: CoolProp's
   saturation temperature and latent heat, the absorbed solar and Earth flux,
   and the outer-surface temperature it produces.
2. **The two-column table.** The trade itself. For each column the code has
   already picked the *best number of MLI layers*:
   - **Passive:** more insulation versus more boil-off. The best layer count
     depends on mission duration.
   - **Zero boil-off (ZBO):** insulation versus cooler, power system and
     radiator. The best layer count doesn't depend on duration.

   `m_lost` is the propellant boiled away over the mission.
3. **"Lighter at N days".** Which column has the smaller total.
4. **Crossover t\*.** The duration at which the two totals are equal, under
   three ullage assumptions: vented, the surface bound at your chosen layer
   size, and the fully mixed (homogeneous) bound.
5. **Hold time.** How long a *closed* tank lasts before it reaches the vent
   pressure. It's computed at the heat leak of the passive design for your
   `--days`. Change `--days` and the best insulation changes, and so does the
   hold time.

**Try it:** run `--days 50` and `--days 65`. Watch the winner flip, and check
that the flip sits at the vented t\* printed in the last block.

### Every option

| Option | Meaning | Default |
|---|---|---|
| `--fluid` | `LH2`, `LO2`, `LCH4` | `LH2` |
| `--env` | `leo`, `deep`, `mars` (add your own in the register) | `leo` |
| `--radius` | tank radius [m] | register (2.0) |
| `--fill` | liquid fill fraction, 0–1 | register (0.9) |
| `--days` | mission duration [days] | 30 |
| `--p-vent` | vent pressure [Pa], SI on purpose | register (3e5) |
| `--surf-fraction` | surface-bound layer size, m_surf/m_liquid, in (0, 1] | 0.01 |
| `--register` | use a different parameter file | `data/parameters.yaml` |
| `--sweep` | vary one parameter (see Level 3) | – |
| `--save NAME` | write the result to `runs/` | not saved |

If a closed-tank result says `LIQUID-FULL`, the liquid swelled to fill the
tank before reaching the vent pressure. That's a real physical limit, not a
bug. Try `--fill 0.95`.

---

## 4. Level 3: ask questions with sweeps

A sweep varies one input and prints one row per value:

```bash
python scripts/run_case.py --sweep radius=0.75:5:10                 # 10 evenly spaced radii
python scripts/run_case.py --sweep surf_fraction=log:1e-4:1:9       # 9 log-spaced layer sizes
python scripts/run_case.py --sweep fill=0.3:0.9:7 --fluid LO2
python scripts/run_case.py --sweep days=1:400:20 --save duration_lh2
```

The format is `PARAM=start:stop:count`; add a `log:` prefix for log spacing.
You can sweep `radius`, `fill`, `days`, `p_vent` and `surf_fraction`.

**Always predict before you run.** Write one sentence in `runs/LOG.md` saying
what you expect and why, *then* run it. When the model surprises you, either
you've learned some physics or you've found a bug. Both are worth a
paragraph in the report.

---

## 5. Level 4: change the assumptions

Most parameters in `data/parameters.yaml` are marked `status: provisional`.
Sourcing them is your register task, and changing them is how you run a
sensitivity study.

**Don't edit the main register to experiment.** Copy it:

```bash
cp data/parameters.yaml runs/register_eta_low.yaml
# edit runs/register_eta_low.yaml: e.g. Hydrogen eta_carnot 0.075 -> 0.05
python scripts/run_case.py --register runs/register_eta_low.yaml --save eta_low
```

When you've *sourced* a value (found it in a paper), edit
`data/parameters.yaml`: change the value, set `status: sourced`, add the
`source:` citation. Then regenerate the report (Level 5). Every figure and
number moves together, because they all read the same file.

A YAML detail: PyYAML reads `3.0e5` as *text*, not a number. The loader
converts every value with `float()`, which handles it, but `300000.0` or
`3.0e+5` is the unambiguous way to write it.

---

## 6. Level 5: capture results and write about them

### Every saved run is reproducible

`--save NAME` writes `runs/<timestamp>_NAME.json`, plus a `.csv` for a sweep.
The JSON `meta` block records:

- the exact command you typed
- the git commit, and whether you had uncommitted changes (`git_dirty`)
- the register file and its SHA-256 hash, so you can prove which parameters
  produced a number
- the CoolProp and Python versions

If `git_dirty` is `true`, commit first and re-run before you cite the
number. A result you can't reproduce is a result you can't defend.

### Keep a log

`runs/LOG.md` is your lab notebook. One entry per question: what you
expected, the command, what happened, and what it means. Most of your
results section will come straight out of this file.

### Plot a sweep yourself

```python
import csv, matplotlib.pyplot as plt
rows = list(csv.DictReader(open("runs/20261001-101500_duration_lh2.csv")))
x = [float(r["days"]) for r in rows]
plt.plot(x, [float(r["m_passive"]) for r in rows], label="passive (vented)")
plt.plot(x, [float(r["m_zbo"]) for r in rows], "--", label="ZBO")
plt.xlabel("Mission duration [days]"); plt.ylabel("Mass penalty [kg]"); plt.legend()
plt.savefig("runs/duration_lh2.png", dpi=200)
```

### Regenerate the report

```bash
python scripts/make_figures.py         # figures, results.json, results_macros.tex
cd report && latexmk -pdf CryoSand_FirstReport.tex
```

Numbers in the report text are LaTeX macros (`\tsH`, `\tchH`, …) defined in
the auto-generated `results_macros.tex`. Never type a result into the `.tex`
by hand. Add a macro in `make_figures.py` instead, so the text can never
drift from the code.

---

## 7. Level 6: changing the physics

| File | What lives there |
|---|---|
| `cryosand/core/constants.py` | σ, solar constant, Earth IR, each with its source |
| `cryosand/core/properties.py` | the only door to CoolProp |
| `cryosand/core/network.py` | series/parallel resistances, harmonic-mean face k |
| `cryosand/layers/l0_environment.py` | absorbed flux, surface temperature |
| `cryosand/layers/l1_insulation.py` | MLI in series, struts in parallel |
| `cryosand/layers/l2_wall.py` | wall resistance (the 2-D port goes here) |
| `cryosand/layers/l3_ullage.py` | vented, homogeneous, two-zone surface closures |
| `cryosand/layers/l4_cooling.py` | Carnot, specific power, radiator, Strobridge mass |
| `cryosand/layers/l5_closure.py` | mass bookkeeping, crossover closed form |
| `cryosand/model.py` | assembles layers; closed-tank energies and pressure histories |
| `cryosand/trade.py` | optimises insulation, finds the crossover |
| `cryosand/scenario.py` | register → model inputs |

The working rules (full list in `PROJECT.md` §3 and §12):

1. **Write the test first** when an analytic limit exists. The ordering test
   caught the surface-closure error *because* it was written before the physics.
2. **SI inside, always.** Convert to days, bar or kJ only when printing.
3. **Never hard-code a fluid property.** Call `cryosand.core.properties`.
4. **Cite any empirical number** in a comment where it's used, and in the register.
5. Run `python -m pytest` before every commit.

---

## 8. First experiments to try

These are questions, not answers. Predict, run, then explain.

1. **The cooler asymmetry.** Compare `--fluid LH2` and `--fluid LO2` at
   `--days 30`. Where does the ZBO mass go for each? Which term dominates, and why?
2. **Does environment matter?** Sweep `radius` under `--env leo`, `deep` and
   `mars`. The report claims t\* is nearly environment-independent. Do you
   see that? Can you explain it from t\* = h_fg · μ_active?
3. **How thin is thin?** `--sweep surf_fraction=log:1e-4:1:13`. At what layer
   size does the closed-tank crossover become twice the vented one? That
   threshold is a candidate headline number.
4. **Fill fraction.** `--sweep fill=0.3:0.9:7`, then `--sweep fill=0.80:0.96:9`,
   for LH2 and LO2. How do the two bounds' hold times respond to fill, and why
   do they move in different directions? Where does `LIQUID-FULL` start, and
   which bound reaches it first?
5. **Sensitivity.** Make a register copy with Hydrogen `eta_carnot` 0.05, and
   another with 0.10 (the provisional range). How far does t\* move? Compare
   that with how far the ullage assumption moves it. Which uncertainty matters
   more? This is objective R in one comparison.
6. **Vent pressure.** Sweep `p_vent` from 1.5e5 to 6e5 Pa. What does a
   stronger (heavier) tank buy you in crossover time?

---

## 9. Troubleshooting

| Symptom | Fix |
|---|---|
| `ModuleNotFoundError: cryosand` | Run from the repo root, with the venv activated. |
| `CoolProp is required` | `pip install -r requirements.txt` inside the venv. |
| `pytest` can't find numpy but `python -m pytest` works | Your `pytest` command belongs to a different Python. Always use `python -m pytest`. |
| `could not convert string to float` when loading the register | A non-numeric value (typo, unit text) where a number belongs in the YAML. |
| `LIQUID-FULL` in a result | Physics, not a bug: lower `--fill` or raise `--p-vent`. |
| LaTeX `Undefined control sequence \tsH` | Run `make_figures.py` first; it writes the macros. |
