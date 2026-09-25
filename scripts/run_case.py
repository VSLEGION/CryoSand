"""Run one CryoSand case, or a sweep, and save the result.

Examples
--------
    python scripts/run_case.py                                  # baseline: LH2, LEO, r = 2 m
    python scripts/run_case.py --fluid LO2 --env deep --days 90
    python scripts/run_case.py --surf-fraction 0.001            # thinner interface layer
    python scripts/run_case.py --sweep radius=0.75:5:10 --save radius_lh2
    python scripts/run_case.py --sweep fill=0.3:0.9:7 --fluid LO2 --save fill_lo2
    python scripts/run_case.py --register my_register.yaml     # alternative parameter set

Every value not given on the command line comes from data/parameters.yaml.
With --save NAME the run is written to runs/<timestamp>_<NAME>.json (and .csv
for a sweep), together with what is needed to reproduce it: the command, the
git commit, the register file and its hash, and the CoolProp version.

Units: SI inside the model. Days, bar and kg/day appear only in the printout.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cryosand import scenario as S  # noqa: E402
from cryosand.core import constants as C  # noqa: E402
from cryosand.core import properties as props  # noqa: E402
from cryosand.core.types import UllageClosure  # noqa: E402
from cryosand.layers import l4_cooling as L4  # noqa: E402
from cryosand.model import ClosedTankSpec, LiquidFullError, stored_energy  # noqa: E402
from cryosand import trade as T  # noqa: E402

DAY = C.SECONDS_PER_DAY
ALIASES = {"lh2": "Hydrogen", "hydrogen": "Hydrogen", "lo2": "Oxygen", "oxygen": "Oxygen",
           "lox": "Oxygen", "lch4": "Methane", "methane": "Methane"}
SWEEPABLE = ("radius", "fill", "days", "p_vent", "surf_fraction")


# ------------------------------------------------------------------ core ---
def run(p: S.Params, fluid: str, env: str, radius: float, fill: float, days: float,
        p_vent: float, surf_fraction: float) -> dict:
    """Evaluate one case. Returns a flat dict of SI results (plus *_days keys)."""
    e = S.environment(p, env)
    sat = S.saturated(fluid, p)
    g = S.tank(p, radius, fill)
    c = S.trade_case(p, fluid, e["T_s"], radius)
    t = days * DAY
    out = dict(fluid=fluid, env=env, radius=radius, fill=fill, days=days, p_vent=p_vent,
               surf_fraction=surf_fraction, T_sat=sat["T_sat"], h_fg=sat["h_fg"],
               q_abs=e["q_abs"], T_surface=e["T_s"], V_tank=g.volume, A_tank=g.surface_area,
               m_liquid=sat["rho_l"] * fill * g.volume)

    # vented passive: best insulation for THIS duration
    m_pas, n_pas = T.best_passive(c, t, UllageClosure.VENTED, 0.0)
    Q_pas = T.Q_leak(c, n_pas)
    out.update(n_passive=n_pas, layers_passive=n_pas * p.unit_layers, Q_leak_passive=Q_pas,
               mdot_boil=Q_pas / sat["h_fg"], m_passive=m_pas,
               **{f"passive_{k}": v for k, v in
                  T.passive_breakdown(c, n_pas, t, UllageClosure.VENTED, 0.0).items()})

    # zero boil-off: best insulation (independent of duration)
    m_zbo, n_zbo = T.best_zbo(c)
    Q_zbo = T.Q_leak(c, n_zbo)
    W = L4.work_input(Q_zbo, c.T_cold, c.T_reject, c.eta_carnot)
    out.update(n_zbo=n_zbo, layers_zbo=n_zbo * p.unit_layers, Q_leak_zbo=Q_zbo, W_in_zbo=W,
               A_rad_zbo=L4.radiator_area(Q_zbo + W, p.eps_rad, p.T_rad, C.T_DEEP_SPACE),
               m_zbo=m_zbo, **{f"zbo_{k}": v for k, v in T.zbo_breakdown(c, n_zbo).items()})
    out["winner_vented"] = "passive" if m_pas < m_zbo else "ZBO"
    out["t_star_vented_days"] = T.t_crossover(c, UllageClosure.VENTED, 0.0) / DAY

    # closed tank, both bounds. Hold times use the passive-optimum leak.
    for tag, closure, f in (("hom", UllageClosure.HOMOGENEOUS, 1.0),
                            ("surf", UllageClosure.SURFACE, surf_fraction)):
        try:
            E = stored_energy(closure, fluid, g, p.P_storage, ClosedTankSpec(p_vent, f))
        except LiquidFullError:
            out.update({f"E_store_{tag}": math.nan, f"t_hold_{tag}_days": math.nan,
                        f"t_star_{tag}_days": math.nan, f"liquid_full_{tag}": True})
            continue
        out.update({f"E_store_{tag}": E, f"t_hold_{tag}_days": E / Q_pas / DAY,
                    f"t_star_{tag}_days": T.t_crossover(c, closure, E) / DAY,
                    f"liquid_full_{tag}": False})
    return out


# ------------------------------------------------------------- printing ---
def show(r: dict) -> None:
    kg = lambda x: f"{x:10.1f} kg"  # noqa: E731
    print(f"\n=== {r['fluid']} | {r['env']} | r = {r['radius']:.2f} m | fill {r['fill']:.0%} | "
          f"{r['days']:g} days | vent {r['p_vent'] / 1e5:.2f} bar ===")
    print(f"Propellant   T_sat {r['T_sat']:.2f} K   h_fg {r['h_fg'] / 1e3:.1f} kJ/kg   "
          f"liquid {r['m_liquid']:,.0f} kg")
    print(f"Environment  q_abs {r['q_abs']:.1f} W/m2   outer surface {r['T_surface']:.1f} K")
    print(f"Tank         V {r['V_tank']:.1f} m3   A {r['A_tank']:.1f} m2")
    print("\n                    PASSIVE (vented)      ZERO BOIL-OFF")
    print(f"MLI layers         {r['layers_passive']:>10d}          {r['layers_zbo']:>10d}")
    print(f"Heat leak          {r['Q_leak_passive']:10.2f} W        {r['Q_leak_zbo']:10.2f} W")
    for k in ("insulation", "cooler", "power", "radiator", "lost"):
        print(f"  m_{k:<14}{kg(r['passive_m_' + k])}      {kg(r['zbo_m_' + k])}")
    print(f"  TOTAL           {kg(r['m_passive'])}      {kg(r['m_zbo'])}")
    print(f"Boil-off (passive) {r['mdot_boil'] * DAY:.2f} kg/day   "
          f"ZBO cooler input {r['W_in_zbo']:.0f} W, radiator {r['A_rad_zbo']:.1f} m2")
    print(f"\nLighter at {r['days']:g} days (vented): {r['winner_vented']}")
    print("\nCrossover t* [days]   vented {:.1f}   surface bound (m_surf/m_liq = {:g}) {}   "
          "homogeneous bound {}".format(
              r["t_star_vented_days"], r["surf_fraction"], _d(r, "t_star_surf_days"),
              _d(r, "t_star_hom_days")))
    print("Hold time to vent [days]           surface {}   homogeneous {}".format(
        _d(r, "t_hold_surf_days"), _d(r, "t_hold_hom_days")))


def _d(r, k):
    return "LIQUID-FULL" if r.get(k.replace("t_star", "liquid_full").replace("t_hold", "liquid_full")
                                  .replace("_days", "")) else f"{r[k]:.1f}"


SWEEP_COLS = ["layers_passive", "Q_leak_passive", "m_passive", "m_zbo", "winner_vented",
              "t_star_vented_days", "t_hold_surf_days", "t_star_surf_days",
              "t_hold_hom_days", "t_star_hom_days"]


def show_sweep(param: str, rows: list[dict]) -> None:
    cols = [param] + SWEEP_COLS
    print("  ".join(f"{c:>16}" for c in cols))
    for r in rows:
        print("  ".join(f"{r[c]:>16.4g}" if isinstance(r[c], float) else f"{str(r[c]):>16}"
                        for c in cols))


# --------------------------------------------------------- provenance ---
def provenance(p: S.Params) -> dict:
    def git(*a):
        try:
            return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True,
                                  check=True).stdout.strip()
        except Exception:  # noqa: BLE001
            return None
    import CoolProp
    reg = Path(p.source)
    return dict(timestamp=dt.datetime.now().isoformat(timespec="seconds"),
                command="python " + " ".join(sys.argv), git_commit=git("rev-parse", "HEAD"),
                git_dirty=bool(git("status", "--porcelain")),
                register=str(reg.relative_to(ROOT)) if reg.is_relative_to(ROOT) else str(reg),
                register_sha256=hashlib.sha256(reg.read_bytes()).hexdigest(),
                coolprop_version=CoolProp.__version__, python=sys.version.split()[0])


def save(name: str, meta: dict, rows: list[dict]) -> Path:
    runs = ROOT / "runs"
    runs.mkdir(exist_ok=True)
    stem = runs / f"{dt.datetime.now():%Y%m%d-%H%M%S}_{name}"
    stem.with_suffix(".json").write_text(json.dumps(dict(meta=meta, results=rows), indent=2,
                                                    default=float))
    if len(rows) > 1:
        with stem.with_suffix(".csv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    return stem


# ------------------------------------------------------------------ CLI ---
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--register", help="parameter register YAML (default data/parameters.yaml)")
    ap.add_argument("--fluid", default="LH2", help="LH2 | LO2 | LCH4 (default LH2)")
    ap.add_argument("--env", default="leo", help="environment key in the register (default leo)")
    ap.add_argument("--radius", type=float, help="tank radius [m] (default: register)")
    ap.add_argument("--fill", type=float, help="fill fraction 0-1 (default: register)")
    ap.add_argument("--days", type=float, default=30.0, help="mission duration [days] (default 30)")
    ap.add_argument("--p-vent", type=float, help="vent pressure [Pa] (default: register)")
    ap.add_argument("--surf-fraction", type=float, default=1e-2,
                    help="surface bound m_surf/m_liquid in (0, 1] (default 0.01)")
    ap.add_argument("--sweep", help=f"PARAM=start:stop:n, PARAM in {SWEEPABLE}; "
                                    "log-spaced if prefixed log:, e.g. surf_fraction=log:1e-4:1:9")
    ap.add_argument("--save", metavar="NAME", help="save to runs/<timestamp>_NAME.json/.csv")
    a = ap.parse_args(argv)

    if not props.coolprop_available():
        sys.exit("CoolProp is required: pip install -r requirements.txt")
    p = S.load_params(a.register)
    fluid = ALIASES.get(a.fluid.lower())
    if fluid is None:
        sys.exit(f"unknown fluid {a.fluid!r}; use LH2, LO2 or LCH4")
    base = dict(fluid=fluid, env=a.env,
                radius=p.radius if a.radius is None else a.radius,
                fill=p.fill if a.fill is None else a.fill, days=a.days,
                p_vent=p.P_vent if a.p_vent is None else a.p_vent,
                surf_fraction=a.surf_fraction)

    if a.sweep:
        param, spec = a.sweep.split("=", 1)
        if param not in SWEEPABLE:
            sys.exit(f"cannot sweep {param!r}; choose from {SWEEPABLE}")
        log = spec.startswith("log:")
        lo, hi, n = spec.removeprefix("log:").split(":")
        vals = (np.logspace(math.log10(float(lo)), math.log10(float(hi)), int(n)) if log
                else np.linspace(float(lo), float(hi), int(n)))
        rows = [run(p, **{**base, param: float(v)}) for v in vals]
        print(f"\nSweep of {param} | {fluid} | {a.env} | register {Path(p.source).name}")
        show_sweep(param, rows)
    else:
        rows = [run(p, **base)]
        show(rows[0])

    if a.save:
        stem = save(a.save, dict(provenance(p), inputs=base, sweep=a.sweep), rows)
        print(f"\nSaved {stem.name}.json" + (" and .csv" if len(rows) > 1 else "")
              + " in runs/. Record what you learned in runs/LOG.md.")


if __name__ == "__main__":
    main()
