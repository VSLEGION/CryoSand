"""Preliminary results for the Fall 2026 first written report.

Uses CoolProp when installed. Where it is not installed, falls back to the PROJECT.md s.9 reference values, printed loudly and
recorded in results.json. Re-run on a machine with CoolProp to replace them:
    pip install CoolProp && python scripts/make_figures.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cryosand.core import constants as C  # noqa: E402
from cryosand.core import properties as props  # noqa: E402
from cryosand.core.types import UllageClosure  # noqa: E402
from cryosand.layers import l0_environment as L0  # noqa: E402
from cryosand.layers import l1_insulation as L1  # noqa: E402,F401
from cryosand.layers import l2_wall as L2  # noqa: E402
from cryosand.layers import l3_ullage as L3  # noqa: E402
from cryosand.layers import l4_cooling as L4  # noqa: E402
from cryosand.trade import TradeCase, best_passive, best_zbo, t_crossover, Q_leak  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "report" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
DAY = C.SECONDS_PER_DAY
P0 = 101325.0

# ------------------------------------------------------------------ fluid ---
FLUIDS = ["Hydrogen", "Methane", "Oxygen"]
LABEL = {"Hydrogen": "LH$_2$", "Methane": "LCH$_4$", "Oxygen": "LO$_2$"}
# PROJECT.md s.9 reference values: no-CoolProp fallback ONLY (not a property source)
REF = {"Hydrogen": dict(T_sat=20.3, rho_l=70.8, h_fg=446e3),
       "Methane": dict(T_sat=111.7, rho_l=422.0, h_fg=511e3),
       "Oxygen": dict(T_sat=90.2, rho_l=1141.0, h_fg=213e3)}
# INDICATIVE liquid cp and k near the normal boiling point, for the closed-tank
# order-of-magnitude estimate only. Approximate handbook values (Barron, Cryogenic
# Heat Transfer); replaced by CoolProp cp_l_sat / conductivity in production.
INDICATIVE = {"Hydrogen": dict(cp=9.7e3, k=0.10, M=2.016e-3),
              "Methane": dict(cp=3.5e3, k=0.18, M=16.04e-3),
              "Oxygen": dict(cp=1.7e3, k=0.15, M=32.00e-3)}

if props.coolprop_available():
    SRC = "CoolProp"
    SAT = {f: dict(T_sat=props.T_sat(f, P0), rho_l=props.rho_l_sat(f, P0),
                   h_fg=props.h_fg(f, P0)) for f in FLUIDS}
else:
    SRC = "PROJECT.md reference values (CoolProp unavailable)"
    SAT = REF
    print("WARNING: CoolProp not installed - using PROJECT.md s.9 reference values")

# ------------------------------------------------------------- parameters ---
ETA = {"Hydrogen": 0.075, "Methane": 0.20, "Oxygen": 0.20}
T_REJECT, T_RAD, EPS_RAD, S_RAD, S_POW = 300.0, 300.0, 0.9, 5.0, 0.025
ALPHA_OUT, EPS_OUT = 0.32, 0.86
# design unit = 10-layer sub-blanket (30 layers = 25 mm => 8.33 mm per 10 layers)
BLANKET_T, K_EFF, RHO_A_BLANKET = 0.025 / 3, 3.0e-5, 10 * 0.02
G_STRUT, T_WALL, K_WALL = 5.0e-3, 0.002, 10.0
FILL, P_MAX = 0.9, 3.0e5


def env_presets():
    F_leo = L0.view_factor_sphere_to_earth(400e3, C.R_EARTH)
    presets = {
        "LEO 400 km": L0.absorbed_flux(ALPHA_OUT, EPS_OUT, C.G_SUN_1AU, C.ALBEDO_EARTH,
                                       C.G_IR_EARTH, F_leo),
        "Deep space, 1 AU": L0.absorbed_flux(ALPHA_OUT, EPS_OUT, C.G_SUN_1AU, 0.0, 0.0, 0.0),
        "Mars transfer, 1.52 AU": L0.absorbed_flux(ALPHA_OUT, EPS_OUT,
                                                   L0.solar_flux(1.52, C.G_SUN_1AU), 0.0, 0.0, 0.0),
    }
    return {k: dict(q_abs=q, T_s=L0.surface_temperature(q, EPS_OUT, C.T_DEEP_SPACE))
            for k, q in presets.items()}, F_leo


def geom(r):
    L = 1.5 * r  # cylinder section; total length 3.5 r
    V = math.pi * r**2 * L + 4 / 3 * math.pi * r**3
    A = 2 * math.pi * r * L + 4 * math.pi * r**2
    return V, A


def case(fluid, T_s, r):
    V, A = geom(r)
    return TradeCase(h_fg=SAT[fluid]["h_fg"], T_cold=SAT[fluid]["T_sat"], T_hot=T_s, A=A,
                     G_strut=G_STRUT,  # held constant across scale (provisional)
                     t_blanket=BLANKET_T, k_eff=K_EFF, rho_A_blanket=RHO_A_BLANKET,
                     R_wall=L2.R_wall_plane(T_WALL, K_WALL, A), eta_carnot=ETA[fluid],
                     T_reject=T_REJECT, s_pow=S_POW, s_rad=S_RAD,
                     q_rad=L4.radiator_flux(EPS_RAD, T_RAD, C.T_DEEP_SPACE))


# ----------------------------------------------------------------- style ---
COL = {"Hydrogen": "#2a78d6", "Methane": "#1baf7a", "Oxygen": "#eb6834"}  # validated slots 1-3
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.family": "serif", "font.size": 9, "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "lines.linewidth": 1.6, "legend.frameon": False, "savefig.dpi": 300})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png")
    plt.close(fig)


results = {"property_source": SRC}

# ======================================================= L0 environments ===
ENV, F_leo = env_presets()
results["environments"] = ENV
results["F_leo"] = F_leo

# ================================================ Fig 1: cooler asymmetry ===
Tc = np.linspace(15, 125, 300)
fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.6))
for eta, ls, lab in [(1.0, ":", "Carnot"), (0.20, "--", "20 % Carnot"), (0.075, "-", "7.5 % Carnot")]:
    phi = np.array([L4.specific_power(t, T_REJECT, eta) for t in Tc])
    ax[0].semilogy(Tc, phi, ls, color=INK2, lw=1.2, label=lab)
    ax[1].semilogy(Tc, (1 + phi) / L4.radiator_flux(EPS_RAD, T_RAD, C.T_DEEP_SPACE) * 1e2,
                   ls, color=INK2, lw=1.2)
for f in FLUIDS:
    T = SAT[f]["T_sat"]
    phi = L4.specific_power(T, T_REJECT, ETA[f])
    arad = (1 + phi) / L4.radiator_flux(EPS_RAD, T_RAD, C.T_DEEP_SPACE) * 1e2
    for a, y in [(ax[0], phi), (ax[1], arad)]:
        a.plot(T, y, "o", ms=6, color=COL[f], mec="white", mew=1.2, zorder=5)
        a.annotate(LABEL[f], (T, y), xytext=(5, 5), textcoords="offset points", color=INK)
    results.setdefault("cooler", {})[f] = dict(phi=phi, A_rad_per_100W=arad)
ax[0].set(xlabel="Cold-end temperature $T_c$ [K]", ylabel="$W_{in}/\\dot Q_{lift}$ [W/W]")
ax[1].set(xlabel="Cold-end temperature $T_c$ [K]", ylabel="Radiator area per 100 W lift [m$^2$]")
ax[0].legend(loc="upper right", fontsize=7.5)
save(fig, "fig1_cooler_asymmetry")

# ======================================== Fig 2: mass vs mission duration ===
r0 = 2.0
T_leo = ENV["LEO 400 km"]["T_s"]
days = np.logspace(0, 3.3, 200)
fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.7), sharey=False)
for a, f in zip(ax, ["Hydrogen", "Oxygen"]):
    c = case(f, T_leo, r0)
    mp = [best_passive(c, d * DAY, UllageClosure.VENTED, 0.0)[0] for d in days]
    mz, nz = best_zbo(c)
    ts = t_crossover(c, UllageClosure.VENTED, 0.0) / DAY
    a.loglog(days, mp, color=COL[f], label="Passive (vented), best $n$")
    a.axhline(mz, color=INK2, ls="--", lw=1.2, label="ZBO, best $n$")
    a.plot(ts, mz, "o", ms=7, color=COL[f], mec="white", mew=1.2, zorder=5)
    a.annotate(f"$t^*$ = {ts:.0f} d" if ts >= 10 else f"$t^*$ = {ts:.1f} d", (ts, mz),
               xytext=(8, -14), textcoords="offset points", color=INK)
    a.set_ylim(bottom=0.55 * min(mp))
    a.set(title=f"{LABEL[f]}, LEO, $r$ = {r0:.0f} m", xlabel="Mission duration [days]",
          ylabel="Mass penalty [kg]")
    V, A = geom(r0)
    results.setdefault("fig2", {})[f] = dict(t_star_days=ts, m_zbo=mz, n_zbo=nz, V=V, A=A,
                                             Q_leak_n1=Q_leak(c, 1), Q_leak_nzbo=Q_leak(c, nz))
ax[0].legend(loc="upper left", fontsize=7.5)
save(fig, "fig2_mass_vs_duration")

# ================================================ Fig 3: crossover map ===
radii = np.linspace(0.75, 5.0, 18)
fig, ax = plt.subplots(figsize=(6.5, 3.0))
styles = {"LEO 400 km": "-", "Deep space, 1 AU": "--", "Mars transfer, 1.52 AU": ":"}
cmap = {}
for f in FLUIDS:
    for e, ls in styles.items():
        ts = [t_crossover(case(f, ENV[e]["T_s"], r), UllageClosure.VENTED, 0.0) / DAY for r in radii]
        ax.semilogy(radii, ts, ls, color=COL[f])
        cmap[f"{f} | {e}"] = dict(zip([round(r, 3) for r in radii], ts))
    ax.annotate(LABEL[f], (radii[-1], cmap[f"{f} | LEO 400 km"][round(radii[-1], 3)]),
                xytext=(4, 0), textcoords="offset points", color=INK, va="center")
for e, ls in styles.items():
    ax.plot([], [], ls, color=INK2, label=e)
ax.legend(loc="upper right", fontsize=7.5, ncol=3, bbox_to_anchor=(1.0, 1.15))
ax.set(xlabel="Tank radius $r$ [m]  (cylinder length 1.5 $r$)",
       ylabel="Crossover duration $t^*$ [days]", xlim=(0.6, 5.5))
save(fig, "fig3_crossover_map")
results["crossover_map_days"] = cmap

# ===================== Fig 4: closed-tank hold-time band ===================
# With CoolProp: E_hom from model.stored_energy (exact (rho,u) path) and E_surf
# from the cp integral to T_sat(P_max). Without it: INDICATIVE cp, k and a
# Clausius-Clapeyron T_sat(P_max), and E_hom is approximated by E_surf(m_total).
from cryosand.core.types import TankGeometry  # noqa: E402
from cryosand.model import ClosedTankSpec, stored_energy  # noqa: E402

fig, ax = plt.subplots(figsize=(6.5, 2.8))
frac = np.logspace(-4, 0, 200)
hold = {}
for f in ["Hydrogen", "Oxygen"]:
    s, ind = SAT[f], INDICATIVE[f]
    c = case(f, T_leo, r0)
    _, n_p = best_passive(c, results["fig2"][f]["t_star_days"] * DAY, UllageClosure.VENTED, 0.0)
    Q = Q_leak(c, n_p)
    V, A = geom(r0)
    m = s["rho_l"] * FILL * V
    if props.coolprop_available():
        g = TankGeometry(r0, 1.5 * r0, T_WALL, "Al6061", FILL)
        Tmax = props.T_sat(f, P_MAX)
        cp, k = props.cp_l_sat(f, s["T_sat"]), props.k_l_sat(f, P0)
        E_s1 = stored_energy(UllageClosure.SURFACE, f, g, P0, ClosedTankSpec(P_MAX, 1.0))
        E_h = stored_energy(UllageClosure.HOMOGENEOUS, f, g, P0, ClosedTankSpec(P_MAX))
    else:
        Tmax = 1.0 / (1.0 / s["T_sat"] - 8.314 / (s["h_fg"] * ind["M"]) * math.log(P_MAX / P0))
        cp, k = ind["cp"], ind["k"]
        E_s1 = m * cp * (Tmax - s["T_sat"])
        E_h = E_s1
    dT = Tmax - s["T_sat"]
    t_h = frac * E_s1 / Q / DAY
    alpha = k / (s["rho_l"] * cp)
    A_int = math.pi * r0**2
    # conduction anchor: solve 0.5 rho A sqrt(pi alpha t) cp dT = Q t for t
    t_cond = (0.5 * s["rho_l"] * A_int * math.sqrt(math.pi * alpha) * cp * dT / Q) ** 2
    m_cond = L3.m_surf_conduction(s["rho_l"], alpha, A_int, t_cond)
    ax.loglog(frac, t_h, color=COL[f])
    ax.plot(1.0, E_h / Q / DAY, "s", ms=6, color=COL[f], mec="white", mew=1.2, zorder=5)
    ax.axhline(results["fig2"][f]["t_star_days"], color=COL[f], ls="--", lw=1.1)
    ax.plot(m_cond / m, t_cond / DAY, "D", ms=6, color=COL[f], mec="white", mew=1.2, zorder=5)
    ax.annotate(f"{LABEL[f]} surface bound", (frac[0], t_h[0]),
                xytext=(4, -13) if f == "Hydrogen" else (4, 10),
                textcoords="offset points", ha="left", color=INK)
    ax.annotate(f"{LABEL[f]} vented $t^*$", (1e-4, results["fig2"][f]["t_star_days"]),
                xytext=(2, 3), textcoords="offset points", color=INK2, fontsize=7.5)
    hold[f] = dict(Q_leak=Q, n=n_p, m_liquid=m, T_max=Tmax, dT=dT, alpha=alpha, cp=cp, k=k,
                   E_surf_full=E_s1, E_hom=E_h, t_hold_homog_days=E_h / Q / DAY,
                   t_hold_cond_days=t_cond / DAY, m_cond_frac=m_cond / m)
ax.plot([], [], "s", color=INK2, label="homogeneous bound")
ax.plot([], [], "D", color=INK2, label="conduction-layer anchor")
ax.legend(loc="lower right", fontsize=7.5)
ax.set(xlabel="Interface-layer mass fraction $m_{surf}/m_{total}$",
       ylabel="Hold time to 3 bar [days]", ylim=(1e-3, 3e4))
save(fig, "fig4_hold_time_band")
results["hold_time"] = hold

# ======================= closed-tank crossover shift (Table 3) ===============
shift = {}
for f in ["Hydrogen", "Oxygen"]:
    h = hold[f]
    c = case(f, T_leo, r0)
    row = {"vented": t_crossover(c, UllageClosure.VENTED, 0.0) / DAY}
    for lab, cl, E in [("s3", UllageClosure.SURFACE, 1e-3 * h["E_surf_full"]),
                       ("cond", UllageClosure.SURFACE, h["m_cond_frac"] * h["E_surf_full"]),
                       ("hom", UllageClosure.HOMOGENEOUS, h["E_hom"])]:
        row[lab] = t_crossover(c, cl, E) / DAY
    shift[f] = row
results["closed_crossover_days"] = shift

# ============================================ scalars for the report text ===
V, A = geom(r0)
R_b1 = BLANKET_T / (K_EFF * A)
results["wall"] = dict(R_wall=L2.R_wall_plane(T_WALL, K_WALL, A), R_blanket_1=R_b1,
                       ratio=L2.R_wall_plane(T_WALL, K_WALL, A) / R_b1)
results["baseline"] = dict(r=r0, V=V, A=A, T_leo=T_leo)
(OUT.parent / "results.json").write_text(json.dumps(results, indent=2, default=float))

# ---------------------------------- LaTeX macros: report numbers stay in sync
def sci(x, d=1):
    m, e = f"{x:.{d}e}".split("e")
    return f"{m}\\times10^{{{int(e)}}}"


def g3(x):
    return f"{x:.3g}" if x < 100 else f"{x:.0f}"


E, F2, CM, HT, SH = ENV, results["fig2"], results["crossover_map_days"], hold, shift
cz = {f: case(f, T_leo, r0) for f in FLUIDS}
nz = {f: best_zbo(cz[f])[1] for f in FLUIDS}
Qz = {f: Q_leak(cz[f], nz[f]) for f in FLUIDS}
macros = {
    "PropSource": "CoolProp" if props.coolprop_available() else
                  "the reference values of Table~\\ref{tab:params} (CoolProp run pending)",
    "FLeo": f"{F_leo:.3f}",
    "qLeo": f"{E['LEO 400 km']['q_abs']:.0f}", "TsLeo": f"{E['LEO 400 km']['T_s']:.0f}",
    "qDeep": f"{E['Deep space, 1 AU']['q_abs']:.0f}", "TsDeep": f"{E['Deep space, 1 AU']['T_s']:.0f}",
    "qMars": f"{E['Mars transfer, 1.52 AU']['q_abs']:.0f}",
    "TsMars": f"{E['Mars transfer, 1.52 AU']['T_s']:.0f}",
    "Vbase": f"{V:.1f}", "Abase": f"{A:.1f}",
    "Rwall": sci(results["wall"]["R_wall"]), "Rblanket": f"{R_b1:.2f}",
    "WallRatio": sci(results["wall"]["ratio"]),
    "PhiH": f"{results['cooler']['Hydrogen']['phi']:.0f}",
    "PhiO": f"{results['cooler']['Oxygen']['phi']:.1f}",
    "PhiM": f"{results['cooler']['Methane']['phi']:.1f}",
    "ArH": f"{results['cooler']['Hydrogen']['A_rad_per_100W']:.1f}",
    "ArO": f"{results['cooler']['Oxygen']['A_rad_per_100W']:.1f}",
    "ArM": f"{results['cooler']['Methane']['A_rad_per_100W']:.1f}",
    "tsH": f"{F2['Hydrogen']['t_star_days']:.0f}", "tsO": f"{F2['Oxygen']['t_star_days']:.1f}",
    "tsM": f"{CM['Methane | LEO 400 km'][2.0]:.1f}",
    "mzH": f"{F2['Hydrogen']['m_zbo']:.0f}", "mzO": f"{F2['Oxygen']['m_zbo']:.0f}",
    "nzH": f"{10 * nz['Hydrogen']}", "QzH": f"{Qz['Hydrogen']:.1f}",
    "StrutShareH": f"{100 * G_STRUT * (T_leo - SAT['Hydrogen']['T_sat']) / Qz['Hydrogen']:.0f}",
    "tsHsmall": f"{CM['Hydrogen | LEO 400 km'][0.75]:.0f}", "tsHbig": f"{CM['Hydrogen | LEO 400 km'][5.0]:.0f}",
    "tsOsmall": f"{CM['Oxygen | LEO 400 km'][0.75]:.1f}", "tsObig": f"{CM['Oxygen | LEO 400 km'][5.0]:.1f}",
    "tsMsmall": f"{CM['Methane | LEO 400 km'][0.75]:.1f}", "tsMbig": f"{CM['Methane | LEO 400 km'][5.0]:.1f}",
    "EnvSpread": f"{max(100 * (max(CM[f'{fl} | {e}'][r] for e in ENV) / min(CM[f'{fl} | {e}'][r] for e in ENV) - 1) for fl in FLUIDS for r in CM[f'{FLUIDS[0]} | LEO 400 km']):.0f}",
}
for f, tag in [("Hydrogen", "H"), ("Oxygen", "O")]:
    h, sh = HT[f], SH[f]
    macros.update({
        f"Qhold{tag}": f"{h['Q_leak']:.1f}", f"mliq{tag}": f"{h['m_liquid']:,.0f}".replace(",", "{,}"),
        f"dT{tag}": f"{h['dT']:.1f}", f"thom{tag}": f"{h['t_hold_homog_days']:.0f}",
        f"tcond{tag}": f"{h['t_hold_cond_days']:.0f}", f"fcond{tag}": f"{h['m_cond_frac']:.3f}",
        f"tcv{tag}": g3(sh['vented']), f"tcs{tag}": g3(sh['s3']),
        f"tcc{tag}": g3(sh['cond']), f"tch{tag}": g3(sh['hom']),
    })
tex = "% AUTO-GENERATED by scripts/make_figures.py - do not edit by hand\n"
tex += "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in macros.items())
(OUT.parent / "results_macros.tex").write_text(tex)
print(tex)
print(json.dumps(shift, indent=1))
