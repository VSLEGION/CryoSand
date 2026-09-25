# CryoSand

Open trade-space sandbox for in-space cryogenic propellant storage. See PROJECT.md for the
physics specification and conventions.

## Setup (laptop)
    pip install numpy scipy matplotlib CoolProp pytest
    pytest -q                        # all tests, including CoolProp-dependent bounds tests
    python scripts/make_figures.py   # regenerates report/figures, results.json, results_macros.tex
    cd report && latexmk -pdf CryoSand_FirstReport.tex

`make_figures.py` falls back to PROJECT.md reference values only when CoolProp is missing,
and says so in the report (\PropSource macro). With CoolProp installed, every number in the
report is regenerated from CoolProp automatically.

## Layout
    cryosand/core      constants, CoolProp wrapper, resistance network, dataclasses
    cryosand/layers    L0-L5 physics (pure functions)
    cryosand/model.py  steady assembly + closed-tank bounds (CoolProp)
    cryosand/trade.py  passive vs ZBO optimisation, crossover
    tests/             test_limits.py (property-free), test_coolprop.py
    data/parameters.yaml   parameter register (sourced / provisional)
    report/            LaTeX first written report
