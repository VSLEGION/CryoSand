# CryoSand

Open trade-space sandbox for in-space cryogenic propellant storage. See PROJECT.md for the
physics specification and conventions.

## Setup (laptop)
    pip install -r requirements.txt
    python -m pytest                 # all tests; CoolProp is required, CI enforces it
    python scripts/make_figures.py   # regenerates report/figures, results.json, results_macros.tex
    cd report && latexmk -pdf CryoSand_FirstReport.tex

Every number in the report comes from `make_figures.py` through `results_macros.tex`, and every
fluid property comes from CoolProp. The script refuses to run without CoolProp.

## Layout
    cryosand/core      constants, CoolProp wrapper, resistance network, dataclasses
    cryosand/layers    L0-L5 physics (pure functions)
    cryosand/model.py  steady assembly + closed-tank bounds (CoolProp)
    cryosand/trade.py  passive vs ZBO optimisation, crossover
    tests/             test_limits.py (property-free), test_coolprop.py
    data/parameters.yaml   parameter register (sourced / provisional)
    report/            LaTeX first written report
