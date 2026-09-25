"""Physical constants. SI base units only.

Every value carries its source. Fluid properties do NOT live here; they come
from CoolProp via cryosand.core.properties.
"""

# Stefan-Boltzmann constant [W/(m^2 K^4)]. CODATA 2018 exact value.
SIGMA = 5.670374419e-8

# Total solar irradiance at 1 AU [W/m^2]. Kopp & Lean (2011), GRL 38, L01706.
G_SUN_1AU = 1361.0

# Earth mean outgoing long-wave radiation [W/m^2]. Gilmore (ed.), Spacecraft
# Thermal Control Handbook, Vol. I, 2nd ed. (2002), ch. 2 (orbit-average ~237).
G_IR_EARTH = 237.0

# Earth mean Bond albedo [-]. Gilmore (2002), ch. 2 (~0.30).
ALBEDO_EARTH = 0.30

# Mean Earth radius [m]. IUGG mean radius.
R_EARTH = 6.371e6

# Effective deep-space sink temperature [K] (cosmic background ~2.7 K; ~4 K
# is the conventional engineering value). Gilmore (2002).
T_DEEP_SPACE = 4.0

# Seconds per day [s/day]. Display-layer helper; never used inside physics.
SECONDS_PER_DAY = 86400.0
