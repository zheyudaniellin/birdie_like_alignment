"""
calculate representative timescales adoping a representative set of numbers
"""
import numpy as np
import matplotlib.pyplot as plt
import doublesphere
import pdb
import natconst

#
# Time units
#
year = 3.1536e7    # Year [s]
hour = 3.6000e3    # Hour [s]
day = 8.6400e4    # Day  [s]
au = 1.496e13

minute = hour / 60

# ==== parameters ====
# gas drift
v = 20e2 # [cm/s]

# gas density
rho_g = 1e-15 # [g / cm^3]

# thermal velocity
vth = 200e2 # [cm/s]

# grain size
l = 1e-1 # [cm]

# grain specific weight
rho_s = 3.0 # [g / cm^3]

# temperature
T = 20 # [K]

# ==== calculations ====
# stopping time
ts = rho_s * l / rho_g / vth 
print('stopping time = {:.1f} years'.format(ts / year))

# characteristic oscillation time
to = np.sqrt(rho_s / rho_g * l**2 / v / vth)
print('characteristic oscillation time = {:f} minute'.format(to/minute))

# characteristic damping ratio
q = np.sqrt(rho_g * vth / rho_s / v)
print('characteristic damping ratio = {:e}'.format(q))

epsilon = 1.01
kappa = 1.0
lam = 1.0

coef = doublesphere.coefficients(l, rho_s, epsilon, kappa, lam, rho_g, vth, v)
time = doublesphere.timescales(coef['I'], coef['D'], coef['P'])
print('hat_zeta = {:f}'.format(time['q'] / q))

# drift threshold
fac = lam**2 * epsilon**2 * (epsilon * kappa -1) / (1+epsilon**2)**2 / (1 + epsilon**2 * kappa)
vthres = np.sqrt(natconst.kk * T) / l**3 / rho_g * np.sqrt(2.3 * natconst.mp) / fac / (4*np.pi/3) / 2
print('vthres = {:f} [cm/s]'.format(vthres))

