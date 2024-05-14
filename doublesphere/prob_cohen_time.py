"""
estimate the timescales copying a real badminton case
"""
import numpy as np
import matplotlib.pyplot as plt
#
# Time units
#
year = 3.1536e7    # Year [s]
hour = 3.6000e3    # Hour [s]
day = 8.6400e4    # Day  [s]
au = 1.496e13

# ==== parameters ====
u = 20e2 # [cm/s]
rho = 1.2 * 1e3 / 1e6 # [g/cm^3]
Mc = 3.0
Mb = 2.0
S = 28.0 # [cm^2]
Cd = 0.44
l_GC = 2.0 # [cm]

# aerodynamic length scale
L = 2 * (Mb + Mc) / rho / S / Cd
print('L = %f'%(L / 1e2))

# stabilizing time
#t_stable = 2 * Mb * (1 + Mb / Mc) / rho / S / Cd / u
t_stable = Mb / Mc * L / u
print('stable = %f'%(t_stable))

# oscillation time
t_osc = 2 * np.pi * np.sqrt(L * l_GC) / u
print('oscillation = %f'%(t_osc))

