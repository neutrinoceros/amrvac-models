import f90nml
import numpy as np
import matplotlib.pyplot as plt

from vtk_vacreader import VacDataSorter


conf = f90nml.read('conf1D.nml')
ds = VacDataSorter('out/pl_drift0000.vtu')

try:
    ri = conf['disk_list']['ref_radius']
except KeyError:
    ri = 1
mstar = conf['disk_list']['central_mass']

G = 4*np.pi**2 * ri**3/mstar

rv = ds.get_ticks()
vK = np.sqrt(G*mstar/rv)
vtg = ds['m2']/ds['rho']
vtd  = ds['m2d1']/ds['rhod1']

fig, ax = plt.subplots()
ax.plot(rv, (vtd-vK)/vK, label='dust')
#ax.plot(rv, (vtg-vK)/(2*np.pi), label='gas')


ax.set_xlabel(r'$r$')
ax.set_ylabel(r'$(v_\varphi - v_K)/2\pi$')
ax.set_xscale('log')
ax.grid()
ax.legend()
fig.savefig(f'{__file__[:-3]}', tight_bbox=True)
