#!/usr/bin/env python3
'''Addition to test suite : get verbatim comments and graphical check ups associated with tests
'''
from os import mkdir
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np


from test_disk_physics import r, r_range, my_model, my_model_prime,\
    ref_length, sample0, sample1, conf, usr_list, DTDisk, evaluate_disk_mass_ratio

parser = ArgumentParser()
parser.add_argument('-d', dest='output_dir', type=str, help='select output directory')
args = parser.parse_args()


out = pathlib.Path(args.output_dir)
if out.exists():
    raise FileExistsError(out)
else:
    mkdir(out)
conf.write(out/'conf.nml')

center = my_model.get_value('r0')
sig = my_model.get_value('sig')

def eval_sample(f, sampler):
    r = sp.symbols('r')
    return np.array(list(map(lambda x: f.subs(r,x), sampler)))

def map_functions(functions:dict, sampler:list):
    ys = {}
    for k, func in functions.items():
        ys.update(
            { k : eval_sample(func, sampler) }
        )
    return ys

print('Quick physical diagnostics of the model')
print('=====================================================')

# From test_disk_phys.py
# ......................

disk_mass = evaluate_disk_mass_ratio(conf, my_model_prime)
Q = my_model.toomre_number.subs(r,center+sig).evalf(3)
print()
print(f'early disk/star mass ratio        {round(disk_mass*100,3)}%')
print(f"typical Toomre's Q                {Q}")
print(f'effective flaring index           {round(get_flaring(),3)}')



# From test_2D_statbility.py
# ..........................
"""make a plot of the LoveLace function"""
fig, axes = plt.subplots(nrows=3)

fs0 = {
    'rayleigh': my_model.rayleigh_function,
    'density' : my_model.surface_density,
    'omega'   : my_model.angular_velocity,
}
ys0 = map_functions(fs0, sample0)

fs1 = {
    'lovelace' : my_model.lovelace_function,
    'dlovelace': my_model.lovelace_function.diff(r),
}
ys1 = map_functions(fs1, sample1)

axes[0].plot(sample0, ys0['density'])
axes[0].set_ylabel(r'$\rho(r)$')

axes[1].plot(sample0, ys0['rayleigh'], label=r'$\partial_r (r^2\Omega)$')
xlims = axes[1].get_xlim()
axes[1].plot(np.linspace(*xlims, 2), np.zeros(2), ls=':', lw=0.5, c='k')
axes[1].set_xlim(xlims)

axes[2].plot(sample1, ys1['lovelace'], label=r'$L(r)$')
axes[2].plot(sample1, ys1['dlovelace'], ls='--', label=r'$\partial_r L(r)$')
for ax in axes:
    #mark the cavity location and add legends
    ylims = ax.get_ylim()
    ax.plot(usr_list['cavity_radius']*np.ones(2), np.linspace(*ylims, 2), ls=':', lw=0.5, c='k')
    ax.set_ylim(*ylims)
    ax.legend()
axes[-1].set_xlabel(r'$r$')
fig.savefig(str(out/'rayleigh_lovelace.eps'), dpi=900, bbox_tight=True)
fig.savefig(str(out/'rayleigh_lovelace.png'), bbox_tight=True)




# From test_dust.py
# .................
def local_upper_particle_size(model, R):
    '''Return maximum particle size for \rho_p = 1g/cm^3 at given radial position R.'''    
    rho = model.midplane_volumic_density
    cs  = model.sound_speed
    OmK = model.keplerian_angular_velocity
    ups = (cs * rho / OmK).subs(DTDisk.r,R) * sp.sqrt(2/sp.pi) / model.values[DTDisk.rhop]
    return ups.evalf()

def global_upper_particle_size(model, r_range):
    #In practice, often useless because for most models, the minimal
    #value of this function is given at either one of the boundaries.
    ups = 1.0
    for r in r_range:
        ups = min(ups, local_upper_particle_size(model, r))
    return ups

fig,ax = plt.subplots()
lups = np.zeros(len(r_range))
for i,rad in enumerate(r_range):
    lups[i] = local_upper_particle_size(my_model, rad)
ax.plot(r_range, lups)
ax.set_yscale('log')
ax.set_xlabel('$r$')
ax.set_ylabel('$s_{p,max}$ (code units)')

#add lines at 1mm, 1µm
xlims = ax.get_xlim()
for n in [1,4]:
    ax.plot(np.linspace(*xlims,2), np.ones(2)*10**(-n)/ref_length, ls='--', c='k', lw=0.5)
ax.set_xlim(*xlims)
ylims = ax.get_ylim()
ax.fill_between(r_range, ylims[1]*np.ones(len(lups)), lups, facecolor='red', alpha=0.2)
ax.set_ylim(*ylims)
ax.set_xlim(r_range.min(), r_range.max())

#twin axis with physical unit
axb = ax.twinx()
axb.set_yscale('log')
axb.set_ylim(*[y*ref_length for y in ax.get_ylim()])
axb.set_ylabel('$s_{p,max}$ (cm)')

fig.suptitle(r'''Maximum particle size allowed for $\rho_p=%.1f$g/cm$^3$''' % conf['usr_dust_list']['intrinsic_grain_density'])

fig.savefig(str(out/'max_grain_size.eps'), dpi=900, bbox_tight=True)
fig.savefig(str(out/'max_grain_size.png'), bbox_tight=True)
spmax = local_upper_particle_size(my_model, usr_list['cavity_radius'])
print(f'maximum grain size at cavity radius {spmax} (code units)')
plt.close(fig)


# Reynolds number for dust ------------------------
fig,ax = plt.subplots()
ax.plot(r_range,
        abs(eval_sample(my_model.dust_reynolds, r_range))
)

ax.set_yscale('log')
fig.savefig(str(out/'dust_reynolds.png'), bbox_tight=True)


# Stokes number -----------------------------------
fig,ax = plt.subplots()
fig.suptitle(r'''Initial Stokes number VS (Hersant 2009) criterion.
the $2\sigma$ wide region around the bump is displayed''')

for grain_size in conf['usr_dust_list']['grain_size']:
    sp = grain_size / ref_length
    St = np.zeros(len(r_range))
    my_model.values[DTDisk.sp] = sp
    for i,rv in enumerate(r_range):
        St[i] = my_model.stokes_number.subs(r, rv)
    ax.plot(r_range, St, label=f'$s_p={grain_size}$cm')
ax.set_yscale('log')
ax.plot(r_range, 0.5*np.ones(len(r_range)), color='k', lw=0.2)

ylims = ax.get_ylim()
ax.fill_between(
    r_range,
    ylims[1]*np.ones(len(r_range)),
    0.5*np.ones(len(r_range)),
    alpha=0.2,
    facecolor='red'
)


#annotate
asymp = dict(ls='--', c='k', lw=0.4)
ax.plot(r_range, r_range/r_range*0.5, **asymp)
sigmm = (center-sig/2)
sigpp = (center+3*sig/2)
ax.plot(sigmm*np.ones(2), np.array([1e-4,1e4]), **asymp)
ax.plot(sigpp*np.ones(2), np.array([1e-4,1e4]), **asymp)
ax.set_ylim(ylims)
ax.annotate(s=r'', xy=(sigmm,1), xytext=(sigpp,1), arrowprops=dict(arrowstyle='<->', shrinkA=0, shrinkB=0))
ax.annotate(s=r'$2\sigma$', xy=(center, 1.2))
ax.set_xlabel(r'$r$')
ax.set_ylabel(r'$\mathrm{St}$')
ax.set_xlim(r_range.min(), r_range.max())
ax.legend()
fig.savefig(str(out/'hersant.png'), bbox_tight=True)


# Toomre number -----------------------------------
fig,ax = plt.subplots()
fig.suptitle('Initial Toomre number $Q$ (SGI excpected where $Q < 1$)')
Q = np.zeros(len(r_range))
for i,rv in enumerate(r_range):
    Q[i] = my_model.toomre_number.subs(r, rv)
ax.plot(r_range, Q, label='$Q$ model')
ax.set_yscale('log')
ylims = ax.get_ylim()
ax.set_xlabel(r'$r$')
ax.set_ylabel(r'$Q$')

ax.fill_between(r_range, np.ones(len(r_range)), 1e-1*np.ones(len(r_range)), alpha=0.2, facecolor='red')
ax2 = ax.twinx()
omega = np.zeros(len(r_range))
kappa = np.zeros(len(r_range))
for i,rv in enumerate(r_range):
    omega[i] = my_model.angular_velocity.subs(r, rv)
    kappa[i] = my_model.epicyclic_frequency.subs(r, rv)
ax2.plot(r_range, omega, ls='--', lw=0)#skip a color
ax2.plot(r_range, omega, ls='--', label=r'$\Omega$')
ax2.plot(r_range, kappa, ls='--', label=r'$\kappa$')
ax2.set_ylabel('frequency')
ax2.legend()
fig.savefig(str(out/'toomre.png'), bbox_tight=True)


print()
print('=====================================================')
print(f'Graphic checks were added to {out}')
