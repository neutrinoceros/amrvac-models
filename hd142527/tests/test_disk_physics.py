'''Analytical tests on basic assumptions such as the disk to star mass ratio.'''

import pathlib

import numpy as np
import sympy as sp

from amrvac_pywrap import merge_configs
from symdisk import Disk, TransitionDisk as TDisk, DustyDisk as DDisk

here = pathlib.Path(__file__).absolute().parent.parent

# units and convertion factors definitions
# ----------------------------------------
au2m = 1.49597870691e11
m2au = au2m**-1
m2cm = 1e2
cm2au = m2cm**-1 * m2au

kg2g = 1e3

s2yr = 3.1536e7

Msun       = 1.988e33                 #in g
ref_length = au2m * m2cm              #1au in cm
ref_time   = s2yr / (2*np.pi)         #1yr/2pi in s
ref_sigma  = Msun / ref_length**2
ref_rho    = ref_sigma / ref_length

h2_viscosity_si   = 2e-6 #in N.s/m^2. ~ tabulated value for molecular hydrogen @ T=50K
h2_viscosity_cgs  = h2_viscosity_si * kg2g / m2au / s2yr
h2_viscosity_code = h2_viscosity_cgs / Msun * ref_time * ref_length

# build the disk model object from configuration files
# ----------------------------------------------------
class DTDisk(DDisk, TDisk):
    '''Simply merge properties from DustyDisk and TransitionDisk'''
    pass


def build_model(conf, model_class):
    usr_list = conf['usr_list']

    if isinstance(conf['usr_dust_list']['grain_size'], list):
        maxsize = max(conf['usr_dust_list']['grain_size'])
    else:
        maxsize = conf['usr_dust_list']['grain_size']
    my_subs = {
        DTDisk.star_mass     : conf['disk_list']['central_mass'],
        DTDisk.sig           : usr_list['cavity_width'],
        DTDisk.r0            : usr_list['cavity_radius'],
        DTDisk.rho0          : usr_list['rhozero'],
        DTDisk.slope         : usr_list['density_slope'],
        DTDisk.aspect_ratio0 : usr_list['aspect_ratio'],
        DTDisk.gamma         : conf['hd_list']['hd_gamma'],
        DTDisk.aspect_ratio0 : usr_list['aspect_ratio'],
        DTDisk.eta           : h2_viscosity_code,
        DTDisk.sp            : maxsize / ref_length,
        DTDisk.rhop          : conf['usr_dust_list']['intrinsic_grain_density'] / Msun * ref_length**3
    }
    return model_class(my_subs)


conf = merge_configs([here/'hd142527.nml', here/'add_dust.par'])

r_range = np.linspace(conf['meshlist']['xprobmin1'], conf['meshlist']['xprobmax1'], num=100)

my_model  = build_model(conf, DTDisk)

#an 'earlier' version of the same disk, without a cavity
my_model_prime = build_model(conf, Disk)

r = sp.symbols('r')
x = sp.symbols('x')

usr = conf['usr_list']
sample0 = np.linspace(conf['meshlist']['xprobmin1'], conf['meshlist']['xprobmax1'], 100)
sample1 = np.linspace(usr['cavity_radius'] - 2*usr['cavity_width'],
                      usr['cavity_radius'] + 2*usr['cavity_width'],
                      100)


# definitions
# -----------
def evaluate_fiducial_mass(conf, model):
    rmin= conf['usr_list']['cavity_radius']
    rmax = conf['meshlist']['xprobmax1']
    integrand = 2*sp.pi*model.surface_density*r
    fid_mass = sp.integrate(integrand, (r, rmin, rmax))
    return fid_mass.evalf()

def evaluate_disk_mass_ratio(conf, model):
    rmin, rmax = conf['meshlist']['xprobmin1'], conf['meshlist']['xprobmax1']
    integrand = 2*sp.pi*model.surface_density*r
    disk_mass = sp.integrate(integrand, (r, rmin, rmax))
    return disk_mass.evalf() / conf['disk_list']['central_mass']

def get_flaring(alpha, gamma):
    return 0.5 * (1 + alpha*(1-gamma))



# basic tests
# -----------
def test_neglible_self_gravity():
    '''Evaluate the Toomre criterion (Toomre 1964) for self-gravity.

    SG is considered neglible when Q = c_s\Omega/(\pi G \Sigma) < 1
    '''
    models = my_model_prime, my_model

    bounds = (conf['meshlist']['xprobmin1'],conf['meshlist']['xprobmax1'])
    sampler = np.linspace(*bounds, 100)
    for model in models:
        for s in sampler:
            assert model.toomre_number.subs(r,s) > 1

def test_disk_mass_ratio():
    '''Check that total disk mass is less than 30% of the mass star.

    Star mass is assumed unity.
    We test the powerlaw distribution itself, without taking the cavity into account.
    This test is somewhat arbitrary and should not be considered too important as long as
    self gravity is negligle (see previous test).
    '''
    assert evaluate_disk_mass_ratio(conf, my_model_prime) < 0.30

def test_fiducial_mass():
    '''Compare an estimation of the disk mass with observations (20% margin)'''
    obs = 1e-1 #disk mass from obervations
    eff = evaluate_fiducial_mass(conf, my_model_prime)
    print('\n\n\n', eff, obs, '\n\n')
    assert abs(eff-obs)/obs < 2e-1

def test_effective_flaring_index():
    params = dict(
        alpha = usr['density_slope'],
        gamma = conf['hd_list']['hd_gamma']
    )
    assert not get_flaring(**params) < 0


# stability tests
# ---------------
def test_rayleigh_stability():
    for x in sample0:
        assert my_model.rayleigh_function.subs(r,x) > 0

def test_rwi_instability():
    '''Check that initial setup is RWI UNstable.'''
    LoveP  = my_model.lovelace_function.diff(r)
    bounds = (
        my_model.values[TDisk.r0] - 2*my_model.values[TDisk.sig],
        my_model.values[TDisk.r0] + 2*my_model.values[TDisk.sig]
    )
    try:
        #search for a maximum of LovelaceFunction near r0
        root = sp.nsolve(LoveP, bounds, solver='bisect')
        delta = (root-my_model.values[TDisk.r0])
        print(f"\n found root at {root}, delta = {delta}", end =', ')
        print(f"ratio = {abs(delta/my_model.values[TDisk.r0])}\n")
    except ValueError:
        #raised when a root can't be found in the interval
        delta = sp.nan
    assert(delta is not sp.nan)
    assert(abs(delta/my_model.values[TDisk.r0]) <= 0.2)#assert the extremum is close to r0 with a 20% margin


# dust related test
# -----------------
def test_no_turbulence_around_grains():
    for rv in r_range:
        assert my_model.dust_reynolds.subs(r, rv) < 1e3

def test_hersant():
    '''Check for (Hersant 2009)'s criterion for pressure_less approximation validity'''
    margin = 0.1
    center = my_model.get_value('r0')
    sig = my_model.get_value('sig')

    rmin = center-sig/2
    rmax = center+sig*3/2
    rs = np.linspace(rmin, rmax, 100)
    St = np.zeros(len(rs))
    for i,rv in enumerate(rs):
        St[i] = my_model.stokes_number.subs(r, rv)
        assert St[i] < 0.5 + margin
