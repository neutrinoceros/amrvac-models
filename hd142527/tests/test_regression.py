"""Some integration tests based on config1D.par.

Run similar but modified versions of the basic simulation.

Warning, requirements
=====================
 * Python3.6+ with pytest.
 * $AMRVAC_DIR/tools/py3tools/ must be part of your $PYTHONPATH.
 * $AMRVAC_TESTING_ARCH must be defined.
"""

import time
import pathlib

import numpy as np
import matplotlib.pyplot as plt

import f90nml
from amrvac_pywrap import ForkSimulation, merge_configs
from vtk_vacreader import VacDataSorter
import pytest

# Globals
# -------
here = pathlib.Path(__file__).absolute().parent
defs_rep = here/'data/defs'
gold1D_pref = str(here/'data/gold/ref1D/cavity')
gold2D_pref = str(here/'data/gold/ref2D/cavity2D')

tests_dir = here/f"LOG_TESTS_{time.strftime('%H%M%S')}"
tests_dir.mkdir()

dirdefs = dict(origin=here.parent, outdir=tests_dir)

# Loading reference snapshots
# ---------------------------
s0 = VacDataSorter(file_name=gold1D_pref+'0000.vtu')
s1 = VacDataSorter(file_name=gold1D_pref+'0010.vtu')
s2D0 = VacDataSorter(file_name=gold2D_pref+'0000.vtu', shape=None) #can't remember the shape :/
s2D1 = VacDataSorter(file_name=gold2D_pref+'0010.vtu', shape=None)

def m2v(d):
    d.fields['v1'] = d['m1']/d['rho']
    d.fields['v2'] = d['m2']/d['rho']

for snap in [s0,s1,s2D0,s2D1]: m2v(snap)
def err(out0, out1, k):
    return np.mean(np.abs((out1.fields[k]-out0.fields[k])/out0.fields[k]))


# 1D tests
# ========
#@pytest.mark.skip(reason='tmp')
def test_rerun():
    conf = f90nml.read(defs_rep/'config1D.par') #Pure integration test : no modifications.
    sim = ForkSimulation(conf=conf, tag='rerun', **dirdefs)
    assert(sim.run() == 0)
    dat = sim.load_data(offset=10)
    m2v(dat)
    for key in ['rho','v2']:
        assert(err(s0,dat,key)<=err(s0,s1,key))
        assert(err(s0,dat,key)<=1e-6)

def test_naive_boundary():
    import os
    conf = [f90nml.read(defs_rep/'config1D.par'), f90nml.read(defs_rep/'naivebc.par')]
    #devnote : I'm using a ready-made file here instead of a
    # dictionnary because I have issues with correctly writting out
    # str values to f90nml.Namelist
    sim = ForkSimulation(conf=conf, tag='naive_bc', **dirdefs)
    assert(sim.run() == 0)
    dat = sim.load_data(offset=10)
    err_ref = err(s0,s1,'rho')
    err_new = err(s0,dat,'rho')
    assert(err_new > err_ref)
    return dat

# 2D tests
# ========
def test_rerun_2D():
    conf = [f90nml.read(defs_rep/'config1D.par'), f90nml.read(defs_rep/'to2D.par')]
    sim = ForkSimulation(conf=conf, tag='rerun2D', dim=2, **dirdefs)
    exitcode = sim.run()
    dat = sim.load_data(offset=10, reshape=False)
    m2v(dat)
    for key in ['rho','v2']:
        assert(err(s2D0,dat,key)<=err(s2D0,s2D1,key))

def test_nosplit_2D():
    baseconf = f90nml.read(defs_rep/'config1D.par')
    custom = {'meshlist':
              {
                  'block_nx1':baseconf['meshlist']['domain_nx1']
              }
    }
    conf = [f90nml.read(defs_rep/'config1D.par'), f90nml.read(defs_rep/'to2D.par'), custom]
    sim = ForkSimulation(conf=conf, tag='nosplit2D', dim=2, **dirdefs)
    assert(sim.run() == 0)
    dat = sim.load_data(offset=10)
    return dat

if __name__=='__main__':
    print("run me with pytest")
    print("==================\n")
