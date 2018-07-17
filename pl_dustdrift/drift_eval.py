#!/usr/bin/env python3

'''Plot the expected and effective drift rates (g/d velocity
discrepency) and different evaluations of the Stokes numbers (with
respect to r)
'''
from argparse import ArgumentParser
import itertools as itt
from vtk_vacreader import VacDataSorter as VDS

import numpy as np
import matplotlib.pyplot as plt
import f90nml

# raw conversion factors
au2cm = 8*60*3e10 #1au (cm)
msun2g = 2e33


class VelVDS(VDS):
    '''Replace moments with velocities'''
    def __init__(self, file_name:str, data_shape:tuple=None):
        super().__init__(file_name, data_shape)
        keys = [k for k in self.fields.keys()]
        for k in keys:
            if 'm' in k:
                tag = k[1:]
                vkey = 'v' + tag
                self.fields.update({vkey: self[k]/self['rho'+tag[1:]]})
                self.fields.pop(k)

class TheoCrusher:
    '''Use numerical data on gas velocity to predict the dust drift
    based on equations 13, 14 in (Chiang & Youdin 2010)'''

    def __init__(self, namelist_file, dataholder):
        conf = f90nml.read(namelist_file)
        grain_sizes_cm = conf['usr_dust_list']['grain_size_cm']
        if isinstance(grain_sizes_cm, list):
            self.grain_sizes = np.array(grain_sizes_cm) / au2cm
        else:
            self.grain_sizes = np.array([grain_sizes_cm]) / au2cm
        self.rho_p = conf['usr_dust_list']['grain_density_gcm3'] / msun2g * au2cm**3
        self.gamma = conf['hd_list']['hd_gamma']
        self.S = conf['hd_list']['hd_adiab']
        self.A = 4/3 *np.sqrt(self.gamma/3)/(3/4)

        try:
            self.ri = conf['disk_list']['ref_radius']
        except KeyError:
            self.ri = 1
        self.mstar = conf['disk_list']['central_mass']
        self.G = 4*np.pi**2 * self.ri**3 / self.mstar

        self.conf = conf
        self.data = dataholder

    def get_sound_speed(self):
        '''Get dimensionless soundspeed from dimensionless SURFACE density'''
        return np.sqrt(self.gamma * self.S * self.data['rho']**(self.gamma-1))

    def get_stopping_time(self, dust_index:int):
        '''dimensionless stopping time'''
        rhog = self.data['rho']
        sp = self.grain_sizes[dust_index]
        cs = self.get_sound_speed()
        ts = self.A * self.rho_p * sp / (rhog * cs)
        return ts

    def get_keplerian_pulsation(self):
        r = self.data.get_ticks(0)
        OmegaK = np.sqrt(self.G * self.mstar / r**3)
        return OmegaK

    def get_stokes(self, dust_index:int, stokes_method='base'):
        if stokes_method == 'base':
            ts = self.get_stopping_time(dust_index)
            OmegaK = self.get_keplerian_pulsation()
            St = OmegaK * ts

        elif stokes_method == 'mod':
            sigmag = self.data['rho']
            sp = self.grain_sizes[dust_index]
            St = np.sqrt(2*np.pi) * self.A * self.rho_p*sp / sigmag
        return St

    def get_theo_r_drift(self, dust_index:int, stokes_method='base'):
        '''from (Chiang & Youdin 2010) eq 13'''
        r = self.data.get_ticks(0)
        Omega_K = self.get_keplerian_pulsation()
        vtg = self.data['v2']
        eta = (vK - vtg) / vK
        St = self.get_stokes(dust_index, stokes_method)
        return -2 * eta * Omega_K * r * St / (1 + St**2)

    def get_theo_phi_drift(self, dust_index:int, stokes_method='base'):
        '''from (Chiang & Youdin 2010) eq 14'''
        St = self.get_stokes(dust_index, stokes_method)
        return self.get_theo_r_drift(dust_index, stokes_method) / (2*St)

# ----------------------------------------------------------------------

ap = ArgumentParser()
ap.add_argument('-o', dest='out', type=str, default='.')
args = ap.parse_args()
out = args.out

offs = (0, 1, 50, 501)
for n in offs:
    dh = VelVDS(file_name=f'{out}/pl_drift{str(n).zfill(4)}.vtu')
    fig, axes = plt.subplots(nrows=3, ncols=4, sharex=True, figsize=(20,10))
    tc = TheoCrusher('conf1D.nml', dh)

    rvect = dh.get_ticks()
    vK = tc.get_keplerian_pulsation()*rvect
    vrg = dh['v1']
    vphig = dh['v2']
    dust_n_species = len([k for k in dh.fields.keys() if 'd' in k and 'rho' in k])
    lss = ['--', ':', '-.', '-', '--']

    for i, ls in zip(range(dust_n_species), lss):
        vrd   = dh[f'v1d{i+1}']
        vphid = dh[f'v2d{i+1}']
        delta_vr = vrd-vrg
        delta_vphi = (vphid - vphig) / rvect
        rhod =  dh[f'rhod{i+1}']
        qties = [
            (r'$v_r$', vrd),
            (r'$\delta v_r$ (d-g)', delta_vr),
            (r'$(v_\varphi - v_K)$', (vphid - vK)),
            (r'$(\delta v_\varphi)/r$ (d-g)', delta_vphi),
            (r'$\dot{r}$ th', tc.get_theo_r_drift(i)),
            (r'$\dot{r}$ th (mod)', tc.get_theo_r_drift(i, 'mod')),
            (r'$\dot{\varphi}$ th', tc.get_theo_phi_drift(i)),
            (r'$\dot{\varphi}$ th (mod)', tc.get_theo_phi_drift(i, 'mod')),
            (r'$\mathrm{St}$ th', tc.get_stokes(i)),
            (r'$\mathrm{St}$ th (mod)', tc.get_stokes(i, 'mod')),
            (r'rho dust', rhod),
        ]

        for ax, (tit, field) in zip(axes.flatten(), qties):
            ax.plot(rvect, field, ls=ls, label=str(i+1))
            ax.set_title(tit)

    axes[0,0].plot(rvect, vrg, c='k')
    axes[0,0].legend()

    axes[0,2].plot(rvect, (vphig-vK), c='k')
    for ax in itt.chain(axes[0:2,2], axes[0:2,3]): #vphi
        ax.set_ylim(-3e-4, 1e-4)

    #vr
    #axes[0,0].set_ylim(-1.2e-7, 1e-8)
    axes[0,1].set_ylim(-6e-8, 1e-8)

    #axes[2,2].set_yscale('linear')
    #axes[2,2].plot(rvect, tc.get_sound_speed()/(rvect*tc.get_keplerian_pulsation()))
    axes[2,2].plot(rvect, dh['rho'], c='k')
    for ax in axes[-1]:
        ax.set_yscale('log')

    fig.savefig(f'diag{n}.png')
