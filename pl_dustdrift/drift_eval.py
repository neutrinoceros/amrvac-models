from vtk_vacreader import VacDataSorter

import numpy as np
import matplotlib.pyplot as plt
import f90nml

au2cm = 8*60*3e10
msun2g = 2e33

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
        #self.S = conf['usr_list']['aspect_ratio']**2 * conf['disk_list']['central_mass'] * conf['usr_list']['rhozero']**(1-self.gamma) / self.gamma
        self.S = conf['hd_list']['hd_adiab']
        self.A = 4/3 *np.sqrt(self.gamma/3)/(3/4)

        self.conf = conf
        self.data = dataholder

    def get_sound_speed(self):
        return np.sqrt(self.gamma * self.S * self.data['rho']**(self.gamma-1))

    def get_stopping_time(self, dust_index:int, Akey:str='amrvac'):
        rhog = self.data['rho']
        sp = self.grain_sizes[dust_index]
        cs = self.get_sound_speed()
        ts = self.A * self.rho_p * sp / (rhog * cs)
        return ts

    def get_keplerian_pulsation(self):
        r = self.data.get_ticks(0)
        OmegaK = np.sqrt(self.conf['disk_list']['central_mass'] / r**3)
        return OmegaK

    def get_stokes(self, dust_index:int, stokes_method='base'):
        if stokes_method == 'base':
            ts = self.get_stopping_time(dust_index)
            OmegaK = self.get_keplerian_pulsation()
            St = OmegaK * ts

        elif stokes_method == 'mod':
            sigmag = self.data['rho']
            sp = self.grain_sizes[dust_index]
            St = np.sqrt(2) * self.A * self.rho_p*sp / sigmag
        return St

    def get_theo_r_drift(self, dust_index:int, stokes_method='base'):
        r = self.data.get_ticks(0)
        Omega_K = self.get_keplerian_pulsation()
        vK = Omega_K * r
        vtg = self.data['m2'] / self.data['rho']
        eta = (vK - vtg)/ vK
        St = self.get_stokes(dust_index, stokes_method)
        return -2 * eta * Omega_K * r * St / (1 + St**2)


# ----------------------------------------------------------------------
fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, figsize=(10,10))

titles = [
    r'$v_r$',
    r'$\delta v_r$ (d-g)',
    r'$\dot{r}$ th',
    r'$\dot{r}$ th (mod)',
    r'$St$ th',
    r'$St$ th (mod)',
]
for ax, tit in zip(axes.flatten(), titles):
    ax.set_title(tit)

dh = VacDataSorter(f'out_1D/pl_drift0010.vtu')
tc = TheoCrusher('conf1D.nml', dh)

vrg = dh['m1']/dh['rho']
rvect = dh.get_ticks()
lss = ['--', ':', '-.', '-', '--']
for i, ls in enumerate(lss):
    vrd = dh[f'm1d{i+1}']/dh[f'rhod{i+1}']
    delta = vrd-vrg
    qties = [
        delta,
        tc.get_theo_r_drift(i),
        tc.get_theo_r_drift(i, 'mod'),
        tc.get_stokes(i),
        tc.get_stokes(i, 'mod'),
    ]
    
    axes[0,0].plot(rvect, vrd, ls=ls, label=str(i+1))
    for ax, qty in zip(axes.flatten()[1:], qties):
        ax.plot(rvect, qty, ls=ls)

axes[0,0].plot(rvect, vrg, c='k')
axes[0,0].legend()
for ax in axes[-1]:
    ax.set_yscale('log')

fig.savefig(f'{__file__}.png')

#dh2D = VacDataSorter(f'out_kwok/pl_drift0070.vtu', data_shape=(512,128))
