from test_disk_physics import my_model, conf, r

def test_rho_min():
    rmin = conf['meshlist']['xprobmin1']
    rho_min = my_model.surface_density.subs(r, rmin)
    assert rho_min > 5*conf['usr_list']['rhomin']
