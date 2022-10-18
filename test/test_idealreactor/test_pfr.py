import numpy as np

import reactord as rd


# Isothermic pfr
def test_fogler_p1_15a():
    """Fogler fourth ed. P1.15a"""

    fa_initial = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 3600  # m3/s

    k = 0.05 / 3600 / 0.001  # mol/s/m3

    v_pfr = 99 * 0.001  # m3

    def volume(temperature, pressure):
        return 1 / (fa_initial / f_volumetric)

    def kinetic(concentrations, temperature):
        return k

    substance_a = rd.Substance.from_thermo_database("methane")
    substance_a.volume_liquid = volume

    substance_b = substance_a

    mixture = rd.mix.IdealSolution([substance_a, substance_b])

    pfr = rd.idealreactor.pfr_classes.PfrHomogStatIsoth(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=np.array([-1, 1]),
        reactor_dims_minmax=[0, v_pfr],
        transversal_area=1,
        pressure=101325,
        reactor_isothermic_temperature=298.15,
        reactor_f_in=[fa_initial, 0],
        reactor_f_out=[np.nan, np.nan],
        kinetic_argument="concentration",
    )

    solution = pfr.simulate()

    fa_final = solution.y[0][-1]

    assert np.allclose(fa_final, 0.01 * fa_initial)

def test_fogler_p1_15b():
    """Fogler fourth ed. P1.15a"""

    fa_initial = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 3600  # m3/s

    k = 0.0001  # 1/s

    v_pfr = 127.9 * 0.001  # m3

    def volume(temperature, pressure):
        return 1 / (fa_initial / f_volumetric)

    def kinetic(concentrations, temperature):
        return k * concentrations[0]

    substance_a = rd.Substance.from_thermo_database("methane")
    substance_a.volume_liquid = volume

    substance_b = substance_a

    mixture = rd.mix.IdealSolution([substance_a, substance_b])

    pfr = rd.idealreactor.pfr_classes.PfrHomogStatIsoth(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=np.array([-1, 1]),
        reactor_dims_minmax=[0, v_pfr],
        transversal_area=1,
        pressure=101325,
        reactor_isothermic_temperature=298.15,
        reactor_f_in=[fa_initial, 0],
        reactor_f_out=[np.nan, np.nan],
        kinetic_argument="concentration",
    )

    solution = pfr.simulate()

    fa_final = solution.y[0][-1]

    assert np.allclose(fa_final, 0.01 * fa_initial)
