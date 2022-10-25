import numpy as np

import pytest

import reactord as rd

# ======================================================================
# Fogler's Elements Of Chemical Reaction Engineering 4th ed. exercises
# p1.15 a, b and c
# ======================================================================

pytest.fixture(scope='module')
def reactive_mixture():
    fa_initial = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 60  # m3/s

    def volume(temperature, pressure):
        return 1 / (fa_initial / f_volumetric)

    substance_a = rd.Substance()
    substance_a.volume_liquid = volume

    substance_b = rd.Substance()
    substance_b.volume_liquid = volume

    mixture = rd.mix.IdealSolution([substance_a, substance_b])

    return mixture

pytest.fixture(scope='module')
def kinetic_a():
    k = 0.05 / 3600 / 0.001  # mol/s/m3
    def kinetic(concentrations, temperature):
        return k
    
    return kinetic

pytest.fixture(scope='module')
def kinetic_b():
    k = 0.0001  # 1/s
    def kinetic(concentrations, temperature):
        return k * concentrations[0]
    return kinetic

pytest.fixture(scope='module')
def kinetic_c():
    k = 3 / 3600 * 0.001  # 1/s
    def kinetic(concentrations, temperature):
        return k * concentrations[0] ** 2
    return kinetic

# Excercise a

def test_fogler_p1_15a_ivp(reactive_mixture, kinetic_a):
    """Fogler fourth ed. P1.15a as initial value problem"""

    fa_initial = 5 / 3600  # mol/s
    v_pfr = 99 * 0.001  # m3

    pfr = rd.idealreactor.pfr_classes.PfrHomogStatIsoth(
        mix=reactive_mixture,
        list_of_reactions=[kinetic_a],
        stoichiometry=np.array([-1, 1]),
        reactor_dims_minmax=[0, v_pfr],
        transversal_area=1,
        pressure=101325,
        reactor_isothermic_temperature=298.15,
        reactor_f_in=[fa_initial, 0],
        reactor_f_out=[np.nan, np.nan],
        kinetic_argument="concentration",
    )

    # Fogler's exact solution

    def fogler(vol):
        concentration = (fa_initial / f_volumetric) - vol * k / f_volumetric

        return concentration

    fogler_concentrations = np.array([])

    # reactord solution
    solution = pfr.simulate()
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[:, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)

def test_fogler_p1_15a_bvp():
    """Fogler fourth ed. P1.15a as border value problem"""

    fa_initial = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 60  # m3/s

    k = 0.05 / 3600 / 0.001  # mol/s/m3

    v_pfr = 99 * 0.001  # m3

    def volume(temperature, pressure):
        return 1 / (fa_initial / f_volumetric)

    def kinetic(concentrations, temperature):
        return k

    substance_a = rd.Substance()
    substance_a.volume_liquid = volume

    substance_b = rd.Substance()
    substance_b.volume_liquid = volume

    mixture = rd.mix.IdealSolution([substance_a, substance_b])

    # Fogler's exact solution

    def fogler(vol):
        concentration = (fa_initial / f_volumetric) - vol * k / f_volumetric
        return concentration

    fogler_concentrations = np.array([])


    fb_out = (fogler(0) - fogler(v_pfr)) * f_volumetric

    pfr = rd.idealreactor.pfr_classes.PfrHomogStatIsoth(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=np.array([-1, 1]),
        reactor_dims_minmax=[0, v_pfr],
        transversal_area=1,
        pressure=101325,
        reactor_isothermic_temperature=298.15,
        reactor_f_in=[fa_initial, np.nan],
        reactor_f_out=[np.nan, fb_out],
        kinetic_argument="concentration",
    )

    # ==================================================================
    # reactord solution border value problem
    # ==================================================================

    solution = pfr.simulate()
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[:, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)

def test_fogler_p1_15b_ivp():
    """Fogler fourth ed. P1.15b"""

    fa_initial = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 60  # m3/s

    k = 0.0001  # 1/s

    v_pfr = 127.9 * 0.001  # m3

    def volume(temperature, pressure):
        return 1 / (fa_initial / f_volumetric)

    def kinetic(concentrations, temperature):
        return k * concentrations[0]

    substance_a = rd.Substance()
    substance_a.volume_liquid = volume

    substance_b = rd.Substance()
    substance_b.volume_liquid = volume

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

    # Fogler's exact solution

    def fogler(vol):
        concentration = (fa_initial / f_volumetric) / np.exp(
            vol * k / f_volumetric
        )
        return concentration

    fogler_concentrations = np.array([])

    # reactord solution
    solution = pfr.simulate()
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[:, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)

def test_fogler_p1_15c_ivp():
    """Fogler fourth ed. P1.15b"""

    fa_initial = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 60  # m3/s

    k = 3 / 3600 * 0.001  # 1/s

    v_pfr = 660 * 0.001  # m3

    def volume(temperature, pressure):
        return 1 / (fa_initial / f_volumetric)

    def kinetic(concentrations, temperature):
        return k * concentrations[0] ** 2

    substance_a = rd.Substance()
    substance_a.volume_liquid = volume

    substance_b = rd.Substance()
    substance_b.volume_liquid = volume

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

    # Fogler's exact solution
    def fogler(vol):
        concentration = 1 / (
            vol * k / f_volumetric + 1 / (fa_initial / f_volumetric)
        )
        return concentration

    fogler_concentrations = np.array([])

    # reactord solution
    solution = pfr.simulate()
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[:, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)
