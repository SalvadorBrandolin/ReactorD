import numpy as np
import pytest

import reactord as rd

# ======================================================================
# Fogler's Elements Of Chemical Reaction Engineering 4th ed.
# p1.15
# ======================================================================


def test_fogler_p1_15a_ivp():
    """Fogler fourth ed. P1.15a as an initial value problem"""

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
    """Fogler fourth ed. P1.15a as a border value problem"""

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

    # flux of B on the reactor's outlet

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


def test_fogler_p1_15a_ivp_outlet():
    """Fogler fourth ed. P1.15a as a initial value problem, but initial
    conditions are in the reactor's outlet"""

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

    # flux of A and B on the reactor's outlet

    fa_out = fogler(v_pfr) * f_volumetric
    fb_out = (fogler(0) - fogler(v_pfr)) * f_volumetric

    pfr = rd.idealreactor.pfr_classes.PfrHomogStatIsoth(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=np.array([-1, 1]),
        reactor_dims_minmax=[0, v_pfr],
        transversal_area=1,
        pressure=101325,
        reactor_isothermic_temperature=298.15,
        reactor_f_in=[np.nan, np.nan],
        reactor_f_out=[fa_out, fb_out],
        kinetic_argument="concentration",
    )

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


def test_fogler_p1_15b():
    """Fogler fourth ed. P1.15b as initial value problem"""

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


def test_fogler_p1_15b_bvp():
    """Fogler fourth ed. P1.15b as initial value problem"""

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

    # Fogler's exact solution

    def fogler(vol):
        concentration = (fa_initial / f_volumetric) / np.exp(
            vol * k / f_volumetric
        )
        return concentration

    fogler_concentrations = np.array([])

    # Out molar flux of B to border condition

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

    # import ipdb
    # ipdb.set_trace()
    # reactord solution
    solution = pfr.simulate(tol=0.00001)
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[:, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)


def test_fogler_p1_15b_ivp_outlet():
    """Fogler fourth ed. P1.15b as initial value problem but initial
    conditions are on the reactor's outlet"""

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

    # Fogler's exact solution

    def fogler(vol):
        concentration = (fa_initial / f_volumetric) / np.exp(
            vol * k / f_volumetric
        )
        return concentration

    fogler_concentrations = np.array([])

    # Out molar flux of A and B to border condition

    fa_out = fogler(v_pfr) * f_volumetric
    fb_out = (fogler(0) - fogler(v_pfr)) * f_volumetric

    pfr = rd.idealreactor.pfr_classes.PfrHomogStatIsoth(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=np.array([-1, 1]),
        reactor_dims_minmax=[0, v_pfr],
        transversal_area=1,
        pressure=101325,
        reactor_isothermic_temperature=298.15,
        reactor_f_in=[np.nan, np.nan],
        reactor_f_out=[fa_out, fb_out],
        kinetic_argument="concentration",
    )

    # import ipdb
    # ipdb.set_trace()
    # reactord solution
    solution = pfr.simulate(tol=0.00001)
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
    """Fogler fourth ed. P1.15c as an initial value problem"""

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


def test_fogler_p1_15c_bvp():
    """Fogler fourth ed. P1.15c as a border value problem"""

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

    # Fogler's exact solution
    def fogler(vol):
        concentration = 1 / (
            vol * k / f_volumetric + 1 / (fa_initial / f_volumetric)
        )
        return concentration

    fogler_concentrations = np.array([])

    # B molar flux in the reactor's outlet for the border conditions

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

    # reactord solution
    solution = pfr.simulate(tol=1e-6)
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[:, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)


def test_fogler_p1_15c_ivp_outlet():
    """Fogler fourth ed. P1.15c as an initial value problem but the
    initial conditions are on the reactor's outlet"""

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

    # Fogler's exact solution
    def fogler(vol):
        concentration = 1 / (
            vol * k / f_volumetric + 1 / (fa_initial / f_volumetric)
        )
        return concentration

    fogler_concentrations = np.array([])

    # A and B molar flux in the reactor's outlet for the border conditions
    fa_out = fogler(v_pfr) * f_volumetric
    fb_out = (fogler(0) - fogler(v_pfr)) * f_volumetric

    pfr = rd.idealreactor.pfr_classes.PfrHomogStatIsoth(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=np.array([-1, 1]),
        reactor_dims_minmax=[0, v_pfr],
        transversal_area=1,
        pressure=101325,
        reactor_isothermic_temperature=298.15,
        reactor_f_in=[np.nan, np.nan],
        reactor_f_out=[fa_out, fb_out],
        kinetic_argument="concentration",
    )

    # reactord solution
    solution = pfr.simulate(tol=1e-6)
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[:, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)


# ======================================================================
# Fundamentals of Chemical Reaction Engineering. Mark E. Davis;
# Robert J. Davis.
# Example 3.4.1 - page 78
# ======================================================================


@pytest.mark.skip
def test_davis_3_4_1():
    return None
