import numpy as np

import reactord as rd

# ======================================================================
# Fogler's Elements Of Chemical Reaction Engineering 4th ed.
# p1.15
# ======================================================================


def test_fogler_p1_15a_ivp():
    """Fogler fourth ed. P1.15a as an initial value problem"""

    def volume(temperature, pressure):
        return 1 / (fa_in / f_volumetric)

    def kinetic(concentrations, temperature):
        return k

    # Fogler's exact solution

    def fogler(vol):
        concentration = (fa_in / f_volumetric) - vol * k / f_volumetric

        return concentration

    fogler_concentrations = np.array([])

    fa_in = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 60  # m3/s

    k = 0.05 / 3600 / 0.001  # mol/s/m3

    v_pfr = 99 * 0.001  # m3

    substance_a = rd.Substance(name="A")
    substance_a.volume_liquid = volume

    substance_b = rd.Substance(name="B")
    substance_b.volume_liquid = volume

    mixture = rd.mix.IdealSolution(A=substance_a, B=substance_b)

    pfr = rd.idealreactor.StationaryPFR.set_isothermic_isobaric(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=[-1, 1],
        kinetic_argument="concentration",
        reactor_dim_minmax=[0, v_pfr],
        transversal_area=1,
        isothermic_temperature=298.15,
        isobaric_pressure=101325,
        molar_flow_in={"A": fa_in, "B": 0},
    )

    # reactord solution
    solution = pfr.simulate(grid_size=100)
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[0:2, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)


def test_fogler_p1_15a_bvp():
    """Fogler fourth ed. P1.15a as a border value problem"""

    def volume(temperature, pressure):
        return 1 / (fa_in / f_volumetric)

    def kinetic(concentrations, temperature):
        return k

    # Fogler's exact solution

    def fogler(vol):
        concentration = (fa_in / f_volumetric) - vol * k / f_volumetric

        return concentration

    fa_in = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 60  # m3/s

    k = 0.05 / 3600 / 0.001  # mol/s/m3

    v_pfr = 99 * 0.001  # m3

    substance_a = rd.Substance(name="A")
    substance_a.volume_liquid = volume

    substance_b = rd.Substance(name="B")
    substance_b.volume_liquid = volume

    mixture = rd.mix.IdealSolution(A=substance_a, B=substance_b)

    fogler_concentrations = np.array([])

    # flux of B on the reactor's outlet

    fb_out = (fogler(0) - fogler(v_pfr)) * f_volumetric

    pfr = rd.idealreactor.StationaryPFR.set_isothermic_isobaric(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=[-1, 1],
        kinetic_argument="concentration",
        reactor_dim_minmax=[0, v_pfr],
        transversal_area=1,
        isothermic_temperature=298.15,
        isobaric_pressure=101325,
        molar_flow_in={"A": fa_in},
        molar_flow_out={"B": fb_out},
    )

    # reactord solution
    solution = pfr.simulate(grid_size=100)
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[0:2, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)


def test_fogler_p1_15a_ivp_outlet():
    """Fogler fourth ed. P1.15a as a initial value problem, but initial
    conditions are in the reactor's outlet"""

    def volume(temperature, pressure):
        return 1 / (fa_in / f_volumetric)

    def kinetic(concentrations, temperature):
        return k

    # Fogler's exact solution

    def fogler(vol):
        concentration = (fa_in / f_volumetric) - vol * k / f_volumetric

        return concentration

    fa_in = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 60  # m3/s

    k = 0.05 / 3600 / 0.001  # mol/s/m3

    v_pfr = 99 * 0.001  # m3

    substance_a = rd.Substance(name="A")
    substance_a.volume_liquid = volume

    substance_b = rd.Substance(name="B")
    substance_b.volume_liquid = volume

    mixture = rd.mix.IdealSolution(A=substance_a, B=substance_b)

    fogler_concentrations = np.array([])

    # flux of A and B on the reactor's outlet

    fa_out = fogler(v_pfr) * f_volumetric
    fb_out = (fogler(0) - fogler(v_pfr)) * f_volumetric

    pfr = rd.idealreactor.StationaryPFR.set_isothermic_isobaric(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=[-1, 1],
        kinetic_argument="concentration",
        reactor_dim_minmax=[0, v_pfr],
        transversal_area=1,
        isothermic_temperature=298.15,
        isobaric_pressure=101325,
        molar_flow_out={"A": fa_out, "B": fb_out},
    )

    # reactord solution
    solution = pfr.simulate(grid_size=100)
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[0:2, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)


def test_fogler_p1_15b():
    """Fogler fourth ed. P1.15b as initial value problem"""

    def volume(temperature, pressure):
        return 1 / (fa_in / f_volumetric)

    def kinetic(concentrations, temperature):
        return k * concentrations[0]

    # Fogler's exact solution

    def fogler(vol):
        concentration = (fa_in / f_volumetric) / np.exp(vol * k / f_volumetric)
        return concentration

    fa_in = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 60  # m3/s

    k = 0.0001  # 1/s

    v_pfr = 127.9 * 0.001  # m3

    substance_a = rd.Substance(name="A")
    substance_a.volume_liquid = volume

    substance_b = rd.Substance(name="B")
    substance_b.volume_liquid = volume

    mixture = rd.mix.IdealSolution(A=substance_a, B=substance_b)

    pfr = rd.idealreactor.StationaryPFR.set_isothermic_isobaric(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=[-1, 1],
        kinetic_argument="concentration",
        reactor_dim_minmax=[0, v_pfr],
        transversal_area=1,
        isothermic_temperature=298.15,
        isobaric_pressure=101325,
        molar_flow_in={"A": fa_in, "B": 0},
    )

    fogler_concentrations = np.array([])

    # reactord solution
    solution = pfr.simulate(grid_size=100)
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[0:2, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)


def test_fogler_p1_15b_bvp():
    """Fogler fourth ed. P1.15b as initial value problem"""

    def volume(temperature, pressure):
        return 1 / (fa_in / f_volumetric)

    def kinetic(concentrations, temperature):
        return k * concentrations[0]

    # Fogler's exact solution

    def fogler(vol):
        concentration = (fa_in / f_volumetric) / np.exp(vol * k / f_volumetric)
        return concentration

    fa_in = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 60  # m3/s

    k = 0.0001  # 1/s

    v_pfr = 127.9 * 0.001  # m3

    substance_a = rd.Substance(name="A")
    substance_a.volume_liquid = volume

    substance_b = rd.Substance(name="B")
    substance_b.volume_liquid = volume

    mixture = rd.mix.IdealSolution(A=substance_a, B=substance_b)

    fogler_concentrations = np.array([])

    # Out molar flux of B to border condition

    fb_out = (fogler(0) - fogler(v_pfr)) * f_volumetric

    pfr = rd.idealreactor.StationaryPFR.set_isothermic_isobaric(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=[-1, 1],
        kinetic_argument="concentration",
        reactor_dim_minmax=[0, v_pfr],
        transversal_area=1,
        isothermic_temperature=298.15,
        isobaric_pressure=101325,
        molar_flow_in={"A": fa_in},
        molar_flow_out={"B": fb_out},
    )

    # reactord solution
    solution = pfr.simulate(grid_size=100, tol=0.00001)
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[0:2, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)


def test_fogler_p1_15b_ivp_outlet():
    """Fogler fourth ed. P1.15b as initial value problem but initial
    conditions are on the reactor's outlet"""

    def volume(temperature, pressure):
        return 1 / (fa_in / f_volumetric)

    def kinetic(concentrations, temperature):
        return k * concentrations[0]

    # Fogler's exact solution

    def fogler(vol):
        concentration = (fa_in / f_volumetric) / np.exp(vol * k / f_volumetric)
        return concentration

    fa_in = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 60  # m3/s

    k = 0.0001  # 1/s

    v_pfr = 127.9 * 0.001  # m3

    substance_a = rd.Substance(name="A")
    substance_a.volume_liquid = volume

    substance_b = rd.Substance(name="B")
    substance_b.volume_liquid = volume

    mixture = rd.mix.IdealSolution(A=substance_a, B=substance_b)

    fogler_concentrations = np.array([])

    # Out molar flux of A and B to border condition

    fa_out = fogler(v_pfr) * f_volumetric
    fb_out = (fogler(0) - fogler(v_pfr)) * f_volumetric

    pfr = rd.idealreactor.StationaryPFR.set_isothermic_isobaric(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=[-1, 1],
        kinetic_argument="concentration",
        reactor_dim_minmax=[0, v_pfr],
        transversal_area=1,
        isothermic_temperature=298.15,
        isobaric_pressure=101325,
        molar_flow_out={"A": fa_out, "B": fb_out},
    )

    # reactord solution
    solution = pfr.simulate(grid_size=100, tol=0.00001)
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[0:2, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)


def test_fogler_p1_15c_ivp():
    """Fogler fourth ed. P1.15c as an initial value problem"""

    def volume(temperature, pressure):
        return 1 / (fa_in / f_volumetric)

    def kinetic(concentrations, temperature):
        return k * concentrations[0] ** 2

    # Fogler's exact solution
    def fogler(vol):
        concentration = 1 / (
            vol * k / f_volumetric + 1 / (fa_in / f_volumetric)
        )
        return concentration

    fa_in = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 60  # m3/s

    k = 3 / 3600 * 0.001  # 1/s

    v_pfr = 660 * 0.001  # m3

    substance_a = rd.Substance(name="A")
    substance_a.volume_liquid = volume

    substance_b = rd.Substance(name="B")
    substance_b.volume_liquid = volume

    mixture = rd.mix.IdealSolution(A=substance_a, B=substance_b)

    pfr = rd.idealreactor.StationaryPFR.set_isothermic_isobaric(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=[-1, 1],
        kinetic_argument="concentration",
        reactor_dim_minmax=[0, v_pfr],
        transversal_area=1,
        isothermic_temperature=298.15,
        isobaric_pressure=101325,
        molar_flow_in={"A": fa_in, "B": 0},
    )

    fogler_concentrations = np.array([])

    # reactord solution
    solution = pfr.simulate(grid_size=100)
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[0:2, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)


def test_fogler_p1_15c_bvp():
    """Fogler fourth ed. P1.15c as a border value problem"""

    def volume(temperature, pressure):
        return 1 / (fa_in / f_volumetric)

    def kinetic(concentrations, temperature):
        return k * concentrations[0] ** 2

    # Fogler's exact solution
    def fogler(vol):
        concentration = 1 / (
            vol * k / f_volumetric + 1 / (fa_in / f_volumetric)
        )
        return concentration

    fa_in = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 60  # m3/s

    k = 3 / 3600 * 0.001  # 1/s

    v_pfr = 660 * 0.001  # m3

    substance_a = rd.Substance(name="A")
    substance_a.volume_liquid = volume

    substance_b = rd.Substance(name="B")
    substance_b.volume_liquid = volume

    mixture = rd.mix.IdealSolution(A=substance_a, B=substance_b)

    fogler_concentrations = np.array([])

    # B molar flux in the reactor's outlet for the border conditions

    fb_out = (fogler(0) - fogler(v_pfr)) * f_volumetric

    pfr = rd.idealreactor.StationaryPFR.set_isothermic_isobaric(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=[-1, 1],
        kinetic_argument="concentration",
        reactor_dim_minmax=[0, v_pfr],
        transversal_area=1,
        isothermic_temperature=298.15,
        isobaric_pressure=101325,
        molar_flow_in={"A": fa_in},
        molar_flow_out={"B": fb_out},
    )

    # reactord solution
    solution = pfr.simulate(grid_size=100, tol=1e-6)
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[0:2, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)


def test_fogler_p1_15c_ivp_outlet():
    """Fogler fourth ed. P1.15c as an initial value problem but the
    initial conditions are on the reactor's outlet"""

    def volume(temperature, pressure):
        return 1 / (fa_in / f_volumetric)

    def kinetic(concentrations, temperature):
        return k * concentrations[0] ** 2

    # Fogler's exact solution
    def fogler(vol):
        concentration = 1 / (
            vol * k / f_volumetric + 1 / (fa_in / f_volumetric)
        )
        return concentration

    fa_in = 5 / 3600  # mol/s

    f_volumetric = 10 * 0.001 / 60  # m3/s

    k = 3 / 3600 * 0.001  # 1/s

    v_pfr = 660 * 0.001  # m3

    substance_a = rd.Substance(name="A")
    substance_a.volume_liquid = volume

    substance_b = rd.Substance(name="B")
    substance_b.volume_liquid = volume

    mixture = rd.mix.IdealSolution(A=substance_a, B=substance_b)

    fogler_concentrations = np.array([])

    # A and B molar flux in the reactor's outlet for the border conditions
    fa_out = fogler(v_pfr) * f_volumetric
    fb_out = (fogler(0) - fogler(v_pfr)) * f_volumetric

    pfr = rd.idealreactor.StationaryPFR.set_isothermic_isobaric(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=[-1, 1],
        kinetic_argument="concentration",
        reactor_dim_minmax=[0, v_pfr],
        transversal_area=1,
        isothermic_temperature=298.15,
        isobaric_pressure=101325,
        molar_flow_out={"A": fa_out, "B": fb_out},
    )

    # reactord solution
    solution = pfr.simulate(grid_size=100, tol=1e-6)
    reactord_concentrations = np.array([])

    # Comparisson

    for i, v in enumerate(solution.x):
        reactord_concentrations = np.append(
            reactord_concentrations,
            mixture.concentrations(solution.y[0:2, i], 298.15, 101325)[0],
        )

        fogler_concentrations = np.append(fogler_concentrations, fogler(v))

    assert np.allclose(reactord_concentrations, fogler_concentrations)


def test_fogler_example_4_4():
    def kinetic(concentrations, temperature):
        return 0

    mixture = rd.mix.IdealGas(n="nitrogen", o="oxygen")

    reactor = rd.idealreactor.StationaryPFR.set_isothermic_noisobaric(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=[0, 0],
        kinetic_argument="partial_pressure",
        reactor_dim_minmax=[0, 18.288],
        transversal_area=0.001313648986,
        isothermic_temperature=260 + 273.15,
        pressure_in_out={"in": 10 * 101325},
        pressure_loss_equation="packed bed reactor",
        packed_bed_porosity=0.45,
        packed_bed_particle_diameter=0.00635,
        molar_flow_in={
            "nitrogen": 0.3560387507669636,
            "oxygen": 0.09983673036871321,
        },
    )

    result = reactor.simulate(
        grid_size=7, tol=0.000001, max_nodes=7, verbose=0
    )

    fogler_pressures = [10, 9.2, 8.3, 7.3, 6.2, 4.7, 2.65]

    assert np.allclose(
        result.y[-2, :] / 101325, fogler_pressures, rtol=1e-05, atol=0.1
    )
