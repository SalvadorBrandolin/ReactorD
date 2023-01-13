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

    pfr = rd.idealreactor.StationaryPFR.from_isothermic_isobaric(
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

    pfr = rd.idealreactor.StationaryPFR.from_isothermic_isobaric(
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

    pfr = rd.idealreactor.StationaryPFR.from_isothermic_isobaric(
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

    pfr = rd.idealreactor.StationaryPFR.from_isothermic_isobaric(
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

    pfr = rd.idealreactor.StationaryPFR.from_isothermic_isobaric(
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

    pfr = rd.idealreactor.StationaryPFR.from_isothermic_isobaric(
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

    pfr = rd.idealreactor.StationaryPFR.from_isothermic_isobaric(
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

    pfr = rd.idealreactor.StationaryPFR.from_isothermic_isobaric(
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

    pfr = rd.idealreactor.StationaryPFR.from_isothermic_isobaric(
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
    # Pressure border condition information is given at the reactor's inlet.
    def kinetic(concentrations, temperature):
        return 0

    mixture = rd.mix.IdealGas(n="nitrogen", o="oxygen")

    reactor = rd.idealreactor.StationaryPFR.from_isothermic_noisobaric(
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

    fogler_z = np.array([0, 10, 20, 30, 40, 50, 60]) * 0.30480370641
    fogler_pressures = np.array([10, 9.2, 8.3, 7.3, 6.2, 4.7, 2.65])
    reactord_pressures = result.sol(fogler_z)[-2] / 101325  # Pa to atm

    assert np.allclose(
        reactord_pressures, fogler_pressures, rtol=1e-05, atol=0.1
    )


def test_fogler_example_4_4_border_condition():
    # Pressure border condition information is given at the reactor's outlet.
    def kinetic(concentrations, temperature):
        return 0

    mixture = rd.mix.IdealGas(n="nitrogen", o="oxygen")

    reactor = rd.idealreactor.StationaryPFR.from_isothermic_noisobaric(
        mix=mixture,
        list_of_reactions=[kinetic],
        stoichiometry=[0, 0],
        kinetic_argument="partial_pressure",
        reactor_dim_minmax=[0, 18.288],
        transversal_area=0.001313648986,
        isothermic_temperature=260 + 273.15,
        pressure_in_out={"out": 2.65 * 101325},
        pressure_loss_equation="packed bed reactor",
        packed_bed_porosity=0.45,
        packed_bed_particle_diameter=0.00635,
        molar_flow_in={
            "nitrogen": 0.3560387507669636,
            "oxygen": 0.09983673036871321,
        },
    )

    result = reactor.simulate(
        grid_size=100, tol=0.0001, max_nodes=1000, verbose=0
    )

    fogler_z = np.array([0, 10, 20, 30, 40, 50, 60]) * 0.30480370641
    fogler_pressures = np.array([10, 9.2, 8.3, 7.3, 6.2, 4.7, 2.65])
    reactord_pressures = result.sol(fogler_z)[-2] / 101325  # Pa to atm

    assert np.allclose(
        reactord_pressures, fogler_pressures, rtol=1e-05, atol=0.1
    )


def test_fogler_example_8_5():
    # 4th edition Fogler example 8-5
    # acetone -> ethenone + methane 

    # Data obtained Graph grabber
    fogler_z_for_temp = np.array(
        [
            0.0111,
            0.0998,
            0.2318,
            0.3896,
            0.6214,
            0.9359,
            1.5254,
            2.8656,
            3.9872,
        ]
    )

    fogler_temp = np.array(
        [
            1005.04,
            976.14,
            955.79,
            941.86,
            930.88,
            924.16,
            916.09,
            905.28,
            901.98,
        ]
    )

    fogler_z_for_conversion = np.array(
        [
            0.0033,
            0.1291,
            0.2857,
            0.5080,
            0.8567,
            1.2188,
            2.0565,
            2.9422,
            3.9763,
        ]
    )

    fogler_conversion = np.array(
        [
            0.0081,
            0.0954,
            0.1434,
            0.1800,
            0.2153,
            0.2305,
            0.2430,
            0.2549,
            0.2667,
        ]
    )

    def reaction_rate(concentrations, temperature):
        k = np.exp(34.34 - 34222 / temperature)
        return k * concentrations[0]

    mix = rd.mix.IdealGas(
        acetone="acetone", anhydride="ethenone", methane="methane"
    )

    pfr = rd.idealreactor.StationaryPFR.from_adiabatic_isobaric(
        mix = mix,
        list_of_reactions=[reaction_rate],
        stoichiometry=[-1, 1, 1],
        kinetic_argument='concentration',
        reactor_dim_minmax=[0, 5/1000],
        transversal_area=1,
        temperature_in_out={'in': 1035},
        isobaric_pressure=162000,
        molar_flow_in={'acetone': 0.0376, 'ethenone': 0, 'methane': 0},
    )

    # pfr.simulate()
