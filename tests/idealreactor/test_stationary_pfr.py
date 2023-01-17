import numpy as np

import reactord as rd

import pytest

# ======================================================================
# Fogler's Elements Of Chemical Reaction Engineering 4th ed.
# p1.15
# ======================================================================

@pytest.mark.sim
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

@pytest.mark.sim
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

@pytest.mark.sim
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

@pytest.mark.sim
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

@pytest.mark.sim
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

@pytest.mark.sim
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

@pytest.mark.sim
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

@pytest.mark.sim
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

@pytest.mark.sim
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

@pytest.mark.sim
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

@pytest.mark.sim
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

@pytest.mark.sim
@pytest.mark.slow
@pytest.mark.skip("need debug")
def test_fogler_example_8_5():
    # 4th edition Fogler example 8-5
    # acetone -> ethenone + methane

    def reaction_rate(concentrations, temperature):
        k = np.exp(34.34 - 34222 / temperature)
        return k * concentrations[0]

    def acetone_cp(temperature, pressure):
        return 163

    def ethenone_cp(temperature, pressure):
        return 83

    def methane_cp(temperature, pressure):
        return 71

    # Data obtained WebPlotDigitizer 4.6
    fogler_z_for_temp = np.array(
        [
            0.031,
            0.148,
            0.309,
            0.768,
            1.502,
            2.446,
            3.288,
            3.949,
        ]
    )

    fogler_temp = np.array(
        [
            1018.128,
            983.496,
            960.397,
            933.96,
            921.222,
            910.653,
            903.944,
            902.205,
        ]
    )

    fogler_z_for_conversion = np.array(
        [
            0.057,
            0.208,
            0.518,
            0.989,
            1.751,
            2.803,
            3.587,
            3.964,
        ]
    )
    fogler_conversion = np.array(
        [
            0.06,
            0.123,
            0.184,
            0.224,
            0.241,
            0.255,
            0.264,
            0.269,
        ]
    )

    acetone = rd.Substance(name="acetone", heat_capacity_gas=acetone_cp)
    ethenone = rd.Substance(name="ethenone", heat_capacity_gas=ethenone_cp)
    methane = rd.Substance(name="methane", heat_capacity_gas=methane_cp)

    mix = rd.mix.IdealGas(acetone=acetone, ethenone=ethenone, methane=methane)

    pfr = rd.idealreactor.StationaryPFR.from_adiabatic_isobaric(
        mix=mix,
        list_of_reactions=[reaction_rate],
        stoichiometry=[-1, 1, 1],
        kinetic_argument="concentration",
        reactor_dim_minmax=[0, 5],
        transversal_area=1,
        temperature_in_out={"in": 1035},
        isobaric_pressure=162000,
        molar_flow_in={"acetone": 37.6, "ethenone": 0, "methane": 0},
        reaction_enthalpies=[80.77 * 1000],
    )

    result = pfr.simulate()

    conversion_reactord = (
        37.6 - result.sol(fogler_z_for_conversion)[0]
    ) / 37.6

    print("reactord")
    print(result.sol(fogler_z_for_temp)[-3])
    print(conversion_reactord)
    print("fogler")
    print(fogler_temp)

    assert np.allclose(result.sol(fogler_z_for_temp)[-3], fogler_temp, atol=4)

    conversion_reactord = (37.6 - result.sol(fogler_z_for_conversion)) / 37.6
    assert np.allclose(conversion_reactord, fogler_conversion)

@pytest.mark.sim
@pytest.mark.slow
def test_fogler_example_11_3():
    """Fogler 5th edition example 11-3

    http://websites.umich.edu/~elements/5e/11chap/live.html
    """

    def cp_butane(temperature, pressure):
        return 141

    def cp_isobutane(temperature, pressure):
        return 141

    def cp_isopentane(emperature, pressure):
        return 161

    def kinetic_law(concentrations, temperature):
        k_360k = 31.1 / 3600  # 1 / s
        e = 65.7 * 1000  # J/mol
        r = 8.31446261815324

        k = k_360k * np.exp(e / r * (1 / 360 - 1 / temperature))

        kc = 3.03 * np.exp(-6900 / r * (1 / (60 + 273.15) - 1 / temperature))

        r = k * (concentrations[0] - concentrations[1] / kc)

        return r

    def volume(temperature, pressure):
        return 1 / 9300  # m3/mol

    butane = rd.Substance(
        name="butane",
        heat_capacity_liquid=cp_butane,
        volume_liquid=volume,
    )

    isobutane = rd.Substance(
        name="isobutane",
        heat_capacity_liquid=cp_isobutane,
        volume_liquid=volume,
    )

    isopentane = rd.Substance(
        name="isopentane",
        heat_capacity_liquid=cp_isopentane,
        volume_liquid=volume,
    )

    mix = rd.mix.IdealSolution(a=butane, b=isobutane, c=isopentane)

    ft0 = 163 * 1000 / 3600

    pfr = rd.idealreactor.StationaryPFR.from_adiabatic_isobaric(
        mix=mix,
        list_of_reactions=[kinetic_law],
        stoichiometry=[-1, 1, 0],
        kinetic_argument="concentration",
        reactor_dim_minmax=[0, 5],
        transversal_area=1,
        temperature_in_out={"in": 330},
        isobaric_pressure=101325,
        molar_flow_in={
            "butane": ft0 * 0.9,
            "isobutane": 0,
            "isopentane": ft0 * 0.1,
        },
        reaction_enthalpies=[-6900],
    )

    result = pfr.simulate()

    v = np.linspace(0, 5, 100)

    fogler_conversion = np.array(
        [
            0.0,
            0.01385284,
            0.02805776,
            0.04262434,
            0.05756174,
            0.07287853,
            0.08858242,
            0.10468003,
            0.12117659,
            0.13807556,
            0.15537822,
            0.17308316,
            0.19118575,
            0.20967752,
            0.22854549,
            0.24777135,
            0.26733079,
            0.2871926,
            0.30731791,
            0.32765948,
            0.34816106,
            0.36875705,
            0.38937237,
            0.40992273,
            0.43031542,
            0.45045061,
            0.47022326,
            0.48952566,
            0.50825043,
            0.52629397,
            0.54356002,
            0.5599632,
            0.5754321,
            0.58991177,
            0.60336534,
            0.61577461,
            0.62713971,
            0.63747776,
            0.64682083,
            0.65521344,
            0.66270969,
            0.66937046,
            0.67526071,
            0.68044711,
            0.68499606,
            0.68897207,
            0.6924366,
            0.69544722,
            0.69805711,
            0.70031487,
            0.70226441,
            0.70394513,
            0.70539208,
            0.70663629,
            0.70770503,
            0.70862225,
            0.70940881,
            0.71008288,
            0.71066022,
            0.71115446,
            0.71157739,
            0.71193916,
            0.71224853,
            0.71251301,
            0.71273907,
            0.71293225,
            0.71309731,
            0.71323832,
            0.71335876,
            0.71346164,
            0.71354949,
            0.71362452,
            0.71368858,
            0.71374329,
            0.71378999,
            0.71382987,
            0.71386391,
            0.71389298,
            0.71391779,
            0.71393897,
            0.71395706,
            0.71397249,
            0.71398567,
            0.71399692,
            0.71400652,
            0.71401472,
            0.71402172,
            0.71402769,
            0.71403279,
            0.71403714,
            0.71404086,
            0.71404403,
            0.71404674,
            0.71404905,
            0.71405102,
            0.7140527,
            0.71405414,
            0.71405537,
            0.71405641,
            0.71405731,
        ]
    )

    fogler_temperature = np.array(
        [
            330.0,
            330.60158142,
            331.21845218,
            331.85102888,
            332.49970931,
            333.16486495,
            333.84683085,
            334.54589492,
            335.26228396,
            335.99614833,
            336.74754368,
            337.51640839,
            338.30254191,
            339.10557642,
            339.92494732,
            340.75986091,
            341.60926024,
            342.47179053,
            343.34576378,
            344.22912826,
            345.11944168,
            346.01385529,
            346.90910772,
            347.8015397,
            348.68712435,
            349.5615265,
            350.42018506,
            351.25842203,
            352.07157461,
            352.85514363,
            353.60494921,
            354.31728304,
            354.98904421,
            355.61784686,
            356.2020893,
            356.74098151,
            357.23452873,
            357.6834746,
            358.08921233,
            358.45367451,
            358.77921106,
            359.0684656,
            359.32425886,
            359.54948647,
            359.74703178,
            359.91969639,
            360.07014901,
            360.20088966,
            360.31422842,
            360.41227506,
            360.49693712,
            360.56992498,
            360.63276112,
            360.68679253,
            360.73320465,
            360.77303622,
            360.80719386,
            360.83646635,
            360.86153803,
            360.88300123,
            360.90136757,
            360.91707822,
            360.93051304,
            360.94199866,
            360.95181568,
            360.96020486,
            360.96737268,
            360.97349611,
            360.97872668,
            360.98319411,
            360.98700941,
            360.99026753,
            360.99304966,
            360.9954252,
            360.99745349,
            360.99918524,
            361.00066367,
            361.00192583,
            361.00300335,
            361.00392323,
            361.00470853,
            361.00537892,
            361.00595119,
            361.00643971,
            361.00685672,
            361.00721269,
            361.00751651,
            361.00777591,
            361.00799734,
            361.00818637,
            361.00834772,
            361.00848544,
            361.00860299,
            361.00870333,
            361.00878899,
            361.0088621,
            361.00892452,
            361.00897779,
            361.00902326,
            361.00906205,
        ]
    )

    reactord_temperature = result.sol(v)[-3]

    reactord_conversion = (ft0 * 0.9 - result.sol(v)[0]) / (ft0 * 0.9)

    assert np.allclose(
        fogler_temperature, reactord_temperature, atol=0.1, rtol=0.1
    )

    assert np.allclose(
        fogler_conversion, reactord_conversion, atol=0.05, rtol=0.1
    )
