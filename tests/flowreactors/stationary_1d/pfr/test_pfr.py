import numpy as np

import reactord as rd
import reactord.flowreactors.stationary_1d.pfr as pfr

from scipy.constants import R


def test_repr():
    def r_rate():
        ...

    a = rd.Substance("A")
    b = rd.Substance("B")
    c = rd.Substance("C")

    mix = rd.mix.IdealGas([a, b, c])

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={"r1": {"eq": a + b > c, "rate": r_rate}},
        kinetic_constants={},
    )

    mb = pfr.mass_balances.MolarFlow(molar_flows_in={"A": 10, "B": 10, "C": 0})

    eb = pfr.energy_balances.Adiabatic(temperature_in_or_out={"in": 303.15})

    pb = pfr.pressure_balances.Isobaric(pressure=101325)

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=10,
        transversal_area=1,
        grid_size=10,
        mass_balance=mb,
        energy_balance=eb,
        pressure_balance=pb,
    )

    text = (
        f"{kinetic.__repr__()}\n"
        f"{mb.__repr__()}\n"
        f"{eb.__repr__()[0]}\n"
        f"{eb.__repr__()[1]}\n"
        f"{pb.__repr__()}\n"
    )

    assert text == reactor.__repr__()


def test_fogler_p1_15a():
    """Fogler fourth ed. P1.15a."""

    def volume(temperature, pressure):
        n = np.size(temperature)
        return np.full(n, 1 / (fa_in / f_volumetric))

    def r_rate(c, t, cons):
        return np.full(np.size(t), cons["k"])

    # Fogler's exact solution
    def fogler(vol):
        concentration = (fa_in / f_volumetric) - vol * k / f_volumetric
        return concentration

    # =========================================================================
    # Initial value problem
    # =========================================================================
    # Problem data
    fa_in = 5 / 3600  # mol/s
    f_volumetric = 10 * 0.001 / 60  # m3/s
    k = 0.05 / 3600 / 0.001  # mol/s/m3
    v_pfr = 99 * 0.001  # m3

    # Substance definition
    a = rd.Substance(
        name="A",
        volume_liquid=volume,
    )
    b = rd.Substance(
        name="B",
        volume_liquid=volume,
    )

    mix = rd.mix.IdealSolution([a, b])

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={"r1": {"eq": a > b, "rate": r_rate}},
        kinetic_constants={"k": k},
    )

    mb = pfr.mass_balances.MolarFlow(molar_flows_in={"A": fa_in, "B": 0})
    eb = pfr.energy_balances.Isothermic(298.15)
    pb = pfr.pressure_balances.Isobaric(101325)

    reactor = pfr.PFR(kinetic, v_pfr, 1, 100, mb, eb, pb)

    # Simulation
    reactor.simulate(verbose=0)

    # Comparisson
    fogler_a_concentration = fogler(reactor.z)
    fogler_b_concentration = fogler_a_concentration[0] - fogler_a_concentration

    reactord_concentrations = reactor.mix.concentrations(
        reactor.mole_fraction_profile,
        reactor.temperature_profile,
        reactor.pressure_profile,
    )

    assert np.allclose(
        reactord_concentrations[0, :], fogler_a_concentration, atol=1e-8
    )
    assert np.allclose(
        reactord_concentrations[1, :], fogler_b_concentration, atol=1e-8
    )

    # =========================================================================
    # Border value problem 1
    # =========================================================================
    fa_out = fogler_a_concentration[-1] * f_volumetric
    fb_out = fogler_b_concentration[-1] * f_volumetric

    mb1 = pfr.mass_balances.MolarFlow(
        molar_flows_in={"A": fa_in}, molar_flows_out={"B": fb_out}
    )
    eb1 = pfr.energy_balances.Isothermic(298.15)
    pb1 = pfr.pressure_balances.Isobaric(101325)

    reactor1 = pfr.PFR(
        kinetic=kinetic,
        reactor_length=v_pfr,
        transversal_area=1,
        grid_size=100,
        mass_balance=mb1,
        energy_balance=eb1,
        pressure_balance=pb1,
    )

    # Simulation
    reactor1.simulate(tol=0.0001, verbose=0, bc_tol=1e-3)

    # Comparisson
    fogler_a_concentration1 = fogler(reactor1.z)
    fogler_b_concentration1 = (
        fogler_a_concentration1[0] - fogler_a_concentration1
    )

    reactord_concentrations1 = reactor1.mix.concentrations(
        reactor1.mole_fraction_profile,
        reactor1.temperature_profile,
        reactor1.pressure_profile,
    )

    assert np.allclose(
        reactord_concentrations1[0, :], fogler_a_concentration1, atol=1e-8
    )
    assert np.allclose(
        reactord_concentrations1[1, :], fogler_b_concentration1, atol=1e-8
    )

    # =========================================================================
    # Border value problem 2
    # =========================================================================
    mb2 = pfr.mass_balances.MolarFlow(
        molar_flows_out={"A": fa_out, "B": fb_out}
    )
    eb2 = pfr.energy_balances.Isothermic(298.15)
    pb2 = pfr.pressure_balances.Isobaric(101325)

    reactor2 = pfr.PFR(
        kinetic=kinetic,
        reactor_length=v_pfr,
        transversal_area=1,
        grid_size=100,
        mass_balance=mb2,
        energy_balance=eb2,
        pressure_balance=pb2,
    )

    # Simulation
    reactor2.simulate(tol=0.0001, verbose=0, bc_tol=1e-3)

    # Comparisson
    fogler_a_concentration2 = fogler(reactor2.z)
    fogler_b_concentration2 = (
        fogler_a_concentration2[0] - fogler_a_concentration2
    )

    reactord_concentrations2 = reactor2.mix.concentrations(
        reactor2.mole_fraction_profile,
        reactor2.temperature_profile,
        reactor2.pressure_profile,
    )

    assert np.allclose(
        reactord_concentrations2[0, :], fogler_a_concentration2, atol=1e-8
    )
    assert np.allclose(
        reactord_concentrations2[1, :], fogler_b_concentration2, atol=1e-8
    )


def test_fogler_p1_15b():
    """Fogler fourth ed. P1.15b."""

    def volume(temperature, pressure):
        n = np.size(temperature)
        return np.full(n, 1 / (fa_in / f_volumetric))

    def r_rate(c, t, cons):
        return cons["k"] * c["A"]

    # Fogler's exact solution
    def fogler(vol):
        concentration = (fa_in / f_volumetric) / np.exp(vol * k / f_volumetric)
        return concentration

    # =========================================================================
    # Initial value problem
    # =========================================================================
    # Problem data
    fa_in = 5 / 3600  # mol/s
    f_volumetric = 10 * 0.001 / 60  # m3/s
    k = 0.0001  # 1/s
    v_pfr = 127.9 * 0.001  # m3

    # Substance definition
    a = rd.Substance(
        name="A",
        volume_liquid=volume,
    )
    b = rd.Substance(
        name="B",
        volume_liquid=volume,
    )

    # Mixture
    mixture = rd.mix.IdealSolution([a, b])

    # Kinetic
    kinetic = rd.Kinetic(
        mix=mixture,
        reactions={"r1": {"eq": a > b, "rate": r_rate}},
        kinetic_constants={"k": k},
    )

    # Reactor
    mb = pfr.mass_balances.MolarFlow(molar_flows_in={"A": fa_in, "B": 0})
    eb = pfr.energy_balances.Isothermic(298.15)
    pb = pfr.pressure_balances.Isobaric(101325)

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=v_pfr,
        transversal_area=1,
        grid_size=20,
        mass_balance=mb,
        energy_balance=eb,
        pressure_balance=pb,
    )

    # Simulation
    reactor.simulate(tol=1e-8, verbose=0)

    # Comparisson
    fogler_a_concentration = fogler(reactor.z)
    fogler_b_concentration = fogler_a_concentration[0] - fogler_a_concentration

    reactord_concentrations = reactor.mix.concentrations(
        reactor.mole_fraction_profile,
        reactor.temperature_profile,
        reactor.pressure_profile,
    )

    assert np.allclose(
        reactord_concentrations[0, :], fogler_a_concentration, atol=1e-8
    )
    assert np.allclose(
        reactord_concentrations[1, :], fogler_b_concentration, atol=1e-8
    )

    # =========================================================================
    # Border value problem 1
    # =========================================================================
    fa_out = fogler_a_concentration[-1] * f_volumetric
    fb_out = fogler_b_concentration[-1] * f_volumetric

    mb1 = pfr.mass_balances.MolarFlow(
        molar_flows_in={"A": fa_in}, molar_flows_out={"B": fb_out}
    )
    eb1 = pfr.energy_balances.Isothermic(298.15)
    pb1 = pfr.pressure_balances.Isobaric(101325)

    reactor1 = pfr.PFR(
        kinetic=kinetic,
        reactor_length=v_pfr,
        transversal_area=1,
        grid_size=100,
        mass_balance=mb1,
        energy_balance=eb1,
        pressure_balance=pb1,
    )

    # Simulation
    reactor1.simulate(tol=1e-10, verbose=0)

    # Comparisson
    fogler_a_concentration1 = fogler(reactor1.z)
    fogler_b_concentration1 = (
        fogler_a_concentration1[0] - fogler_a_concentration1
    )

    reactord_concentrations1 = reactor1.mix.concentrations(
        reactor1.mole_fraction_profile,
        reactor1.temperature_profile,
        reactor1.pressure_profile,
    )

    assert np.allclose(
        reactord_concentrations1[0, :], fogler_a_concentration1, atol=1e-8
    )
    assert np.allclose(
        reactord_concentrations1[1, :], fogler_b_concentration1, atol=1e-8
    )

    # =========================================================================
    # Border value problem 2
    # =========================================================================
    mb2 = pfr.mass_balances.MolarFlow(
        molar_flows_out={"A": fa_out, "B": fb_out}
    )
    eb2 = pfr.energy_balances.Isothermic(298.15)
    pb2 = pfr.pressure_balances.Isobaric(101325)

    reactor2 = pfr.PFR(
        kinetic=kinetic,
        reactor_length=v_pfr,
        transversal_area=1,
        grid_size=100,
        mass_balance=mb2,
        energy_balance=eb2,
        pressure_balance=pb2,
    )

    # Simulation
    reactor2.simulate(tol=1e-8, verbose=0, bc_tol=1e-3)

    # Comparisson
    fogler_a_concentration2 = fogler(reactor2.z)
    fogler_b_concentration2 = (fa_in / f_volumetric) - fogler_a_concentration2

    reactord_concentrations2 = reactor2.mix.concentrations(
        reactor2.mole_fraction_profile,
        reactor2.temperature_profile,
        reactor2.pressure_profile,
    )

    assert np.allclose(
        reactord_concentrations2[0, :], fogler_a_concentration2, atol=1e-8
    )
    assert np.allclose(
        reactord_concentrations2[1, :], fogler_b_concentration2, atol=1e-8
    )


def test_fogler_p1_15c():
    """Fogler fourth ed. P1.15c."""

    def volume(temperature, pressure):
        n = np.size(temperature)
        return np.full(n, 1 / (fa_in / f_volumetric))

    def r_rate(c, t, cons):
        return cons["k"] * c["A"] ** 2

    # Fogler's exact solution
    def fogler(vol):
        concentration = 1 / (
            vol * k / f_volumetric + 1 / (fa_in / f_volumetric)
        )
        return concentration

    # =========================================================================
    # Initial value problem
    # =========================================================================
    # Problem data
    fa_in = 5 / 3600  # mol/s
    f_volumetric = 10 * 0.001 / 60  # m3/s
    k = 3 / 3600 * 0.001  # 1/s
    v_pfr = 660 * 0.001  # m3

    # Substance definition
    a = rd.Substance(
        name="A",
        volume_liquid=volume,
    )
    b = rd.Substance(
        name="B",
        volume_liquid=volume,
    )

    # Mixture
    mix = rd.mix.IdealSolution([a, b])

    # Kinetic
    kinetic = rd.Kinetic(
        mix=mix,
        reactions={"r1": {"eq": a > b, "rate": r_rate}},
        kinetic_constants={"k": k},
    )

    # Reactor
    mb = pfr.mass_balances.MolarFlow(molar_flows_in={"A": fa_in, "B": 0})
    eb = pfr.energy_balances.Isothermic(298.15)
    pb = pfr.pressure_balances.Isobaric(101325)

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=v_pfr,
        transversal_area=1,
        grid_size=20,
        mass_balance=mb,
        energy_balance=eb,
        pressure_balance=pb,
    )

    # Simulation
    reactor.simulate(tol=1e-8, verbose=0)

    # Comparisson
    fogler_a_concentration = fogler(reactor.z)
    fogler_b_concentration = fogler_a_concentration[0] - fogler_a_concentration

    reactord_concentrations = reactor.mix.concentrations(
        reactor.mole_fraction_profile,
        reactor.temperature_profile,
        reactor.pressure_profile,
    )

    assert np.allclose(
        reactord_concentrations[0, :], fogler_a_concentration, atol=1e-8
    )
    assert np.allclose(
        reactord_concentrations[1, :], fogler_b_concentration, atol=1e-8
    )

    # =========================================================================
    # Border value problem 1
    # =========================================================================
    fa_out = fogler_a_concentration[-1] * f_volumetric
    fb_out = fogler_b_concentration[-1] * f_volumetric

    mb1 = pfr.mass_balances.MolarFlow(
        molar_flows_in={"A": fa_in}, molar_flows_out={"B": fb_out}
    )
    eb1 = pfr.energy_balances.Isothermic(298.15)
    pb1 = pfr.pressure_balances.Isobaric(101325)

    reactor1 = pfr.PFR(
        kinetic=kinetic,
        reactor_length=v_pfr,
        transversal_area=1,
        grid_size=20,
        mass_balance=mb1,
        energy_balance=eb1,
        pressure_balance=pb1,
    )

    # Simulation
    reactor1.simulate(
        tol=1e-10,
        verbose=0,
    )

    # Comparisson
    fogler_a_concentration1 = fogler(reactor1.z)
    fogler_b_concentration1 = (
        fogler_a_concentration1[0] - fogler_a_concentration1
    )

    reactord_concentrations1 = reactor1.mix.concentrations(
        reactor1.mole_fraction_profile,
        reactor1.temperature_profile,
        reactor1.pressure_profile,
    )

    assert np.allclose(
        reactord_concentrations1[0, :], fogler_a_concentration1, atol=1e-8
    )
    assert np.allclose(
        reactord_concentrations1[1, :], fogler_b_concentration1, atol=1e-8
    )

    # =========================================================================
    # Border value problem 2
    # =========================================================================
    mb2 = pfr.mass_balances.MolarFlow(
        molar_flows_out={"A": fa_out, "B": fb_out}
    )
    eb2 = pfr.energy_balances.Isothermic(298.15)
    pb2 = pfr.pressure_balances.Isobaric(101325)

    reactor2 = pfr.PFR(
        kinetic=kinetic,
        reactor_length=v_pfr,
        transversal_area=1,
        grid_size=20,
        mass_balance=mb2,
        energy_balance=eb2,
        pressure_balance=pb2,
    )

    # Simulation
    reactor2.simulate(tol=1e-10, verbose=0, bc_tol=1e-3)

    # Comparisson
    fogler_a_concentration2 = fogler(reactor2.z)
    fogler_b_concentration2 = (fa_in / f_volumetric) - fogler_a_concentration2

    reactord_concentrations2 = reactor2.mix.concentrations(
        reactor2.mole_fraction_profile,
        reactor2.temperature_profile,
        reactor2.pressure_profile,
    )

    assert np.allclose(
        reactord_concentrations2[0, :], fogler_a_concentration2, atol=1e-8
    )
    assert np.allclose(
        reactord_concentrations2[1, :], fogler_b_concentration2, atol=1e-8
    )


def test_fogler_example_4_4():
    # Pressure border condition information is given at the reactor's inlet.
    def ra(c, t, cons):
        return np.zeros(np.size(t))

    n = rd.Substance.from_thermo_database("N", "nitrogen")
    o = rd.Substance.from_thermo_database("O", "oxygen")

    mix = rd.mix.IdealGas([n, o])

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={"r1": {"eq": n + o, "rate": ra}},
        kinetic_constants={},
    )

    mb = pfr.mass_balances.MolarFlow(
        molar_flows_in={
            "N": 0.3560387507669636,
            "O": 0.09983673036871321,
        }
    )
    eb = pfr.energy_balances.Isothermic(260 + 273.15)
    pb = pfr.pressure_balances.Ergun(
        pressure={"in": 10 * 101325}, porosity=0.45, particle_diameter=0.00635
    )

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=18.288,
        transversal_area=0.001313648986,
        grid_size=20,
        mass_balance=mb,
        energy_balance=eb,
        pressure_balance=pb,
    )

    reactor.simulate(tol=0.001)

    fogler_z = np.array([0, 10, 20, 30, 40, 50, 60]) * 0.30480370641
    fogler_pressures = np.array([10, 9.2, 8.3, 7.3, 6.2, 4.7, 2.65])
    reactord_pressures = reactor.ode_solution.sol(fogler_z)[-1] / 101325

    assert np.allclose(reactord_pressures, fogler_pressures, atol=0.1)

    # Pressure border condition information is given at the reactor's outlet.
    pb_out = pfr.pressure_balances.Ergun(
        pressure={"out": 2.65 * 101325},
        porosity=0.45,
        particle_diameter=0.00635,
    )

    reactor2 = pfr.PFR(
        kinetic=kinetic,
        reactor_length=18.288,
        transversal_area=0.001313648986,
        grid_size=20,
        mass_balance=mb,
        energy_balance=eb,
        pressure_balance=pb_out,
    )

    reactor2.simulate(tol=0.001)

    fogler_z = np.array([0, 10, 20, 30, 40, 50, 60]) * 0.30480370641
    fogler_pressures = np.array([10, 9.2, 8.3, 7.3, 6.2, 4.7, 2.65])
    reactord_pressures = reactor2.ode_solution.sol(fogler_z)[-1] / 101325

    assert np.allclose(reactord_pressures, fogler_pressures, atol=0.1)


def test_fogler_example_11_3():
    """Fogler 6th ed. example 11.3"""

    # temperature as initial value
    def volume(temperature, pressure):
        f_mol = 163 * 1000 / 3600  # mol / s
        f_volumetric = 100_000 * 0.00378541 / 24 / 3600  # m3 / s

        rho = f_mol / f_volumetric
        v = 1 / rho  # m3 / mol
        return np.full(np.size(temperature), v)

    def cp_butane(temperature, pressure):
        return np.full(np.size(temperature), 141)  # J / mol / K

    def cp_pentane(temperature, pressure):
        return np.full(np.size(temperature), 161)  # J / mol / K

    def r_rate(c, t, cons):
        k360, e, keq60, dh = cons["k360"], cons["e"], cons["keq60"], cons["dh"]

        kd_t = k360 * np.exp(e / R * (1 / 360 - 1 / t))
        keq_t = keq60 * np.exp(dh / R * (1 / (60 + 273.15) - 1 / t))

        rd = kd_t * c["but"]
        ri = kd_t / keq_t * c["i-but"]

        return rd - ri

    dh = -6900  # J / mol

    but = rd.Substance(
        "but", volume_liquid=volume, heat_capacity_liquid=cp_butane
    )

    ibut = rd.Substance(
        "i-but", volume_liquid=volume, heat_capacity_liquid=cp_butane
    )

    ipen = rd.Substance(
        "i-pen", volume_liquid=volume, heat_capacity_liquid=cp_pentane
    )

    mix = rd.mix.IdealSolution([but, ibut, ipen])

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={"r1": {"eq": but > ibut, "rate": r_rate, "DH": dh}},
        kinetic_constants={
            "k360": 31.1 / 3600,
            "e": 65.7 * 1000,
            "keq60": 3.03,
            "dh": dh,
        },
        rates_argument="concentration",
    )

    f_mol = 163 * 1000 / 3600  # mol / s

    mb = pfr.mass_balances.MolarFlow(
        molar_flows_in={"but": f_mol * 0.9, "i-but": 0, "i-pen": f_mol * 0.1},
    )
    eb = pfr.energy_balances.Adiabatic(temperature_in_or_out={"in": 330})
    pb = pfr.pressure_balances.Isobaric(101325)

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=5,
        transversal_area=1,
        grid_size=100,
        mass_balance=mb,
        energy_balance=eb,
        pressure_balance=pb,
    )

    reactor.simulate(1e-10)

    # Fogler's data WebPlotDigitizer v4.6
    vol_temps = np.array(
        [
            0.11963,
            0.30367,
            0.49701,
            0.6535,
            0.98498,
            1.20608,
            1.42704,
            1.64792,
            1.93304,
            2.54816,
            3.07124,
            3.81441,
            4.47502,
            4.97954,
        ]
    )
    temps = np.array(
        [
            331.69216,
            334.0963,
            337.03582,
            339.35166,
            344.60768,
            348.61764,
            352.0026,
            355.03043,
            358.14525,
            360.08756,
            360.69388,
            360.9352,
            361.17946,
            360.89359,
        ]
    )

    vol_x = np.array(
        [
            0.09158,
            0.32051,
            0.57692,
            0.78755,
            0.99817,
            1.26374,
            1.52015,
            1.77656,
            2.28938,
            2.71978,
            3.26923,
        ]
    )
    x = np.array(
        [
            0.02487,
            0.09767,
            0.18291,
            0.25216,
            0.32674,
            0.43504,
            0.53092,
            0.60551,
            0.68551,
            0.70163,
            0.70892,
        ]
    )

    rd_temps = reactor.ode_solution.sol(vol_temps)[-2]
    rd_fbut = reactor.ode_solution.sol(vol_x)[0]
    rd_x = (f_mol * 0.9 - rd_fbut) / (f_mol * 0.9)

    assert np.allclose(temps, rd_temps, atol=0.5)
    assert np.allclose(x, rd_x, atol=0.013)

    # with temperature as border value
    eb = pfr.energy_balances.Adiabatic(
        temperature_in_or_out={"out": 360.89359}
    )

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=4.97954,
        transversal_area=1,
        grid_size=100,
        mass_balance=mb,
        energy_balance=eb,
        pressure_balance=pb,
    )

    reactor.simulate(1e-10)

    rd_temps = reactor.ode_solution.sol(vol_temps)[-2]
    rd_fbut = reactor.ode_solution.sol(vol_x)[0]
    rd_x = (f_mol * 0.9 - rd_fbut) / (f_mol * 0.9)

    assert np.allclose(temps, rd_temps, atol=0.7)
    assert np.allclose(x, rd_x, atol=0.013)


def test_fogler_example_12_2_case2():
    """Fogler 6th ed. example 12.2"""

    # temperature as in value
    def cpa(t, p):
        return np.full(len(t), 163)

    def cpb(t, p):
        return np.full(len(t), 83)

    def cpc(t, p):
        return np.full(len(t), 71)

    def int_cpa(t1, t2, p):
        return 163 * (t2 - t1)

    def int_cpb(t1, t2, p):
        return 83 * (t2 - t1)

    def int_cpc(t1, t2, p):
        return 71 * (t2 - t1)

    def r_rate(c, t, cons):
        k = 8.2e14 * np.exp(-34222 / t)  # 1/s

        return k * c["acetone"]

    a = rd.Substance(
        "acetone",
        formation_enthalpy_ig=-216.67 * 1000,
        heat_capacity_gas=cpa,
        heat_capacity_gas_dt_integral=int_cpa,
    )

    b = rd.Substance(
        "anhydride",
        formation_enthalpy_ig=-61.09 * 1000,
        heat_capacity_gas=cpb,
        heat_capacity_gas_dt_integral=int_cpb,
    )

    c = rd.Substance(
        "methane",
        formation_enthalpy_ig=-74.81 * 1000,
        heat_capacity_gas=cpc,
        heat_capacity_gas_dt_integral=int_cpc,
    )

    mix = rd.mix.IdealGas([a, b, c])

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={"r1": {"eq": a > b + c, "rate": r_rate}},
        kinetic_constants={},
        rates_argument="concentration",
    )

    mb = pfr.mass_balances.MolarFlow(
        molar_flows_in={"acetone": 0.0376, "anhydride": 0, "methane": 0},
    )
    eb = pfr.energy_balances.NoIsothermicAllConstant(
        temperature_in_or_out={"in": 1035},
        refrigerant_in_temperature=1250,
        heat_exchange_coefficient=110,
    )
    pb = pfr.pressure_balances.Isobaric(162 * 1000)

    area = np.pi * (4 / (16500 / 110)) ** 2 / 4

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=1 / 1000 / area,
        transversal_area=area,
        grid_size=100,
        mass_balance=mb,
        energy_balance=eb,
        pressure_balance=pb,
    )

    reactor.simulate(1e-8, bc_tol=1e-5, max_nodes=10_000)

    # Fogler data
    t_z = (
        np.array(
            [
                0.0159,
                0.0398,
                0.0677,
                0.0956,
                0.1275,
                0.1753,
                0.2191,
                0.259,
                0.3028,
                0.3506,
                0.3904,
                0.4382,
                0.4821,
                0.5498,
                0.6215,
                0.6773,
                0.749,
                0.8287,
                0.8924,
                0.9522,
                0.996,
            ]
        )
        / 1000
        / area
    )
    t = np.array(
        [
            1033.93,
            1030.44,
            1030.44,
            1031.31,
            1032.62,
            1035.24,
            1037.86,
            1040.49,
            1043.54,
            1047.04,
            1050.10,
            1053.59,
            1057.09,
            1062.33,
            1068.45,
            1073.25,
            1080.24,
            1088.98,
            1097.28,
            1106.46,
            1113.88,
        ]
    )

    x_z = (
        np.array(
            [
                0.0121,
                0.0684,
                0.1206,
                0.193,
                0.2814,
                0.3537,
                0.426,
                0.4983,
                0.5705,
                0.6428,
                0.719,
                0.8073,
                0.9036,
                0.9958,
            ]
        )
        / 1000
        / area
    )
    x = np.array(
        [
            0.0194,
            0.0922,
            0.1553,
            0.2427,
            0.3398,
            0.4223,
            0.4951,
            0.568,
            0.6311,
            0.6942,
            0.7524,
            0.8204,
            0.8883,
            0.9466,
        ]
    )

    rd_temps = reactor.ode_solution.sol(t_z)[-3]
    rd_facetone = reactor.ode_solution.sol(x_z)[0]
    rd_x = (0.0376 - rd_facetone) / (0.0376)

    assert np.allclose(t, rd_temps, atol=2)
    assert np.allclose(x, rd_x, atol=0.02)

    # temperature as out value
    eb2 = pfr.energy_balances.NoIsothermicAllConstant(
        temperature_in_or_out={"out": 1114.093},
        refrigerant_in_temperature=1250,
        heat_exchange_coefficient=110,
    )

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=1 / 1000 / area,
        transversal_area=area,
        grid_size=100,
        mass_balance=mb,
        energy_balance=eb2,
        pressure_balance=pb,
    )

    reactor.simulate(1e-8, bc_tol=1e-5, max_nodes=10_000)

    rd_temps = reactor.ode_solution.sol(t_z)[-3]
    rd_facetone = reactor.ode_solution.sol(x_z)[0]
    rd_x = (0.0376 - rd_facetone) / (0.0376)

    assert np.allclose(t, rd_temps, atol=2)
    assert np.allclose(x, rd_x, atol=0.02)


def test_fogler_example_12_2_case3():
    """Fogler 6th ed. example 12.2"""

    # temperature as in value
    def cpa(t, p):
        return np.full(len(t), 163)

    def cpb(t, p):
        return np.full(len(t), 83)

    def cpc(t, p):
        return np.full(len(t), 71)

    def cp_ref(t, p):
        return np.full(len(t), 34.5)

    def int_cpa(t1, t2, p):
        return 163 * (t2 - t1)

    def int_cpb(t1, t2, p):
        return 83 * (t2 - t1)

    def int_cpc(t1, t2, p):
        return 71 * (t2 - t1)

    def r_rate(c, t, cons):
        k = 8.2e14 * np.exp(-34222 / t)  # 1/s

        return k * c["acetone"]

    a = rd.Substance(
        "acetone",
        formation_enthalpy_ig=-216.67 * 1000,
        heat_capacity_gas=cpa,
        heat_capacity_gas_dt_integral=int_cpa,
    )

    b = rd.Substance(
        "anhydride",
        formation_enthalpy_ig=-61.09 * 1000,
        heat_capacity_gas=cpb,
        heat_capacity_gas_dt_integral=int_cpb,
    )

    c = rd.Substance(
        "methane",
        formation_enthalpy_ig=-74.81 * 1000,
        heat_capacity_gas=cpc,
        heat_capacity_gas_dt_integral=int_cpc,
    )

    ref = rd.Substance("refrigerant", heat_capacity_gas=cp_ref)

    mix = rd.mix.IdealGas([a, b, c])
    mix_ref = rd.mix.IdealGas([ref])

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={"r1": {"eq": a > b + c, "rate": r_rate}},
        kinetic_constants={},
        rates_argument="concentration",
    )

    mb = pfr.mass_balances.MolarFlow(
        molar_flows_in={"acetone": 0.0376, "anhydride": 0, "methane": 0},
    )
    eb = pfr.energy_balances.NoIsothermicUConstant(
        temperature_in_or_out={"in": 1035},
        refrigerant_in_temperature=1250,
        heat_exchange_coefficient=110,
        refrigerant_mix=mix_ref,
        refrigerant_composition=np.array([1]),
        refrigerant_pressure=101325,
        refrigerant_molar_flow=0.111,
        heat_exchange_arrangement="co-current",
    )
    pb = pfr.pressure_balances.Isobaric(162 * 1000)

    area = np.pi * (4 / (16500 / 110)) ** 2 / 4

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=1 / 1000 / area,
        transversal_area=area,
        grid_size=100,
        mass_balance=mb,
        energy_balance=eb,
        pressure_balance=pb,
    )

    reactor.simulate(1e-8, bc_tol=1e-5, max_nodes=10_000)

    # Fogler data
    tr_z = (
        np.array(
            [
                0.002445,
                0.070905,
                0.163814,
                0.286064,
                0.427873,
                0.586797,
                0.792176,
                0.982885,
            ]
        )
        / 1000
        / area
    )
    tr = np.array(
        [
            1246.768,
            1187.812,
            1136.959,
            1090.980,
            1053.909,
            1025.744,
            1004.892,
            996.154,
        ]
    )

    t_z = np.array([0.029, 0.203, 0.367, 0.550, 0.782, 0.988]) / 1000 / area
    t = np.array([1031.760, 1018.966, 1009.398, 1000.654, 991.141, 984.840])

    x_z = (
        np.array(
            [
                0.0196,
                0.0685,
                0.1443,
                0.2518,
                0.3888,
                0.5330,
                0.7237,
                0.9218,
                0.9927,
            ]
        )
        / 1000
        / area
    )
    x = np.array(
        [
            0.0275,
            0.0879,
            0.1579,
            0.2307,
            0.2993,
            0.3529,
            0.4023,
            0.4325,
            0.4394,
        ]
    )

    rd_ref_temps = reactor.ode_solution.sol(tr_z)[-2]
    rd_temps = reactor.ode_solution.sol(t_z)[-3]
    rd_facetone = reactor.ode_solution.sol(x_z)[0]
    rd_x = (0.0376 - rd_facetone) / (0.0376)

    assert np.allclose(tr, rd_ref_temps, atol=5)
    assert np.allclose(t, rd_temps, atol=2)
    assert np.allclose(x, rd_x, atol=0.017)

    # temperature as out value
    eb2 = pfr.energy_balances.NoIsothermicUConstant(
        temperature_in_or_out={"out": 984.8171},
        refrigerant_in_temperature=1250,
        heat_exchange_coefficient=110,
        refrigerant_mix=mix_ref,
        refrigerant_composition=np.array([1]),
        refrigerant_pressure=101325,
        refrigerant_molar_flow=0.111,
        heat_exchange_arrangement="co-current",
    )

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=1 / 1000 / area,
        transversal_area=area,
        grid_size=1000,
        mass_balance=mb,
        energy_balance=eb2,
        pressure_balance=pb,
    )

    reactor.simulate(1e-8, bc_tol=1e-5, max_nodes=10_000)

    rd_ref_temps = reactor.ode_solution.sol(tr_z)[-2]
    rd_temps = reactor.ode_solution.sol(t_z)[-3]
    rd_facetone = reactor.ode_solution.sol(x_z)[0]
    rd_x = (0.0376 - rd_facetone) / (0.0376)

    assert np.allclose(tr, rd_ref_temps, atol=5)
    assert np.allclose(t, rd_temps, atol=4)
    assert np.allclose(x, rd_x, atol=0.017)


def test_fogler_example_12_2_case4():
    """Fogler 6th ed. example 12.2"""

    # temperature as in value
    def cpa(t, p):
        return np.full(len(t), 163)

    def cpb(t, p):
        return np.full(len(t), 83)

    def cpc(t, p):
        return np.full(len(t), 71)

    def cp_ref(t, p):
        return np.full(len(t), 34.5)

    def int_cpa(t1, t2, p):
        return 163 * (t2 - t1)

    def int_cpb(t1, t2, p):
        return 83 * (t2 - t1)

    def int_cpc(t1, t2, p):
        return 71 * (t2 - t1)

    def r_rate(c, t, cons):
        k = 8.2e14 * np.exp(-34222 / t)  # 1/s

        return k * c["acetone"]

    a = rd.Substance(
        "acetone",
        formation_enthalpy_ig=-216.67 * 1000,
        heat_capacity_gas=cpa,
        heat_capacity_gas_dt_integral=int_cpa,
    )

    b = rd.Substance(
        "anhydride",
        formation_enthalpy_ig=-61.09 * 1000,
        heat_capacity_gas=cpb,
        heat_capacity_gas_dt_integral=int_cpb,
    )

    c = rd.Substance(
        "methane",
        formation_enthalpy_ig=-74.81 * 1000,
        heat_capacity_gas=cpc,
        heat_capacity_gas_dt_integral=int_cpc,
    )

    ref = rd.Substance("refrigerant", heat_capacity_gas=cp_ref)

    mix = rd.mix.IdealGas([a, b, c])
    mix_ref = rd.mix.IdealGas([ref])

    kinetic = rd.Kinetic(
        mix=mix,
        reactions={"r1": {"eq": a > b + c, "rate": r_rate}},
        kinetic_constants={},
        rates_argument="concentration",
    )

    mb = pfr.mass_balances.MolarFlow(
        molar_flows_in={"acetone": 0.0376, "anhydride": 0, "methane": 0},
    )
    eb = pfr.energy_balances.NoIsothermicUConstant(
        temperature_in_or_out={"in": 1035},
        refrigerant_in_temperature=1250,
        heat_exchange_coefficient=110,
        refrigerant_mix=mix_ref,
        refrigerant_composition=np.array([1]),
        refrigerant_pressure=101325,
        refrigerant_molar_flow=0.111,
        heat_exchange_arrangement="counter-current",
    )
    pb = pfr.pressure_balances.Isobaric(162 * 1000)

    area = np.pi * (4 / (16500 / 110)) ** 2 / 4

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=1 / 1000 / area,
        transversal_area=area,
        grid_size=100,
        mass_balance=mb,
        energy_balance=eb,
        pressure_balance=pb,
    )

    reactor.simulate(1e-8, bc_tol=1e-5, max_nodes=10_000)

    # Fogler data
    tr_z = (
        np.array(
            [
                0.020,
                0.107,
                0.173,
                0.261,
                0.322,
                0.410,
                0.515,
                0.615,
                0.720,
                0.812,
                0.907,
                0.978,
            ]
        )
        / 1000
        / area
    )
    tr = np.array(
        [
            993.134,
            985.970,
            984.776,
            985.970,
            990.746,
            1001.493,
            1021.791,
            1049.254,
            1083.881,
            1124.478,
            1178.209,
            1231.940,
        ]
    )

    t_z = (
        np.array(
            [
                0.0146,
                0.0927,
                0.1537,
                0.2732,
                0.4073,
                0.5098,
                0.6293,
                0.7220,
                0.8366,
                0.9415,
                0.9878,
            ]
        )
        / 1000
        / area
    )
    t = np.array(
        [
            1024.1791,
            995.5224,
            984.7761,
            974.0299,
            972.8358,
            976.4179,
            984.7761,
            993.1343,
            1007.4627,
            1024.1791,
            1032.5373,
        ]
    )

    x_z = (
        np.array(
            [
                0.0073,
                0.0636,
                0.1540,
                0.2225,
                0.3105,
                0.4352,
                0.5452,
                0.6577,
                0.7531,
                0.8289,
                0.8875,
                0.9584,
                0.9878,
            ]
        )
        / 1000
        / area
    )
    x = np.array(
        [
            0.0143,
            0.0537,
            0.0955,
            0.1146,
            0.1313,
            0.1528,
            0.1731,
            0.1982,
            0.2245,
            0.2496,
            0.2782,
            0.3212,
            0.3415,
        ]
    )

    rd_ref_temps = reactor.ode_solution.sol(tr_z)[-2]
    rd_temps = reactor.ode_solution.sol(t_z)[-3]
    rd_facetone = reactor.ode_solution.sol(x_z)[0]
    rd_x = (0.0376 - rd_facetone) / (0.0376)

    assert np.allclose(tr, rd_ref_temps, atol=5)
    assert np.allclose(t, rd_temps, atol=3)
    assert np.allclose(x, rd_x, atol=0.01)

    # temperature as out value
    eb2 = pfr.energy_balances.NoIsothermicUConstant(
        temperature_in_or_out={"out": 1034.475},
        refrigerant_in_temperature=1250,
        heat_exchange_coefficient=110,
        refrigerant_mix=mix_ref,
        refrigerant_composition=np.array([1]),
        refrigerant_pressure=101325,
        refrigerant_molar_flow=0.111,
        heat_exchange_arrangement="counter-current",
    )

    reactor = pfr.PFR(
        kinetic=kinetic,
        reactor_length=1 / 1000 / area,
        transversal_area=area,
        grid_size=1000,
        mass_balance=mb,
        energy_balance=eb2,
        pressure_balance=pb,
    )

    reactor.simulate(1e-8, bc_tol=1e-10, max_nodes=10_000)

    rd_ref_temps = reactor.ode_solution.sol(tr_z)[-2]
    rd_temps = reactor.ode_solution.sol(t_z)[-3]
    rd_facetone = reactor.ode_solution.sol(x_z)[0]
    rd_x = (0.0376 - rd_facetone) / (0.0376)

    assert np.allclose(tr, rd_ref_temps, atol=5)
    assert np.allclose(t, rd_temps, atol=4)
    assert np.allclose(x, rd_x, atol=0.015)
