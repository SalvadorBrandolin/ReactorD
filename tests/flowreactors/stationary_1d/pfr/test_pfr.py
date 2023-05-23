import numpy as np

import reactord as rd
import reactord.flowreactors.stationary_1d.pfr as pfr

from scipy.constants import R


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


def test_fogler_example_11_3():
    """Fogler 6th ed. example 11.3"""

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
    rd_fbut = reactor.ode_solution.sol(vol_x)[-0]
    rd_x = (f_mol * 0.9 - rd_fbut) / (f_mol * 0.9)

    print(eb.__repr__())

    assert np.allclose(temps, rd_temps, atol=0.5)
    assert np.allclose(x, rd_x, atol=0.013)