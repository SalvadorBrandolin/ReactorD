import reactord as rd
import numpy as np
import reactord.flowreactors.stationary_1d.pfr as pfr
from scipy.constants import R
import matplotlib.pyplot as plt

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

reactor.simulate(1e-5)

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
results = reactor.sim_df    
print(results)
print(results.columns)


# ============================================================
#       Results plotting
# ============================================================
fig1, (ax1, ax2, ax3) = plt.subplots (1,3)

ax1.set_xlabel ("Reactor Volume [$m^{3}$]")
ax1.set_ylabel ("Concentrations (mol/$m^{3}$)")
ax1.set_title ("Concentrations")
ax1.plot (results.z, results["but"], "-k", label="butane", linewidth=2.5)
ax1.plot (results.z, results["i-but"], "-g", label="iso-butane",  linewidth=2.5)
ax1.plot (results.z, results["i-pen"], "--c", label="iso-pentane", linewidth=1.8)
ax1.legend()
         

ax2.set_xlabel ("Reactor Volume [$m^{3}$]")
ax2.set_ylabel ("Temperature (K)")
ax2.set_title ("Temperature")
ax2.plot (results.z, results.temperature, "-r")

ax3.set_xlabel ("Reactor Volume [$m^{3}$]")
ax3.set_ylabel ("Pressure (K)")
ax3.set_title ("Pressure")
ax3.plot (results.z, results.pressure, "-b",)

plt.show()