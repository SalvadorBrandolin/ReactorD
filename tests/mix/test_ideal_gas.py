import numpy as np

import reactord as rd


def test_one_substance_mix():
    methane = rd.Substance.from_thermo_database("methane", "methane")

    oxygen = rd.Substance.from_thermo_database("oxygen", "oxygen")

    hydrogen = rd.Substance.from_thermo_database("hydrogen", "hydrogen")

    mix1 = rd.mix.IdealGas([methane])
    mix2 = rd.mix.IdealGas([oxygen])
    mix3 = rd.mix.IdealGas([hydrogen])

    r = 8.31446261815324  # m3⋅Pa/K/mol

    compositions = np.array([[1], [10], [100], [1000], [10000]])
    temperature = np.array([300, 400, 500, 600])
    pressure = np.array([101325, 150000, 200000, 300000])

    assert mix1.formation_enthalpies_ig == methane.formation_enthalpy_ig
    assert mix2.formation_enthalpies_ig == oxygen.formation_enthalpy_ig
    assert mix3.formation_enthalpies_ig == hydrogen.formation_enthalpy_ig

    for t, p in zip(temperature, pressure):
        assert (
            mix1.formation_enthalpies_correction(t, p)
            == methane.heat_capacity_gas_dt_integral(298.15, t, p)
        ).all()
        assert (
            mix2.formation_enthalpies_correction(t, p)
            == oxygen.heat_capacity_gas_dt_integral(298.15, t, p)
        ).all()
        assert (
            mix3.formation_enthalpies_correction(t, p)
            == hydrogen.heat_capacity_gas_dt_integral(298.15, t, p)
        ).all()

    for n in compositions:
        assert mix1.mole_fractions(n) == 1.0
        assert mix2.mole_fractions(n) == 1.0
        assert mix3.mole_fractions(n) == 1.0

        for t, p in zip(temperature, pressure):
            mole_fraction = mix1.mole_fractions(n)

            np.allclose(mix1.concentrations(mole_fraction, t, p), p / (r * t))
            np.allclose(
                mix1.concentrations(mole_fraction, t, p),
                mix2.concentrations(mole_fraction, t, p),
            )
            np.allclose(
                mix2.concentrations(mole_fraction, t, p),
                mix3.concentrations(mole_fraction, t, p),
            )

            np.allclose(r * t / p, mix2.volume(mole_fraction, t, p))
            np.allclose(
                mix2.volume(mole_fraction, t, p),
                mix3.volume(mole_fraction, t, p),
            )
            np.allclose(
                mix1.volume(mole_fraction, t, p),
                mix2.volume(mole_fraction, t, p),
            )

            assert mix1.mix_heat_capacity(mole_fraction, t, p) == (
                methane.heat_capacity_gas(t, p)
            )

            assert mix2.mix_heat_capacity(mole_fraction, t, p) == (
                oxygen.heat_capacity_gas(t, p)
            )

            assert mix3.mix_heat_capacity(mole_fraction, t, p) == (
                hydrogen.heat_capacity_gas(t, p)
            )

            assert mix1.partial_pressures(mole_fraction, t, p) == p
            assert mix2.partial_pressures(mole_fraction, t, p) == p
            assert mix3.partial_pressures(mole_fraction, t, p) == p


def test_three_substances_mix():
    co2 = rd.Substance.from_thermo_database("co2", "carbon dioxide")
    ethane = rd.Substance.from_thermo_database("ethane", "ethane")
    chlorine = rd.Substance.from_thermo_database("chlorine", "chlorine")

    mixture = rd.mix.IdealGas([co2, ethane, chlorine])

    r = 8.31446261815324  # m3⋅Pa/K/mol

    compositions = np.array(
        [
            [1, 5, 8],
            [10, 15, 20],
            [100, 50, 30],
            [1000, 8000, 500],
            [10000, 500, 4000],
        ]
    )

    temperature = np.array([300, 400, 500, 600])
    pressure = np.array([101325, 150000, 200000, 300000])

    for t, p in zip(temperature, pressure):
        raw_formation_enthalpies_correction = np.array(
            [
                co2.heat_capacity_gas_dt_integral(298.15, t, p),
                ethane.heat_capacity_gas_dt_integral(298.15, t, p),
                chlorine.heat_capacity_gas_dt_integral(298.15, t, p),
            ]
        )

        raw_heat_capacities = np.array(
            [
                co2.heat_capacity_gas(t, p),
                ethane.heat_capacity_gas(t, p),
                chlorine.heat_capacity_gas(t, p),
            ]
        )

        raw_densities = p / (r * t)

        for moles in compositions:
            raw_mol_fractions = np.divide(moles, np.sum(moles, axis=0))
            raw_concentrations = np.multiply(raw_mol_fractions, raw_densities)

            # Test of concentrations method
            assert np.allclose(
                raw_concentrations,
                mixture.concentrations(raw_mol_fractions, t, p),
            )

            # Test of volume method
            raw_vol = r * t / p
            np.allclose(mixture.volume(raw_mol_fractions, t, p), raw_vol)

            # Test of mix_heat_capacity method
            raw_mix_heat_capacity = np.dot(
                raw_heat_capacities, raw_mol_fractions
            )
            np.allclose(
                raw_mix_heat_capacity,
                mixture.mix_heat_capacity(raw_mol_fractions, t, p),
            )

            # Test of formation_enthalpies_set method
            raw_enthalpies_set = np.array(
                [
                    co2.formation_enthalpy_ig,
                    ethane.formation_enthalpy_ig,
                    chlorine.formation_enthalpy_ig,
                ]
            )

            assert (
                mixture.formation_enthalpies_ig == raw_enthalpies_set
            ).all()  # OKAY

            # Test of formation_enthalpies_correction method
            assert (
                raw_formation_enthalpies_correction
                == mixture.formation_enthalpies_correction(t, p)
            ).all()  # OKAY


def test_broadcast():
    def cp1(temperature, pressure):
        cp = np.ones(np.size(temperature)) * 10
        return cp

    def cp2(temperature, pressure):
        cp = np.ones(np.size(temperature)) * 20
        return cp

    def cp3(temperature, pressure):
        cp = np.ones(np.size(temperature)) * 30
        return cp

    substance1 = rd.Substance(
        name="sus1",
        molecular_weight=10,
        heat_capacity_gas=cp1,
    )

    substance2 = rd.Substance(
        name="sus2",
        molecular_weight=20,
        heat_capacity_gas=cp2,
    )

    substance3 = rd.Substance(
        name="sus3",
        molecular_weight=30,
        heat_capacity_gas=cp3,
    )

    mixture = rd.mix.IdealGas([substance1, substance2, substance3])

    moles = np.array(
        [[1, 2, 3, 4, 5, 6], [6, 5, 4, 3, 2, 1], [0, 2, 5, 2, 1, 5]]
    )

    # mole fractions
    mole_fractions = mixture.mole_fractions(moles)

    raw_mole_fractions = np.array(
        [
            [
                1 / (1 + 6),
                2 / (2 + 5 + 2),
                3 / (3 + 4 + 5),
                4 / (4 + 3 + 2),
                5 / (5 + 2 + 1),
                6 / (6 + 1 + 5),
            ],
            [
                6 / (1 + 6),
                5 / (2 + 5 + 2),
                4 / (3 + 4 + 5),
                3 / (4 + 3 + 2),
                2 / (5 + 2 + 1),
                1 / (6 + 1 + 5),
            ],
            [
                0 / (1 + 6),
                2 / (2 + 5 + 2),
                5 / (3 + 4 + 5),
                2 / (4 + 3 + 2),
                1 / (5 + 2 + 1),
                5 / (6 + 1 + 5),
            ],
        ]
    )

    assert np.allclose(mole_fractions, raw_mole_fractions)

    # mixture molecular weight
    mix_molecular_weight = mixture.mix_molecular_weight(mole_fractions)

    pure_molecular_weights = np.array([10, 20, 30])

    raw_mix_molecular_weight = np.array([])
    for composition in mole_fractions.T:
        raw_mix_molecular_weight = np.append(
            raw_mix_molecular_weight,
            np.dot(composition, pure_molecular_weights),
        )

    assert np.allclose(mix_molecular_weight, raw_mix_molecular_weight)

    # concentrations
    temperatures = np.array([300, 400, 500, 600, 700, 800])
    pressures = np.array([100000, 200000, 300000, 400000, 500000, 600000])

    concentrations = mixture.concentrations(
        mole_fractions, temperatures, pressures
    )

    raw_concentrations = np.ones((6, 3))

    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        r = 8.31446261815324  # m3⋅Pa/K/mol
        density = p / (r * t)

        conc = np.multiply(mole_fractions.T[idx], density)
        raw_concentrations[idx, :] = conc

    raw_concentrations = raw_concentrations.T

    assert np.allclose(concentrations, raw_concentrations)

    # partial pressures
    partial_pressures = mixture.partial_pressures(
        mole_fractions, temperatures, pressures
    )

    raw_partial_pressures = np.ones((6, 3))
    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        partial_pres = np.multiply(mole_fractions.T[idx], p)

        raw_partial_pressures[idx, :] = partial_pres

    assert np.allclose(partial_pressures, raw_partial_pressures.T)

    # volume
    volumes = mixture.volume(mole_fractions, temperatures, pressures)

    raw_volumes = np.ones(6)

    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        r = 8.31446261815324  # m3⋅Pa/K/mol
        volume = 1 / (p / (r * t))

        raw_volumes[idx] = volume

    assert np.allclose(volumes, raw_volumes)

    # molar density
    molar_densities = mixture.molar_density(
        mole_fractions, temperatures, pressures
    )

    raw_molar_densities = np.ones(6)

    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        r = 8.31446261815324  # m3⋅Pa/K/mol
        molar_density = p / (r * t)

        raw_molar_densities[idx] = molar_density

    assert np.allclose(molar_densities, raw_molar_densities)

    # mass density
    mass_densities = mixture.mass_density(
        mole_fractions, temperatures, pressures
    )

    raw_mass_densities = np.ones(6)

    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        r = 8.31446261815324  # m3⋅Pa/K/mol
        mass_density = (
            np.sum(
                p / (r * t) * mole_fractions.T[idx] * np.array([10, 20, 30])
            )
            / 1000
        )

        raw_mass_densities[idx] = mass_density

    assert np.allclose(mass_densities, raw_mass_densities)

    # mix heat capacity
    heat_capacities = mixture.mix_heat_capacity(
        mole_fractions, temperatures, pressures
    )

    raw_heat_capacities = np.ones(6)

    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        heat_capacity = np.sum(mole_fractions.T[idx] * np.array([10, 20, 30]))

        raw_heat_capacities[idx] = heat_capacity

    assert np.allclose(heat_capacities, raw_heat_capacities)


def test_broadcast_thermo():
    substance1 = rd.Substance.from_thermo_database(
        name="sus1", thermo_identification="acetic acid"
    )

    substance2 = rd.Substance.from_thermo_database(
        name="sus2", thermo_identification="water"
    )

    substance3 = rd.Substance.from_thermo_database(
        name="sus3", thermo_identification="ethanol"
    )

    mixture = rd.mix.IdealGas([substance1, substance2, substance3])

    moles = np.array(
        [[1, 2, 3, 4, 5, 6], [6, 5, 4, 3, 2, 1], [0, 2, 5, 2, 1, 5]]
    )

    # mole fractions
    mole_fractions = mixture.mole_fractions(moles)

    raw_mole_fractions = np.array(
        [
            [
                1 / (1 + 6),
                2 / (2 + 5 + 2),
                3 / (3 + 4 + 5),
                4 / (4 + 3 + 2),
                5 / (5 + 2 + 1),
                6 / (6 + 1 + 5),
            ],
            [
                6 / (1 + 6),
                5 / (2 + 5 + 2),
                4 / (3 + 4 + 5),
                3 / (4 + 3 + 2),
                2 / (5 + 2 + 1),
                1 / (6 + 1 + 5),
            ],
            [
                0 / (1 + 6),
                2 / (2 + 5 + 2),
                5 / (3 + 4 + 5),
                2 / (4 + 3 + 2),
                1 / (5 + 2 + 1),
                5 / (6 + 1 + 5),
            ],
        ]
    )

    assert np.allclose(mole_fractions, raw_mole_fractions)

    # mixture molecular weight
    mix_molecular_weight = mixture.mix_molecular_weight(mole_fractions)

    pure_molecular_weights = np.array(
        [substance.molecular_weight for substance in mixture.substances]
    )

    raw_mix_molecular_weight = np.array([])
    for composition in mole_fractions.T:
        raw_mix_molecular_weight = np.append(
            raw_mix_molecular_weight,
            np.dot(composition, pure_molecular_weights),
        )

    assert np.allclose(mix_molecular_weight, raw_mix_molecular_weight)

    # concentrations
    temperatures = np.array([300, 400, 500, 600, 700, 800])
    pressures = np.array([100000, 200000, 300000, 400000, 500000, 600000])

    concentrations = mixture.concentrations(
        mole_fractions, temperatures, pressures
    )

    raw_concentrations = np.ones((6, 3))

    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        r = 8.31446261815324  # m3⋅Pa/K/mol
        density = p / (r * t)

        conc = np.multiply(mole_fractions.T[idx], density)
        raw_concentrations[idx, :] = conc

    raw_concentrations = raw_concentrations.T

    assert np.allclose(concentrations, raw_concentrations)

    # partial pressures
    partial_pressures = mixture.partial_pressures(
        mole_fractions, temperatures, pressures
    )

    raw_partial_pressures = np.ones((6, 3))
    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        partial_pres = np.multiply(mole_fractions.T[idx], p)

        raw_partial_pressures[idx, :] = partial_pres

    assert np.allclose(partial_pressures, raw_partial_pressures.T)

    # volume
    volumes = mixture.volume(mole_fractions, temperatures, pressures)

    raw_volumes = np.ones(6)

    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        r = 8.31446261815324  # m3⋅Pa/K/mol
        volume = 1 / (p / (r * t))

        raw_volumes[idx] = volume

    assert np.allclose(volumes, raw_volumes)

    # molar density
    molar_densities = mixture.molar_density(
        mole_fractions, temperatures, pressures
    )

    raw_molar_densities = np.ones(6)

    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        r = 8.31446261815324  # m3⋅Pa/K/mol
        molar_density = p / (r * t)

        raw_molar_densities[idx] = molar_density

    assert np.allclose(molar_densities, raw_molar_densities)

    # mass density
    mass_densities = mixture.mass_density(
        mole_fractions, temperatures, pressures
    )

    raw_mass_densities = np.ones(6)

    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        r = 8.31446261815324  # m3⋅Pa/K/mol
        mass_density = (
            np.sum(
                p / (r * t) * mole_fractions.T[idx] * pure_molecular_weights
            )
            / 1000
        )

        raw_mass_densities[idx] = mass_density

    assert np.allclose(mass_densities, raw_mass_densities)

    # mix heat capacity
    heat_capacities = mixture.mix_heat_capacity(
        mole_fractions, temperatures, pressures
    )

    raw_heat_capacities = np.ones(6)

    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        pure_heat_capacities = np.array(
            [
                substance.heat_capacity_gas(t, p)
                for substance in mixture.substances
            ]
        )
        heat_capacity = np.sum(mole_fractions.T[idx] * pure_heat_capacities)

        raw_heat_capacities[idx] = heat_capacity

    assert np.allclose(heat_capacities, raw_heat_capacities)
