import numpy as np

import reactord as rd


def test_one_substance_mix():
    hexane = rd.Substance.from_thermo_database("hexane", "hexane")
    lauric_acid = rd.Substance.from_thermo_database(
        "lauric acid", "lauric acid"
    )

    mixture = rd.mix.IdealSolution([hexane])
    mixture2 = rd.mix.IdealSolution([lauric_acid])

    compositions = np.array([[1], [10], [100], [1000], [10000]])
    temperature = np.array([300, 400, 500, 600])
    pressure = 101325

    # Lauric acid normal boiling point
    nbp_lauric_acid = lauric_acid.normal_melting_point

    for t in temperature:
        assert mixture.formation_enthalpies_correction(
            t, pressure
        ) == hexane.heat_capacity_liquid_dt_integral(298.15, t, pressure)

        cpdt_solid = lauric_acid.heat_capacity_solid_dt_integral(
            298.15, nbp_lauric_acid, pressure
        )
        dh_fus = lauric_acid.fusion_enthalpy(nbp_lauric_acid)
        cpdt_liquid = lauric_acid.heat_capacity_liquid_dt_integral(
            nbp_lauric_acid, t, pressure
        )

        assert mixture2.formation_enthalpies_correction(t, pressure) == (
            cpdt_solid + dh_fus + cpdt_liquid
        )

    for z in compositions:
        assert mixture.mole_fractions(z) == 1.0

        for t in temperature:
            assert mixture.volume(mixture.mole_fractions(z), t, pressure) == (
                hexane.volume_liquid(t, pressure),
            )

            assert mixture.concentrations(
                mixture.mole_fractions(z), t, pressure
            ) == ((1 / hexane.volume_liquid(t, pressure)))

            assert mixture.mix_heat_capacity(
                mixture.mole_fractions(z), t, pressure
            ) == (hexane.heat_capacity_liquid(t, pressure))


def test_three_substance_mix():
    hexane = rd.Substance.from_thermo_database("hexane", "hexane")
    toluene = rd.Substance.from_thermo_database("toluene", "toluene")
    butanol = rd.Substance.from_thermo_database("butanol", "butanol")

    mixture = rd.mix.IdealSolution([hexane, toluene, butanol])

    # Formation enthalpies
    fh = np.array(
        [
            hexane.formation_enthalpy,
            toluene.formation_enthalpy,
            butanol.formation_enthalpy,
        ]
    )

    assert (fh == mixture.get_formation_enthalpies()).all()

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
    pressure = 101325

    for t in temperature:
        heat_cap_correction = np.array(
            [
                hexane.heat_capacity_liquid_dt_integral(298.15, t, pressure),
                toluene.heat_capacity_liquid_dt_integral(298.15, t, pressure),
                butanol.heat_capacity_liquid_dt_integral(298.15, t, pressure),
            ]
        )

        raw_heat_capacities = np.array(
            [
                hexane.heat_capacity_liquid(t, pressure),
                toluene.heat_capacity_liquid(t, pressure),
                butanol.heat_capacity_liquid(t, pressure),
            ]
        )

        volumes = np.array(
            [
                hexane.volume_liquid(t, pressure),
                toluene.volume_liquid(t, pressure),
                butanol.volume_liquid(t, pressure),
            ]
        )

        # Test of formation_enthalpies correction method
        assert (
            mixture.formation_enthalpies_correction(t, pressure)
            == heat_cap_correction
        ).all()  # OKAY

    for moles in compositions:
        raw_mol_fraction = moles / sum(moles)

        # Test of mol_fractions method
        assert (
            mixture.mole_fractions(moles) == raw_mol_fraction
        ).all()  # OKAY

        # Test of volume method
        vol_mix = mixture.volume(raw_mol_fraction, t, pressure)
        assert vol_mix == np.dot(volumes, raw_mol_fraction)  # OKAY

        # Test of concentrations method
        total_molar_vol = np.dot(raw_mol_fraction, volumes)
        raw_concentrations = np.divide(raw_mol_fraction, total_molar_vol)
        assert np.allclose(
            mixture.concentrations(moles, t, pressure), raw_concentrations
        )

        # Test of mix_heat_capacity method
        raw_mix_heat_capacity = np.dot(raw_heat_capacities, raw_mol_fraction)
        assert np.allclose(
            mixture.mix_heat_capacity(raw_mol_fraction, t, pressure),
            raw_mix_heat_capacity,
        )  # OKAY


def test_broadcast():
    def cp1(temperature, pressure):
        cp = np.ones(np.size(temperature)) * 10
        return cp

    def cp2(temperature, pressure):
        cp = np.ones(np.size(temperature)) * 20
        return cp

    def cp3(temperature, pressure):
        vol = np.ones(np.size(temperature)) * 30
        return vol

    def vol1(temperature, pressure):
        vol = np.ones(np.size(temperature)) * 10
        return vol

    def vol2(temperature, pressure):
        vol = np.ones(np.size(temperature)) * 20
        return vol

    def vol3(temperature, pressure):
        vol = np.ones(np.size(temperature)) * 30
        return vol

    substance1 = rd.Substance(
        name="sus1",
        molecular_weight=10,
        heat_capacity_liquid=cp1,
        volume_liquid=vol1,
    )

    substance2 = rd.Substance(
        name="sus2",
        molecular_weight=20,
        heat_capacity_liquid=cp2,
        volume_liquid=vol2,
    )

    substance3 = rd.Substance(
        name="sus3",
        molecular_weight=30,
        heat_capacity_liquid=cp3,
        volume_liquid=vol3,
    )

    mixture = rd.mix.IdealSolution([substance1, substance2, substance3])

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
        volume = np.dot(mole_fractions.T[idx], [10, 20, 30])
        density = 1 / volume

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
        volume = np.dot(mole_fractions.T[idx], [10, 20, 30])

        raw_volumes[idx] = volume

    assert np.allclose(volumes, raw_volumes)

    # molar density
    molar_densities = mixture.molar_density(
        mole_fractions, temperatures, pressures
    )

    raw_molar_densities = np.ones(6)

    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        volume = np.dot(mole_fractions.T[idx], [10, 20, 30])
        molar_density = 1 / volume

        raw_molar_densities[idx] = molar_density

    assert np.allclose(molar_densities, raw_molar_densities)

    # mass density
    mass_densities = mixture.mass_density(
        mole_fractions, temperatures, pressures
    )

    raw_mass_densities = np.ones(6)

    for idx, (t, p) in enumerate(zip(temperatures, pressures)):
        volume = np.dot(mole_fractions.T[idx], [10, 20, 30])
        molar_density = 1 / volume
        mass_density = (
            np.sum(
                molar_density * mole_fractions.T[idx] * np.array([10, 20, 30])
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
