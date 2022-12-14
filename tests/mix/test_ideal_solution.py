import numpy as np

import reactord as rd


def test_one_substance_mix():
    hexane = rd.Substance.from_thermo_database("hexane")
    lauric_acid = rd.Substance.from_thermo_database("lauric acid")

    mixture = rd.mix.IdealSolution(A=hexane)
    mixture2 = rd.mix.IdealSolution(B=lauric_acid)

    compositions = np.array([[1], [10], [100], [1000], [10000]])
    temperature = np.array([300, 400, 500, 600])
    pressure = 101325

    # Lauric acid normal boiling point
    nbp_lauric_acid = lauric_acid.normal_melting_point

    for t in temperature:

        assert mixture.formation_enthalpies_correction(
            t
        ) == hexane.heat_capacity_liquid_dt_integral(298.15, t)

        cpdt_solid = lauric_acid.heat_capacity_solid_dt_integral(
            298.15, nbp_lauric_acid
        )
        dh_fus = lauric_acid.fusion_enthalpy(nbp_lauric_acid)
        cpdt_liquid = lauric_acid.heat_capacity_liquid_dt_integral(
            nbp_lauric_acid, t
        )

        assert mixture2.formation_enthalpies_correction(t) == (
            cpdt_solid + dh_fus + cpdt_liquid
        )

    for z in compositions:
        assert mixture.mol_fractions(z) == 1.0

        for t in temperature:
            assert mixture.volume(z, t, pressure) == (
                hexane.volume_liquid(t, pressure),
            )

            assert mixture.concentrations(z, t, pressure) == (
                (1 / hexane.volume_liquid(t, pressure))
            )

            assert mixture.mix_heat_capacity(z, t, pressure) == (
                hexane.heat_capacity_liquid(t)
            )


def test_three_substance_mix():
    hexane = rd.Substance.from_thermo_database("hexane")
    toluene = rd.Substance.from_thermo_database("toluene")
    butanol = rd.Substance.from_thermo_database("butanol")

    mixture = rd.mix.IdealSolution(A=hexane, B=toluene, C=butanol)

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
                hexane.heat_capacity_liquid_dt_integral(298.15, t),
                toluene.heat_capacity_liquid_dt_integral(298.15, t),
                butanol.heat_capacity_liquid_dt_integral(298.15, t),
            ]
        )

        raw_heat_capacities = np.array(
            [
                hexane.heat_capacity_liquid(t),
                toluene.heat_capacity_liquid(t),
                butanol.heat_capacity_liquid(t),
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
            mixture.formation_enthalpies_correction(t) == heat_cap_correction
        ).all()  # OKAY

    for moles in compositions:
        raw_mol_fraction = moles / sum(moles)

        # Test of mol_fractions method
        assert (mixture.mol_fractions(moles) == raw_mol_fraction).all()  # OKAY

        # Test of volume method
        vol_mix = mixture.volume(moles, t, pressure)
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
            mixture.mix_heat_capacity(moles, t, pressure),
            raw_mix_heat_capacity,
        )  # OKAY


def test_formation_enthalpies():
    substance1 = rd.Substance(formation_enthalpy=1000)
    substance2 = rd.Substance(formation_enthalpy=2000)
    substance3 = rd.Substance(formation_enthalpy=3000)

    mix = rd.mix.IdealSolution(a=substance1, b=substance2, c=substance3)

    enthalpies_mix = mix._formation_enthalpies_set()

    assert enthalpies_mix[0] == 1000

    assert enthalpies_mix[1] == 2000

    assert enthalpies_mix[2] == 3000
