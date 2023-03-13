import numpy as np

import pytest

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
        assert mixture.mol_fractions(z) == 1.0

        for t in temperature:
            np.allclose(
                mixture.volume(z, t, pressure),
                (hexane.volume_liquid(t, pressure)),
            )

            assert mixture.concentrations(z, t, pressure) == (
                (1 / hexane.volume_liquid(t, pressure))
            )

            assert mixture.mix_heat_capacity(z, t, pressure) == (
                hexane.heat_capacity_liquid(t, pressure)
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
        raw_total_moles = sum(moles)
        raw_mol_fraction = moles / raw_total_moles

        # Test of mol_fractions method
        assert (mixture.mol_fractions(moles) == raw_mol_fraction).all()  # OKAY

        # Test of volume method
        vol_mix = mixture.volume(moles, t, pressure)
        assert vol_mix == np.dot(volumes, moles)  # OKAY

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

        # Test of mixture_molecular_weight method
        raw_molecular_weight = [
            substance.molecular_weight for substance in mixture.substances
        ]
        raw_mixture_molecular_weight = np.dot(
            raw_mol_fraction, raw_molecular_weight
        )
        assert np.allclose(
            raw_mixture_molecular_weight,
            mixture.mixture_molecular_weight(moles),
        )
        # OKAY

        # Test of molar_density method
        for t in temperature:
            raw_pure_molar_volumes = np.array(
                [
                    substance.volume_liquid(t, pressure)
                    for substance in mixture.substances
                ]
            )

            raw_molar_density = raw_total_moles / np.dot(
                raw_pure_molar_volumes, moles
            )

            assert np.allclose(
                raw_molar_density,
                mixture.molar_density(moles, t, pressure),
            )
            # OKAY

            # Test of mass_density method
            raw_mass_density = (
                raw_molar_density * raw_mixture_molecular_weight / 1000
            )
            assert np.allclose(
                raw_mass_density, mixture.mass_density(moles, t, pressure)
            )


def test_formation_enthalpies():
    substance1 = rd.Substance(formation_enthalpy=1000)
    substance2 = rd.Substance(formation_enthalpy=2000)
    substance3 = rd.Substance(formation_enthalpy=3000)

    mix = rd.mix.IdealSolution(a=substance1, b=substance2, c=substance3)

    enthalpies_mix = mix._formation_enthalpies_set()

    assert enthalpies_mix[0] == 1000

    assert enthalpies_mix[1] == 2000

    assert enthalpies_mix[2] == 3000


def test_not_implemented():
    hexane = rd.Substance.from_thermo_database("hexane")
    toluene = rd.Substance.from_thermo_database("toluene")
    butanol = rd.Substance.from_thermo_database("butanol")

    mixture = rd.mix.IdealSolution(A=hexane, B=toluene, C=butanol)

    with pytest.raises(NotImplementedError):
        mixture.mixture_viscosity(
            temperature=293, pressure=101300, moles=[1, 1, 1]
        )
