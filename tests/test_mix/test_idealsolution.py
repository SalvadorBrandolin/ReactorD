import numpy as np

import reactord as rd


def test_one_substance_mix():
    hexane = rd.Substance.from_thermo_database("hexane")
    lauric_acid = rd.Substance.from_thermo_database("lauric acid")

    mixture = rd.mix.IdealSolution([hexane])
    mixture2 = rd.mix.IdealSolution([lauric_acid])

    compositions = np.array([[1], [10], [100], [1000], [10000]])
    temperature = np.array([300, 400, 500, 600])
    pressure = 101325

    #------------ Tests Adrian Abajo -------------------------------
    #
    #--------------------------------------------------------------- 
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
    #------------ Tests Adrian Arriba-------------------------------
    #
    #--------------------------------------------------------------- 

    for z in compositions:
        assert mixture.mol_fracations(z) == 1.0

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
