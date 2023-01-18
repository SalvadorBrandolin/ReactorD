import numpy as np

import reactord as rd


def test_constants():
    substance = rd.Substance(
        name="some substance",
        molecular_weight=20,
        normal_boiling_point=500,
        normal_melting_point=210,
        critical_pressure=900,
        acentric_factor=0.34,
        formation_enthalpy=-100000,
        formation_enthalpy_ig=-110000,
        formation_gibbs_ig=-50000,
        vectorize_functions=False,
    )

    assert substance.name == "some substance"
    assert substance.molecular_weight == 20
    assert substance.normal_boiling_point == 500
    assert substance.normal_melting_point == 210
    assert substance.critical_pressure == 900
    assert substance.acentric_factor == 0.34
    assert substance.formation_enthalpy == -100000
    assert substance.formation_enthalpy_ig == -110000
    assert substance.formation_gibbs_ig == -50000
    assert substance.vectorize_functions == False


def test_functions_vectorize_false():
    def vaporization_enthalpy(temperature):
        return np.exp(10 / temperature)

    def sublimation_enthalpy(temperature):
        return np.exp(20 / temperature)

    def volume_solid(temperature, pressure):
        return 6.31 * temperature / pressure

    def volume_liquid(temperature, pressure):
        return 7.31 * temperature / pressure

    def volume_gas(temperature, pressure):
        return 8.31 * temperature / pressure

    def heat_capacity_solid(temperature, pressure):
        return 6.31 * temperature + 3 * pressure

    def heat_capacity_liquid(temperature, pressure):
        return 7.31 * temperature + 3 * pressure

    def heat_capacity_gas(temperature, pressure):
        return 8.31 * temperature + 3 * pressure

    def thermal_conductivity_liquid(temperature, pressure):
        return 7.31 * temperature + 10 * pressure

    def thermal_conductivity_gas(temperature, pressure):
        return 8.31 * temperature + 10 * pressure

    def viscosity_liquid(temperature, pressure):
        return 7.31 * temperature + 20 * pressure

    def viscosity_gas(temperature, pressure):
        return 8.31 * temperature + 20 * pressure

    def heat_capacity_solid_dt_integral(temperature1, temperature2, pressure):
        integral = 6.31 / 2 * (
            temperature2**2 - temperature1**2
        ) + 3 * pressure * (temperature2 - temperature1)
        return integral

    def heat_capacity_liquid_dt_integral(temperature1, temperature2, pressure):
        integral = 7.31 / 2 * (
            temperature2**2 - temperature1**2
        ) + 3 * pressure * (temperature2 - temperature1)
        return integral

    def heat_capacity_gas_dt_integral(temperature1, temperature2, pressure):
        integral = 8.31 / 2 * (
            temperature2**2 - temperature1**2
        ) + 3 * pressure * (temperature2 - temperature1)
        return integral

    substance = rd.Substance(
        name="some substance",
        vaporization_enthalpy=vaporization_enthalpy,
        sublimation_enthalpy=sublimation_enthalpy,
        volume_solid=volume_solid,
        volume_liquid=volume_liquid,
        volume_gas=volume_gas,
        heat_capacity_solid=heat_capacity_solid,
        heat_capacity_liquid=heat_capacity_liquid,
        heat_capacity_gas=heat_capacity_gas,
        thermal_conductivity_liquid=thermal_conductivity_liquid,
        thermal_conductivity_gas=thermal_conductivity_gas,
        viscosity_liquid=viscosity_liquid,
        viscosity_gas=viscosity_gas,
        heat_capacity_solid_dt_integral=heat_capacity_solid_dt_integral,
        heat_capacity_liquid_dt_integral=heat_capacity_liquid_dt_integral,
        heat_capacity_gas_dt_integral=heat_capacity_gas_dt_integral,
        vectorize_functions=False,
    )

    assert not isinstance(substance._vaporization_enthalpy, np.vectorize)
    assert not isinstance(substance._sublimation_enthalpy, np.vectorize)
    assert not isinstance(substance._volume_solid, np.vectorize)
    assert not isinstance(substance._volume_liquid, np.vectorize)
    assert not isinstance(substance._volume_gas, np.vectorize)
    assert not isinstance(substance._heat_capacity_solid, np.vectorize)
    assert not isinstance(substance._heat_capacity_liquid, np.vectorize)
    assert not isinstance(substance._heat_capacity_gas, np.vectorize)
    assert not isinstance(substance._thermal_conductivity_liquid, np.vectorize)
    assert not isinstance(substance._thermal_conductivity_gas, np.vectorize)
    assert not isinstance(substance._viscosity_liquid, np.vectorize)
    assert not isinstance(substance._viscosity_gas, np.vectorize)
    assert not isinstance(
        substance._heat_capacity_solid_dt_integral, np.vectorize
    )
    assert not isinstance(
        substance._heat_capacity_liquid_dt_integral, np.vectorize
    )
    assert not isinstance(
        substance._heat_capacity_gas_dt_integral, np.vectorize
    )

    temperature1_scalar = 300
    temperature2_scalar = 400
    pressure_scalar = 101325

    temperature1 = np.array([300, 400, 500, 600])
    temperature2 = np.array([400, 500, 600, 700])
    pressure = np.array([101325, 201325, 301325, 401325])

    # Scalar evaluation
    assert substance.vaporization_enthalpy(
        temperature1_scalar
    ) == vaporization_enthalpy(temperature1_scalar)
    assert substance.sublimation_enthalpy(
        temperature1_scalar
    ) == sublimation_enthalpy(temperature1_scalar)
    assert substance.fusion_enthalpy(
        temperature1_scalar
    ) == sublimation_enthalpy(temperature1_scalar) - vaporization_enthalpy(
        temperature1_scalar
    )
    assert substance.volume_solid(
        temperature1_scalar, pressure_scalar
    ) == volume_solid(temperature1_scalar, pressure_scalar)
    assert substance.volume_liquid(
        temperature1_scalar, pressure_scalar
    ) == volume_liquid(temperature1_scalar, pressure_scalar)
    assert substance.volume_gas(
        temperature1_scalar, pressure_scalar
    ) == volume_gas(temperature1_scalar, pressure_scalar)
    assert substance.heat_capacity_solid(
        temperature1_scalar, pressure_scalar
    ) == heat_capacity_solid(temperature1_scalar, pressure_scalar)
    assert substance.heat_capacity_liquid(
        temperature1_scalar, pressure_scalar
    ) == heat_capacity_liquid(temperature1_scalar, pressure_scalar)
    assert substance.heat_capacity_gas(
        temperature1_scalar, pressure_scalar
    ) == heat_capacity_gas(temperature1_scalar, pressure_scalar)
    assert substance.thermal_conductivity_liquid(
        temperature1_scalar, pressure_scalar
    ) == thermal_conductivity_liquid(temperature1_scalar, pressure_scalar)
    assert substance.thermal_conductivity_gas(
        temperature1_scalar, pressure_scalar
    ) == thermal_conductivity_gas(temperature1_scalar, pressure_scalar)
    assert substance.viscosity_liquid(
        temperature1_scalar, pressure_scalar
    ) == viscosity_liquid(temperature1_scalar, pressure_scalar)
    assert substance.viscosity_gas(
        temperature1_scalar, pressure_scalar
    ) == viscosity_gas(temperature1_scalar, pressure_scalar)
    assert substance.heat_capacity_solid_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    ) == heat_capacity_solid_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    )
    assert substance.heat_capacity_liquid_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    ) == heat_capacity_liquid_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    )
    assert substance.heat_capacity_gas_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    ) == heat_capacity_gas_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    )

    # vector evaluation
    assert (
        substance.vaporization_enthalpy(temperature1)
        == vaporization_enthalpy(temperature1)
    ).all()
    assert (
        substance.sublimation_enthalpy(temperature1)
        == sublimation_enthalpy(temperature1)
    ).all()
    assert (
        substance.fusion_enthalpy(temperature1)
        == sublimation_enthalpy(temperature1)
        - vaporization_enthalpy(temperature1)
    ).all()
    assert (
        substance.volume_solid(temperature1, pressure)
        == volume_solid(temperature1, pressure)
    ).all()
    assert (
        substance.volume_liquid(temperature1, pressure)
        == volume_liquid(temperature1, pressure)
    ).all()
    assert (
        substance.volume_gas(temperature1, pressure)
        == volume_gas(temperature1, pressure)
    ).all()
    assert (
        substance.heat_capacity_solid(temperature1, pressure)
        == heat_capacity_solid(temperature1, pressure)
    ).all()
    assert (
        substance.heat_capacity_liquid(temperature1, pressure)
        == heat_capacity_liquid(temperature1, pressure)
    ).all()
    assert (
        substance.heat_capacity_gas(temperature1, pressure)
        == heat_capacity_gas(temperature1, pressure)
    ).all()
    assert (
        substance.thermal_conductivity_liquid(temperature1, pressure)
        == thermal_conductivity_liquid(temperature1, pressure)
    ).all()
    assert (
        substance.thermal_conductivity_gas(temperature1, pressure)
        == thermal_conductivity_gas(temperature1, pressure)
    ).all()
    assert (
        substance.viscosity_liquid(temperature1, pressure)
        == viscosity_liquid(temperature1, pressure)
    ).all()
    assert (
        substance.viscosity_gas(temperature1, pressure)
        == viscosity_gas(temperature1, pressure)
    ).all()
    assert (
        substance.heat_capacity_solid_dt_integral(
            temperature1, temperature2, pressure
        )
        == heat_capacity_solid_dt_integral(
            temperature1, temperature2, pressure
        )
    ).all()
    assert (
        substance.heat_capacity_liquid_dt_integral(
            temperature1, temperature2, pressure
        )
        == heat_capacity_liquid_dt_integral(
            temperature1, temperature2, pressure
        )
    ).all()
    assert (
        substance.heat_capacity_gas_dt_integral(
            temperature1, temperature2, pressure
        )
        == heat_capacity_gas_dt_integral(temperature1, temperature2, pressure)
    ).all()


def test_functions_vectorize_true():
    def vaporization_enthalpy(temperature):
        return np.exp(10 / temperature)

    def sublimation_enthalpy(temperature):
        return np.exp(20 / temperature)

    def volume_solid(temperature, pressure):
        return 6.31 * temperature / pressure

    def volume_liquid(temperature, pressure):
        return 7.31 * temperature / pressure

    def volume_gas(temperature, pressure):
        return 8.31 * temperature / pressure

    def heat_capacity_solid(temperature, pressure):
        return 6.31 * temperature + 3 * pressure

    def heat_capacity_liquid(temperature, pressure):
        return 7.31 * temperature + 3 * pressure

    def heat_capacity_gas(temperature, pressure):
        return 8.31 * temperature + 3 * pressure

    def thermal_conductivity_liquid(temperature, pressure):
        return 7.31 * temperature + 10 * pressure

    def thermal_conductivity_gas(temperature, pressure):
        return 8.31 * temperature + 10 * pressure

    def viscosity_liquid(temperature, pressure):
        return 7.31 * temperature + 20 * pressure

    def viscosity_gas(temperature, pressure):
        return 8.31 * temperature + 20 * pressure

    def heat_capacity_solid_dt_integral(temperature1, temperature2, pressure):
        integral = 6.31 / 2 * (
            temperature2**2 - temperature1**2
        ) + 3 * pressure * (temperature2 - temperature1)
        return integral

    def heat_capacity_liquid_dt_integral(temperature1, temperature2, pressure):
        integral = 7.31 / 2 * (
            temperature2**2 - temperature1**2
        ) + 3 * pressure * (temperature2 - temperature1)
        return integral

    def heat_capacity_gas_dt_integral(temperature1, temperature2, pressure):
        integral = 8.31 / 2 * (
            temperature2**2 - temperature1**2
        ) + 3 * pressure * (temperature2 - temperature1)
        return integral

    substance = rd.Substance(
        name="some substance",
        vaporization_enthalpy=vaporization_enthalpy,
        sublimation_enthalpy=sublimation_enthalpy,
        volume_solid=volume_solid,
        volume_liquid=volume_liquid,
        volume_gas=volume_gas,
        heat_capacity_solid=heat_capacity_solid,
        heat_capacity_liquid=heat_capacity_liquid,
        heat_capacity_gas=heat_capacity_gas,
        thermal_conductivity_liquid=thermal_conductivity_liquid,
        thermal_conductivity_gas=thermal_conductivity_gas,
        viscosity_liquid=viscosity_liquid,
        viscosity_gas=viscosity_gas,
        heat_capacity_solid_dt_integral=heat_capacity_solid_dt_integral,
        heat_capacity_liquid_dt_integral=heat_capacity_liquid_dt_integral,
        heat_capacity_gas_dt_integral=heat_capacity_gas_dt_integral,
        vectorize_functions=True,
    )

    assert isinstance(substance._vaporization_enthalpy, np.vectorize)
    assert isinstance(substance._sublimation_enthalpy, np.vectorize)
    assert isinstance(substance._volume_solid, np.vectorize)
    assert isinstance(substance._volume_liquid, np.vectorize)
    assert isinstance(substance._volume_gas, np.vectorize)
    assert isinstance(substance._heat_capacity_solid, np.vectorize)
    assert isinstance(substance._heat_capacity_liquid, np.vectorize)
    assert isinstance(substance._heat_capacity_gas, np.vectorize)
    assert isinstance(substance._thermal_conductivity_liquid, np.vectorize)
    assert isinstance(substance._thermal_conductivity_gas, np.vectorize)
    assert isinstance(substance._viscosity_liquid, np.vectorize)
    assert isinstance(substance._viscosity_gas, np.vectorize)
    assert isinstance(substance._heat_capacity_solid_dt_integral, np.vectorize)
    assert isinstance(
        substance._heat_capacity_liquid_dt_integral, np.vectorize
    )
    assert isinstance(substance._heat_capacity_gas_dt_integral, np.vectorize)

    temperature1_scalar = 300
    temperature2_scalar = 400
    pressure_scalar = 101325

    temperature1 = np.array([300, 400, 500, 600])
    temperature2 = np.array([400, 500, 600, 700])
    pressure = np.array([101325, 201325, 301325, 401325])

    # Scalar evaluation
    assert substance.vaporization_enthalpy(
        temperature1_scalar
    ) == vaporization_enthalpy(temperature1_scalar)
    assert substance.sublimation_enthalpy(
        temperature1_scalar
    ) == sublimation_enthalpy(temperature1_scalar)
    assert substance.fusion_enthalpy(
        temperature1_scalar
    ) == sublimation_enthalpy(temperature1_scalar) - vaporization_enthalpy(
        temperature1_scalar
    )
    assert substance.volume_solid(
        temperature1_scalar, pressure_scalar
    ) == volume_solid(temperature1_scalar, pressure_scalar)
    assert substance.volume_liquid(
        temperature1_scalar, pressure_scalar
    ) == volume_liquid(temperature1_scalar, pressure_scalar)
    assert substance.volume_gas(
        temperature1_scalar, pressure_scalar
    ) == volume_gas(temperature1_scalar, pressure_scalar)
    assert substance.heat_capacity_solid(
        temperature1_scalar, pressure_scalar
    ) == heat_capacity_solid(temperature1_scalar, pressure_scalar)
    assert substance.heat_capacity_liquid(
        temperature1_scalar, pressure_scalar
    ) == heat_capacity_liquid(temperature1_scalar, pressure_scalar)
    assert substance.heat_capacity_gas(
        temperature1_scalar, pressure_scalar
    ) == heat_capacity_gas(temperature1_scalar, pressure_scalar)
    assert substance.thermal_conductivity_liquid(
        temperature1_scalar, pressure_scalar
    ) == thermal_conductivity_liquid(temperature1_scalar, pressure_scalar)
    assert substance.thermal_conductivity_gas(
        temperature1_scalar, pressure_scalar
    ) == thermal_conductivity_gas(temperature1_scalar, pressure_scalar)
    assert substance.viscosity_liquid(
        temperature1_scalar, pressure_scalar
    ) == viscosity_liquid(temperature1_scalar, pressure_scalar)
    assert substance.viscosity_gas(
        temperature1_scalar, pressure_scalar
    ) == viscosity_gas(temperature1_scalar, pressure_scalar)
    assert substance.heat_capacity_solid_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    ) == heat_capacity_solid_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    )
    assert substance.heat_capacity_liquid_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    ) == heat_capacity_liquid_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    )
    assert substance.heat_capacity_gas_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    ) == heat_capacity_gas_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    )

    # vector evaluation
    assert (
        substance.vaporization_enthalpy(temperature1)
        == vaporization_enthalpy(temperature1)
    ).all()
    assert (
        substance.sublimation_enthalpy(temperature1)
        == sublimation_enthalpy(temperature1)
    ).all()
    assert (
        substance.fusion_enthalpy(temperature1)
        == sublimation_enthalpy(temperature1)
        - vaporization_enthalpy(temperature1)
    ).all()
    assert (
        substance.volume_solid(temperature1, pressure)
        == volume_solid(temperature1, pressure)
    ).all()
    assert (
        substance.volume_liquid(temperature1, pressure)
        == volume_liquid(temperature1, pressure)
    ).all()
    assert (
        substance.volume_gas(temperature1, pressure)
        == volume_gas(temperature1, pressure)
    ).all()
    assert (
        substance.heat_capacity_solid(temperature1, pressure)
        == heat_capacity_solid(temperature1, pressure)
    ).all()
    assert (
        substance.heat_capacity_liquid(temperature1, pressure)
        == heat_capacity_liquid(temperature1, pressure)
    ).all()
    assert (
        substance.heat_capacity_gas(temperature1, pressure)
        == heat_capacity_gas(temperature1, pressure)
    ).all()
    assert (
        substance.thermal_conductivity_liquid(temperature1, pressure)
        == thermal_conductivity_liquid(temperature1, pressure)
    ).all()
    assert (
        substance.thermal_conductivity_gas(temperature1, pressure)
        == thermal_conductivity_gas(temperature1, pressure)
    ).all()
    assert (
        substance.viscosity_liquid(temperature1, pressure)
        == viscosity_liquid(temperature1, pressure)
    ).all()
    assert (
        substance.viscosity_gas(temperature1, pressure)
        == viscosity_gas(temperature1, pressure)
    ).all()
    assert (
        substance.heat_capacity_solid_dt_integral(
            temperature1, temperature2, pressure
        )
        == heat_capacity_solid_dt_integral(
            temperature1, temperature2, pressure
        )
    ).all()
    assert (
        substance.heat_capacity_liquid_dt_integral(
            temperature1, temperature2, pressure
        )
        == heat_capacity_liquid_dt_integral(
            temperature1, temperature2, pressure
        )
    ).all()
    assert (
        substance.heat_capacity_gas_dt_integral(
            temperature1, temperature2, pressure
        )
        == heat_capacity_gas_dt_integral(temperature1, temperature2, pressure)
    ).all()


def test_create_pickle_vectorized_false():
    def vaporization_enthalpy(temperature):
        return np.exp(10 / temperature)

    def sublimation_enthalpy(temperature):
        return np.exp(20 / temperature)

    def volume_solid(temperature, pressure):
        return 6.31 * temperature / pressure

    def volume_liquid(temperature, pressure):
        return 7.31 * temperature / pressure

    def volume_gas(temperature, pressure):
        return 8.31 * temperature / pressure

    def heat_capacity_solid(temperature, pressure):
        return 6.31 * temperature + 3 * pressure

    def heat_capacity_liquid(temperature, pressure):
        return 7.31 * temperature + 3 * pressure

    def heat_capacity_gas(temperature, pressure):
        return 8.31 * temperature + 3 * pressure

    def thermal_conductivity_liquid(temperature, pressure):
        return 7.31 * temperature + 10 * pressure

    def thermal_conductivity_gas(temperature, pressure):
        return 8.31 * temperature + 10 * pressure

    def viscosity_liquid(temperature, pressure):
        return 7.31 * temperature + 20 * pressure

    def viscosity_gas(temperature, pressure):
        return 8.31 * temperature + 20 * pressure

    def heat_capacity_solid_dt_integral(temperature1, temperature2, pressure):
        integral = 6.31 / 2 * (
            temperature2**2 - temperature1**2
        ) + 3 * pressure * (temperature2 - temperature1)
        return integral

    def heat_capacity_liquid_dt_integral(temperature1, temperature2, pressure):
        integral = 7.31 / 2 * (
            temperature2**2 - temperature1**2
        ) + 3 * pressure * (temperature2 - temperature1)
        return integral

    def heat_capacity_gas_dt_integral(temperature1, temperature2, pressure):
        integral = 8.31 / 2 * (
            temperature2**2 - temperature1**2
        ) + 3 * pressure * (temperature2 - temperature1)
        return integral

    original_substance = rd.Substance(
        name="some substance",
        molecular_weight=20,
        normal_boiling_point=500,
        normal_melting_point=210,
        critical_pressure=900,
        acentric_factor=0.34,
        formation_enthalpy=-100000,
        formation_enthalpy_ig=-110000,
        formation_gibbs_ig=-50000,
        vaporization_enthalpy=vaporization_enthalpy,
        sublimation_enthalpy=sublimation_enthalpy,
        volume_solid=volume_solid,
        volume_liquid=volume_liquid,
        volume_gas=volume_gas,
        heat_capacity_solid=heat_capacity_solid,
        heat_capacity_liquid=heat_capacity_liquid,
        heat_capacity_gas=heat_capacity_gas,
        thermal_conductivity_liquid=thermal_conductivity_liquid,
        thermal_conductivity_gas=thermal_conductivity_gas,
        viscosity_liquid=viscosity_liquid,
        viscosity_gas=viscosity_gas,
        heat_capacity_solid_dt_integral=heat_capacity_solid_dt_integral,
        heat_capacity_liquid_dt_integral=heat_capacity_liquid_dt_integral,
        heat_capacity_gas_dt_integral=heat_capacity_gas_dt_integral,
        vectorize_functions=False,
    )

    rd.Substance.to_pickle(original_substance,"substance_file")

    substance = rd.Substance.from_pickle("substance_file")

    assert substance.name == "some substance"
    assert substance.molecular_weight == 20
    assert substance.normal_boiling_point == 500
    assert substance.normal_melting_point == 210
    assert substance.critical_pressure == 900
    assert substance.acentric_factor == 0.34
    assert substance.formation_enthalpy == -100000
    assert substance.formation_enthalpy_ig == -110000
    assert substance.formation_gibbs_ig == -50000
    assert substance.vectorize_functions == False

    assert not isinstance(substance._vaporization_enthalpy, np.vectorize)
    assert not isinstance(substance._sublimation_enthalpy, np.vectorize)
    assert not isinstance(substance._volume_solid, np.vectorize)
    assert not isinstance(substance._volume_liquid, np.vectorize)
    assert not isinstance(substance._volume_gas, np.vectorize)
    assert not isinstance(substance._heat_capacity_solid, np.vectorize)
    assert not isinstance(substance._heat_capacity_liquid, np.vectorize)
    assert not isinstance(substance._heat_capacity_gas, np.vectorize)
    assert not isinstance(substance._thermal_conductivity_liquid, np.vectorize)
    assert not isinstance(substance._thermal_conductivity_gas, np.vectorize)
    assert not isinstance(substance._viscosity_liquid, np.vectorize)
    assert not isinstance(substance._viscosity_gas, np.vectorize)
    assert not isinstance(
        substance._heat_capacity_solid_dt_integral, np.vectorize
    )
    assert not isinstance(
        substance._heat_capacity_liquid_dt_integral, np.vectorize
    )
    assert not isinstance(
        substance._heat_capacity_gas_dt_integral, np.vectorize
    )

    temperature1_scalar = 300
    temperature2_scalar = 400
    pressure_scalar = 101325

    temperature1 = np.array([300, 400, 500, 600])
    temperature2 = np.array([400, 500, 600, 700])
    pressure = np.array([101325, 201325, 301325, 401325])

    # Scalar evaluation
    assert substance.vaporization_enthalpy(
        temperature1_scalar
    ) == vaporization_enthalpy(temperature1_scalar)
    assert substance.sublimation_enthalpy(
        temperature1_scalar
    ) == sublimation_enthalpy(temperature1_scalar)
    assert substance.fusion_enthalpy(
        temperature1_scalar
    ) == sublimation_enthalpy(temperature1_scalar) - vaporization_enthalpy(
        temperature1_scalar
    )
    assert substance.volume_solid(
        temperature1_scalar, pressure_scalar
    ) == volume_solid(temperature1_scalar, pressure_scalar)
    assert substance.volume_liquid(
        temperature1_scalar, pressure_scalar
    ) == volume_liquid(temperature1_scalar, pressure_scalar)
    assert substance.volume_gas(
        temperature1_scalar, pressure_scalar
    ) == volume_gas(temperature1_scalar, pressure_scalar)
    assert substance.heat_capacity_solid(
        temperature1_scalar, pressure_scalar
    ) == heat_capacity_solid(temperature1_scalar, pressure_scalar)
    assert substance.heat_capacity_liquid(
        temperature1_scalar, pressure_scalar
    ) == heat_capacity_liquid(temperature1_scalar, pressure_scalar)
    assert substance.heat_capacity_gas(
        temperature1_scalar, pressure_scalar
    ) == heat_capacity_gas(temperature1_scalar, pressure_scalar)
    assert substance.thermal_conductivity_liquid(
        temperature1_scalar, pressure_scalar
    ) == thermal_conductivity_liquid(temperature1_scalar, pressure_scalar)
    assert substance.thermal_conductivity_gas(
        temperature1_scalar, pressure_scalar
    ) == thermal_conductivity_gas(temperature1_scalar, pressure_scalar)
    assert substance.viscosity_liquid(
        temperature1_scalar, pressure_scalar
    ) == viscosity_liquid(temperature1_scalar, pressure_scalar)
    assert substance.viscosity_gas(
        temperature1_scalar, pressure_scalar
    ) == viscosity_gas(temperature1_scalar, pressure_scalar)
    assert substance.heat_capacity_solid_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    ) == heat_capacity_solid_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    )
    assert substance.heat_capacity_liquid_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    ) == heat_capacity_liquid_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    )
    assert substance.heat_capacity_gas_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    ) == heat_capacity_gas_dt_integral(
        temperature1_scalar, temperature2_scalar, pressure_scalar
    )

    # vector evaluation
    assert (
        substance.vaporization_enthalpy(temperature1)
        == vaporization_enthalpy(temperature1)
    ).all()
    assert (
        substance.sublimation_enthalpy(temperature1)
        == sublimation_enthalpy(temperature1)
    ).all()
    assert (
        substance.fusion_enthalpy(temperature1)
        == sublimation_enthalpy(temperature1)
        - vaporization_enthalpy(temperature1)
    ).all()
    assert (
        substance.volume_solid(temperature1, pressure)
        == volume_solid(temperature1, pressure)
    ).all()
    assert (
        substance.volume_liquid(temperature1, pressure)
        == volume_liquid(temperature1, pressure)
    ).all()
    assert (
        substance.volume_gas(temperature1, pressure)
        == volume_gas(temperature1, pressure)
    ).all()
    assert (
        substance.heat_capacity_solid(temperature1, pressure)
        == heat_capacity_solid(temperature1, pressure)
    ).all()
    assert (
        substance.heat_capacity_liquid(temperature1, pressure)
        == heat_capacity_liquid(temperature1, pressure)
    ).all()
    assert (
        substance.heat_capacity_gas(temperature1, pressure)
        == heat_capacity_gas(temperature1, pressure)
    ).all()
    assert (
        substance.thermal_conductivity_liquid(temperature1, pressure)
        == thermal_conductivity_liquid(temperature1, pressure)
    ).all()
    assert (
        substance.thermal_conductivity_gas(temperature1, pressure)
        == thermal_conductivity_gas(temperature1, pressure)
    ).all()
    assert (
        substance.viscosity_liquid(temperature1, pressure)
        == viscosity_liquid(temperature1, pressure)
    ).all()
    assert (
        substance.viscosity_gas(temperature1, pressure)
        == viscosity_gas(temperature1, pressure)
    ).all()
    assert (
        substance.heat_capacity_solid_dt_integral(
            temperature1, temperature2, pressure
        )
        == heat_capacity_solid_dt_integral(
            temperature1, temperature2, pressure
        )
    ).all()
    assert (
        substance.heat_capacity_liquid_dt_integral(
            temperature1, temperature2, pressure
        )
        == heat_capacity_liquid_dt_integral(
            temperature1, temperature2, pressure
        )
    ).all()
    assert (
        substance.heat_capacity_gas_dt_integral(
            temperature1, temperature2, pressure
        )
        == heat_capacity_gas_dt_integral(temperature1, temperature2, pressure)
    ).all()
