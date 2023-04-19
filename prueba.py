import reactord as rd

a = rd.Substance.from_thermo_database(
    name="a", thermo_identification="methane"
)
b = rd.Substance.from_thermo_database(
    name="b", thermo_identification="methane"
)
c = rd.Substance.from_thermo_database(
    name="c", thermo_identification="methane"
)
d = rd.Substance.from_thermo_database(
    name="d", thermo_identification="methane"
)
e = rd.Substance.from_thermo_database(
    name="e", thermo_identification="methane"
)

mix = rd.mix.IdealGas([a, b, c, d, e])

kinetic = rd.Kinetics(
    mix,
    {
        "r1": {"eq": a + 2 * b > c, "rate": None},
        "r2": {"eq": 3 * b + 4 * c > 2 * d, "rate": None},
    },
)
print(kinetic.stoichiometry)
