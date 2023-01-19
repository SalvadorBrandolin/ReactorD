# import the reactord package
import reactord as rd

# The reaction is as follows:
# A ---> B

# The substance objects are instantiated:
A = rd.Substance("water", 18, 373, 273)
Ah = rd.Substance.create_substance_file(A, "name_file")

# The substance objects are called

Ar = rd.Substance.load_file("name_file")

