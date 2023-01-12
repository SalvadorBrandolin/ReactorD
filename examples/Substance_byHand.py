# import the reactord package
import reactord as rd

# The reaction is as follows:
# A ---> B

# The substance objects are instantiated:
A = rd.Substance("water", 18, 373, 273)

# The substance object is save in a file
rd.Substance.create_substance_file(A, "name_file")

# The substance object is called
read_A = rd.Substance.load_file("name_file")

# Check
try:
    type(A) == type(read_A)
except:
    print("Pickle fail")
