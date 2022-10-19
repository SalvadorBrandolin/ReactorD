from setuptools import setup

authors = (
    "Brandolín, Salvador Eduardo",
    "Parodi, Adrián",
    "Rovezzi, Juan Pablo",
    "Santos, Maricel Del Valle",
    "Scilipoti, José Antonio",
)

emails = (
    "salvadorebrandolin@mi.unc.edu.ar",
    "parodiadrian@hotmail.com",
    "juan.rovezzi@mi.unc.edu.ar",
    "maricel.santos@mi.unc.edu.ar",
    "jscilipoti@unc.edu.ar",
)

setup(
    name="reactord",
    version="0.0.1",
    description="Package for the simulation and design of chemical reactors",
    author=authors,
    author_email=emails,
    url="https://github.com/SalvadorBrandolin/reactord",
    packages=["reactord"],
)
