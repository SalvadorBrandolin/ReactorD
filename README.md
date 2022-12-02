# ReactorD

![logo](https://raw.githubusercontent.com/SalvadorBrandolin/ReactorD/readthedocs/logo.png)


<a href="https://github.com/SalvadorBrandolin/ReactorD/actions/workflows/CI.yml">
<img src="https://github.com/SalvadorBrandolin/ReactorD/actions/workflows/CI.yml/badge.svg">
</a> 
<a href='https://reactord.readthedocs.io/en/latest/?badge=latest'>
<img src='https://readthedocs.org/projects/reactord/badge/?version=latest'
alt='Documentation Status'/></a> <a href="https://github.com/leliel12/diseno_sci_sfw">
<img src="https://camo.githubusercontent.com/69644832889fa9dfcdb974614129be2fda8e4591989fd713a983a21e7fd8d1ad/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f4469536f6674436f6d7043692d46414d41462d666664613030"></a>
<a href='https://pypi.org/project/reactord/'>
<img src='https://img.shields.io/pypi/v/reactord'>
</a>

PyForFluids (Python-Fortran-Fluids) is a Python package focused in the
calculation of Fluid properties based on Ecuations of State (EoS). It provides
a simple interface to work from Python but also exploits the high performance
Fortran code for the more heavy calculations.

It’s designed with modularity in mind, in a way that new thermodyinamic models
are easy to add, they even can be written either in Python or Fortran.

- Multifluid equations:
	- GERG-2008 [Paper link](https://pubs.acs.org/doi/10.1021/je300655b)

- Cubic EoS:
	- PengRobinson
	- SoaveRedlichKwong
	- Mixing Rules:
		- Quadratic (Classic Van der Waals)
		- Constant $k_{ij}$ and $l_{ij}$

## Available properties
- Reduced Temperature and Density
- Ideal Helmholtz Energy (Ao)
- Residual Helmholtz Energy (Ar)
- Compresibility Factor (Z)
- Isochoric Heat (Cv)
- Isobaric Heat (Cp)
- Speed of sound (w)
- Isothermal throttling coefficent (δ)
- Pressure derivatives:
	- Temperature
	- Density
	- Volume
- Pressure (P)
- Entropy (S)
- Gibbs Free Energy (G)
- Enthalpy (H)
- Joule-Thompson coefficent
- Isoentropic exponent
- Virial Terms:
	- B
	- C

## Motivation
While nowadays there are a lot of tools for calculation of thermodynamic
properties of fluids, most of them either are hard to mantain and don't have an
integrated testing system or are embeded to other softwares (as spredsheat
software) limiting the things that can be done to that enviroment.

PyForFluids aims to be a tool:

- With high performance, since most of it's calculations are done in Fortran
- Easy to scale due to it's modular design using the power of Python objects.
- Continuosly tested (at every `push`)to spot any problems as soon as possible.

## Instalation
For installing _ReactorD_ you just need to:

```sh
pip install reactord
```

Make sure to check the requirements first!

### Requirements



##### Windows
We recommended using the Windows Subsystem for Linux 
[WSL](https://www.windowscentral.com/install-windows-subsystem-linux-windows-10)

If WSL ain't being used, the native Windows wheels will be download instead,
so no need to worry!

## Authors
Brandolín, Salvador Eduardo 
(<a href=salvadorbrandolin@mi.unc.edu.ar>salvadorbrandolin@mi.unc.edu.ar</a>)
Parodi, Adrián
(<a href=adrian.parodi@mi.unc.edu.ar>adrian.parodi@mi.unc.edu.ar</a>)
Rovezzi, Juan Pablo
(<a href=juan.rovezzi@mi.unc.edu.ar>juan.rovezzi@mi.unc.edu.ar</a>)
Santos, Maricel Del Valle
(<a href=maricel.santos@mi.unc.edu.ar>maricel.santos@mi.unc.edu.ar</a>)
Scilipoti, José Antonio
(<a href=jscilipoti@mi.unc.edu.ar>jscilipoti@mi.unc.edu.ar</a>)














