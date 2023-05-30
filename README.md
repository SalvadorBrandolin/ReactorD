# ReactorD

![logo](https://raw.githubusercontent.com/SalvadorBrandolin/ReactorD/main/logo.png)  
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SalvadorBrandolin/ReactorD/HEAD)
<a href="https://codeclimate.com/github/SalvadorBrandolin/ReactorD/maintainability"><img src="https://api.codeclimate.com/v1/badges/a864fbe28176d9a5d410/maintainability" /></a>
<a href="https://github.com/SalvadorBrandolin/ReactorD/actions/workflows/ci.yml">
<img src="https://github.com/SalvadorBrandolin/ReactorD/actions/workflows/ci.yml/badge.svg">
</a> 
<a href='https://reactord.readthedocs.io/en/latest/?badge=latest'>
<img src='https://readthedocs.org/projects/reactord/badge/?version=latest'
alt='Documentation Status'/></a> <a href="https://github.com/leliel12/diseno_sci_sfw">
<img src="https://camo.githubusercontent.com/69644832889fa9dfcdb974614129be2fda8e4591989fd713a983a21e7fd8d1ad/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f4469536f6674436f6d7043692d46414d41462d666664613030"></a>
<a href='https://pypi.org/project/reactord/'>
<img src='https://img.shields.io/pypi/v/reactord'>
</a>
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://tldrlegal.com/license/mit-license)
![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue)
[![codecov](https://codecov.io/gh/SalvadorBrandolin/ReactorD/branch/main/graph/badge.svg?token=6E0U0F9AYU)](https://codecov.io/gh/SalvadorBrandolin/ReactorD)

ReactorD (Reactor Design) is a `Python` package whose proposal is to simulate 
and design reactors for multiple-reaction systems. The intention is to solve 
the following reactor types in stationary or not-stationary conditions: Plug 
flow (PFR) and Stirred tank (STR) 

## Available in version 0.0.1b4
- PFR (Border and initial value problems)
    - Mass Balances
        - MolarFlow
    - Energy Balances
        - Isothermic
        - Adiabatic
        - No isothermic (constant U coefficient and refrigerant temperature)
        - No isothermic (constant U coefficient)
    - Pressure Balances
        - Isobaric
        - Ergun


## Motivation
Chemical reaction engineering has as its main objective the study and 
optimization of reactive processes, usually, with a chemical reactor as the 
protagonist equipment. To design a chemical reactor, it is necessary to 
consider several physical and chemical phenomena simultaneously, such as the 
inlet and outlet molar flow of chemical substances, mass transfer, heat 
transfer, and reaction kinetics. All these contributions to the system's 
complexity, commonly lead to coupled non-linear algebraic problems, coupled 
differential equations, or either both coupled algebraic-differential 
equations that must be solved by numeric algorithms. ReactorD provides an 
interphase to configure the necessary information for the simulation of the 
chemical reactors. Also, ReactorD implements the mathematical representations 
of mass and energy balances of specific reactors for a numerical resolution.

## Requirements
You need Python 3.8+ to run ReactorD.

## Instalation
For installing _ReactorD_ you just need to:

```sh
pip install reactord
```

## Authors
Brandolín, Salvador Eduardo
(<a href=salvadorebrandolin@mi.unc.edu.ar>salvadorebrandolin@mi.unc.edu.ar</a>)
Parodi, Adrián
(<a href=adrian.parodi@mi.unc.edu.ar>adrian.parodi@mi.unc.edu.ar</a>)
Rovezzi, Juan Pablo
(<a href=juan.rovezzi@mi.unc.edu.ar>juan.rovezzi@mi.unc.edu.ar</a>)
Santos, Maricel Del Valle
(<a href=maricel.santos@mi.unc.edu.ar>maricel.santos@mi.unc.edu.ar</a>)
Scilipoti, José Antonio
(<a href=jscilipoti@mi.unc.edu.ar>jscilipoti@mi.unc.edu.ar</a>)
