.. ReactorD documentation master file, created by
   sphinx-quickstart on Mon Nov 28 11:57:31 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. title:: ReactorD Documentation


.. image:: _static/logo.png
    :align: center
    :scale: 100 %

ReactorD
===========

|image1| |image2| |image3| |image4| |image5|

.. |image1| image:: https://api.codeclimate.com/v1/badges/3551471cd4cdf37e226f/maintainability
   :target: https://codeclimate.com/github/SalvadorBrandolin/ReactorD/maintainability
   :alt: Maintanibility
.. |image2| image:: https://github.com/SalvadorBrandolin/ReactorD/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/SalvadorBrandolin/ReactorD/actions/workflows/ci.yml
   :alt: Tests
.. |image3| image:: https://readthedocs.org/projects/reactord/badge/?version=latest
   :target: https://reactord.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. |image4| image:: https://camo.githubusercontent.com/69644832889fa9dfcdb974614129be2fda8e4591989fd713a983a21e7fd8d1ad/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f4469536f6674436f6d7043692d46414d41462d666664613030
   :target: https://camo.githubusercontent.com/69644832889fa9dfcdb974614129be2fda8e4591989fd713a983a21e7fd8d1ad/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f4469536f6674436f6d7043692d46414d41462d666664613030
   :alt: DiSoftCompCi
.. |image5| image:: https://img.shields.io/pypi/v/reactord
   :target: https://pypi.org/project/reactord/
   :alt: Pypi

----

Release
===========
.. rst-class:: release

Ver. |release|

----

Contents
--------

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

.. toctree::
   :maxdepth: 1
   :caption: Modules

   modules

**ReactorD** (Reactor Design) is a python package whose proposal is to simulate and design reactors for 
multiple-reaction systems. The intention is to solve the following reactor types in stationary or 
not-stationary conditions: Plug flow (PFR) and Stirred tank (STR) 

According to requirements, the operation settings can change as follows;

- Mass Balance:
   - Homogeneous
   - Heterogeneous
   - Continuous
   - Discontinuous

- Energy Balance:
   - Isothermic
   - Non-isothermic
   - Adiabatic

- Pressure Balance:
   - Isobaric
   - Non-isobaric

Available in version 0.0.1a1
----------------------------

- Stationary PFR Isothermic - Isobaric Operation

Motivation
----------

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
