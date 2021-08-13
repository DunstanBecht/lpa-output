<div align="center">
  <img width="250" src="https://dunstan.becht.network/views/signatures/mines.svg" alt="Mines Saint-Etienne">
</div>

# Line Profile Analysis - Output

This repository is related to the analysis of crystals containing dislocations by X-ray diffraction. It is part of a project conducted during a research internship at the laboratory of material and structural sciences of the *École Nationale Supérieure des Mines de Saint-Étienne*. Three python packages have been developed to conduct line profile analyses based on simulation results:
* [`lpa.input`](https://github.com/DunstanBecht/lpa-input) (line profile analysis input generator)
* [`lpa.xrd`](https://github.com/DunstanBecht/lpa-xrd) (line profile analysis x-ray diffraction simulation program)
* [`lpa.output`](https://github.com/DunstanBecht/lpa-output) (line profile analysis output analyzer)

# Features

The package `lpa.output` can be used to:
* average the simulation output files
* export figures presenting the Fourier amplitudes for each harmonic
* fit different model for the calculation of the dislocation density
* export files and graphics containing information on fits

# Installation

The package is indexed on [PyPI](https://pypi.org/project/lpa-input/) and installable directly via pip:
```bash
pip install -U lpa-output
```

# Abbreviations

Some abbreviations are used in the figures:

* **G**: Groma
* **K**: Kamminga
* **W**: Wilkens

# User guide

The directory `tests/` contains several examples of package module usage. The docstrings are carefully written and it is recommended to refer to the documentation with the `help()` command.
