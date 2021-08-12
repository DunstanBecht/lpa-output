<div align="center">
  <img width="250" src="https://dunstan.becht.network/views/signatures/mines.svg" alt="Mines Saint-Etienne">
</div>

# Project

This repository is related to the analysis of crystals containing dislocations by X-ray diffraction. It is part of a larger project conducted during a research internship at the laboratory of material and structural sciences of the *École Nationale Supérieure des Mines de Saint-Étienne*.

# Features

The tools developed can be used to:
* average the simulation output files
* export figures presenting the Fourier amplitudes for each harmonic
* fit different model for the calculation of the dislocation density
* export files and graphics containing information on fits

# Abbreviations

Some abbreviations are used in the figures:

* **G**: Groma
* **K**: Kamminga
* **W**: Wilkens

# User guide

### Installation

The project is indexed on [PyPI](https://pypi.org/project/lpa-output/) and installable directly via pip.
```bash
pip install -U lpa-output
```

### Analysis

To export in the directory `output-fits` the fit information on the simulation output `500_circle_1e3nm_rrdde_1e15m-2_screw` located in the directory `output-data`:

```python
from lpa.output import analyze
s = '500_circle_1e3nm_rrdde_1e15m-2_screw'
analyze.export(s, i='output-data', o='output-fits')
```
