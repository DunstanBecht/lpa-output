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

The package is indexed on [PyPI](https://pypi.org/project/lpa-output/) and installable directly via pip:
```bash
pip install -U lpa-output
```

# Examples

### Simulation output plot
![Output plot](https://raw.githubusercontent.com/DunstanBecht/lpa-output/da72a1c908881ded27c5285b7113f2d40edc94ed/tests/fits/10_rho5e14m-2_square_2000nm_RDD_d5e-4nm-2_edge_PBCR2_S0/10_rho5e14m-2_square_2000nm_RDD_d5e-4nm-2_edge_PBCR2_S0.svg)

### Fits information
```
harmonic of g;fit L max [nm];rho [m-2];Re [nm];error
1;11.1;4.988005e+14;2.858696e+03;4.069797e-13
1;14.8;4.979192e+14;2.890819e+03;5.930914e-13
1;18.5;4.965587e+14;2.940338e+03;1.819189e-12
1;22.2;4.947938e+14;3.005100e+03;4.697750e-12
1;25.9;4.929201e+14;3.074824e+03;8.668341e-12
1;29.6;4.916342e+14;3.123224e+03;1.011788e-11
1;33.3;4.919718e+14;3.110586e+03;8.914320e-12
1;37.0;4.964844e+14;2.949933e+03;5.894697e-11
1;40.7;5.113804e+14;2.497500e+03;6.995298e-10
2;11.1;4.933120e+14;6.907357e+02;2.384938e-11
2;14.8;4.909804e+14;7.093047e+02;5.555701e-11
2;18.5;4.990565e+14;6.490353e+02;8.397701e-10
```

### Plot of fits
![Groma harmonic 1](https://raw.githubusercontent.com/DunstanBecht/lpa-output/da72a1c908881ded27c5285b7113f2d40edc94ed/tests/fits/10_rho5e14m-2_square_2000nm_RDD_d5e-4nm-2_edge_PBCR2_S0/Groma/j1_033nm.svg)
![Groma harmonic 2](https://raw.githubusercontent.com/DunstanBecht/lpa-output/da72a1c908881ded27c5285b7113f2d40edc94ed/tests/fits/10_rho5e14m-2_square_2000nm_RDD_d5e-4nm-2_edge_PBCR2_S0/Groma/j2_018nm.svg)

# Abbreviations

Some abbreviations are used in the figures:

* **G**: Groma
* **K**: Kamminga
* **W**: Wilkens

# User guide

The directory `tests/` contains several examples of package module usage. The docstrings are carefully written and it is recommended to refer to the documentation with the `help()` command.
