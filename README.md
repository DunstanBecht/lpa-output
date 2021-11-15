<div align="center">
  <img width="250" src="https://dunstan.becht.network/views/signatures/mines.svg" alt="Mines Saint-Etienne">
</div>

# Line Profile Analysis - Output

This project is related to the analysis of crystals containing dislocations by X-ray diffraction. It was developed and used for a study conducted during a research internship at the laboratory of material and structural sciences of the *École Nationale Supérieure des Mines de Saint-Étienne*. This repository contains the distribution of one of the three published python packages that have been proposed to conduct line profile analyses based on simulation results:
* [`lpa-input`](https://github.com/DunstanBecht/lpa-input) (line profile analysis input generator)
* [`lpa-xrd`](https://github.com/DunstanBecht/lpa-xrd) (line profile analysis x-ray diffraction simulation program)
* [`lpa-output`](https://github.com/DunstanBecht/lpa-output) (line profile analysis output analyzer)

The repository [`lpa-workspace`](https://github.com/DunstanBecht/lpa-workspace) contains the parameters and the scripts for the generation of the data used in the study. You can then easily replicate the results obtained or use it as inspiration to take the code in hand and conduct your own calculations. The software is placed in the public domain and you can use it as you wish. However, users are encouraged to contribute to the development and report issues.

# Features

The package `lpa.output` can be used to:
* average the simulation output files
* export figures presenting the Fourier amplitudes for each harmonic
* fit different model for the calculation of the dislocation density and the outer cut-off radius
* export files and graphics containing information on fits

# Installation

The package is indexed on [PyPI](https://pypi.org/project/lpa-output/) and installable directly via pip:
```bash
pip install -U lpa-output
```

# Examples

### Simulation output plot
![Output plot](https://raw.githubusercontent.com/DunstanBecht/lpa-output/da72a1c908881ded27c5285b7113f2d40edc94ed/tests/fits/10_rho5e14m-2_square_2000nm_RDD_d5e-4nm-2_edge_PBCR2_S0/10_rho5e14m-2_square_2000nm_RDD_d5e-4nm-2_edge_PBCR2_S0.svg)

### Fits data
```
harmonic Lmax[nm]    error          density[nm-2]     cut-off-radius[nm]            fluctuation                 R0[nm]
       1     35.4  8.9e-05  4.931659460846900e-05  7.472138403441026e+03  0.000000000000000e+00  5.883427243689241e+01
       1     47.2  2.9e-12  4.983086213130913e-05  7.254472719030799e+03  3.100220882576882e+00  6.297994403501055e+01
       1     59.0  4.6e-05  4.693751614779486e-05  1.092709702407375e+04  3.480052893181717e+00  1.138047295308443e+02
       1     70.8  4.8e-05  4.799068711216424e-05  9.401578785830581e+03  3.299205342783069e+00  9.937673264031574e+01
       1     82.6  2.1e-04  4.252550232364194e-05  2.196161571174455e+04  4.649562053713223e+00  1.593358985902378e+02
       1     94.4  4.4e-04  3.543686366828022e-05  9.999999993106799e+04  9.718114911478580e+00  1.924773313654600e+02
       1    106.2  1.1e-03  3.597209845604873e-05  9.999999971094300e+04  1.229906485683035e+01  1.825603498009725e+02
       1    118.0  2.2e-03  3.679878316573691e-05  9.999999993077158e+04  1.549017367413750e+01  1.768543671078823e+02
       1    129.8  4.0e-03  3.796148586094684e-05  9.999999989514385e+04  1.895067236025119e+01  1.746819177247687e+02
       1    141.6  6.0e-03  3.910999934813934e-05  9.999999999802289e+04  2.154243228709221e+01  1.748122478811258e+02
       3     35.4  2.2e-08  3.034175344407458e-05  9.999999999999945e+04  1.632469502814023e+01  6.699281975707740e+01
       3     47.2  1.1e-02  3.335972959068456e-05  9.999999999999991e+04  2.917611718983809e+01  6.224025404952162e+01
```

### Fits plot
![Groma harmonic 1](https://raw.githubusercontent.com/DunstanBecht/lpa-output/da72a1c908881ded27c5285b7113f2d40edc94ed/tests/fits/10_rho5e14m-2_square_2000nm_RDD_d5e-4nm-2_edge_PBCR2_S0/Groma/j1_033nm.svg)
![Groma harmonic 2](https://raw.githubusercontent.com/DunstanBecht/lpa-output/da72a1c908881ded27c5285b7113f2d40edc94ed/tests/fits/10_rho5e14m-2_square_2000nm_RDD_d5e-4nm-2_edge_PBCR2_S0/Groma/j2_018nm.svg)

# Abbreviations

Some abbreviations are used:

* **GUW**: Groma Ungár Wilkens model
* **WS**: Wilkens simplified model
* **WC**: Wilkens complete model

# User guide

The directory `tests/` contains several examples of package module usage. The docstrings are carefully written and it is recommended to refer to the documentation with the `help()` command.

The installation from PyPI does not allow the modification of the code. To edit the package and contribute to the development use the following commands in your working directory.
```bash
pip uninstall lpa-output
git clone https://github.com/DunstanBecht/lpa-output.git
pip install -e lpa-output
cd lpa-output
git branch <name_of_your_new_branch>
```
