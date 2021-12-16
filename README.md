[![DOI](https://zenodo.org/badge/394321358.svg)](https://zenodo.org/badge/latestdoi/394321358)
[![License: CC0-1.0](https://img.shields.io/badge/License-CC0_1.0-lightgrey.svg)](http://creativecommons.org/publicdomain/zero/1.0/)

<div align="center"><br>
  <img width="250" src="https://dunstan.becht.network/permanent/mines.svg" alt="Mines Saint-Etienne">
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
![Output plot](https://raw.githubusercontent.com/DunstanBecht/lpa-output/6a8a310af33dfa4b833e4e82e90909c980c57c3f/tests/fits/10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output_analysis/output_plot.svg)

### Fits data
```
   0.9.5 # v: lpa-ouput version
5.00E+13 # d: dislocation density [m^-2]
    real # t: A(L) transformation
    GUW1 # m: model function
      F1 # f: filter
# j Lmax[nm]    error          density[nm-2]     cut-off-radius[nm]            fluctuation                 R0[nm]
  1     35.4  4.4e-06  4.959125385642077e-05  6.376602738042500e+03  0.000000000000000e+00  1.316554740805069e+02
  1     47.2  8.9e-13  4.883976751207240e-05  7.114562714994695e+03  8.493519339482557e-01  1.179553069207781e+02
  1     59.0  1.6e-05  4.740756763852912e-05  8.746508138807403e+03  1.726095603293172e+00  1.413902939645923e+02
  1     70.8  2.2e-05  4.650493258142388e-05  9.987113814341028e+03  2.137344464847795e+00  1.524449455165236e+02
  1     82.6  3.5e-05  4.552192493633277e-05  1.156473794499751e+04  2.488045675781454e+00  1.637735429495155e+02
  1     94.4  1.2e-04  4.267361034503640e-05  1.823206412358884e+04  3.490470038845312e+00  1.892231683102219e+02
  1    106.2  1.9e-04  3.971819044190383e-05  3.112923942998335e+04  4.706332257364845e+00  2.082445507851757e+02
  1    118.0  2.7e-04  3.656755462555112e-05  6.017193821300680e+04  6.285221597958992e+00  2.243923646692267e+02
  1    129.8  4.0e-04  3.449505768970746e-05  9.999999978443695e+04  7.758891404680341e+00  2.322279216758103e+02
  1    141.6  6.9e-04  3.465085335416217e-05  9.999999974049546e+04  8.302331153570329e+00  2.283763930916437e+02
  1    153.4  1.1e-03  3.483998341727976e-05  9.999999987845041e+04  8.903746813438747e+00  2.250429415661032e+02
  1    165.2  1.6e-03  3.502275893010202e-05  9.999999994564944e+04  9.435010747230601e+00  2.227462336095755e+02
  1    177.0  2.0e-03  3.515772025088381e-05  9.999999999085916e+04  9.794724975000429e+00  2.215270940508053e+02
  2     35.4  1.6e-09  3.110370168290108e-05  9.999999999435680e+04  1.267111199310741e+01  1.119894818816409e+02
  2     47.2  4.5e-04  3.074972857616549e-05  9.999999999161207e+04  9.754220872172349e+00  1.263012448166878e+02
  2     59.0  3.7e-04  3.072511166067763e-05  9.999999998454320e+04  9.614541435388718e+00  1.269981641006365e+02
  2     70.8  8.1e-04  3.096257003324675e-05  9.999999996616039e+04  1.068644505085711e+01  1.230060672885724e+02
  2     82.6  1.9e-03  3.128830805159731e-05  9.999999997517341e+04  1.192409626918514e+01  1.198585237211971e+02
  3     35.4  1.2e-08  2.905979115215341e-05  9.999999997762759e+04  1.005357657801308e+01  8.913340084800805e+01
  3     47.2  9.7e-04  2.937455094018407e-05  9.999999998469619e+04  1.177062124851267e+01  8.446722269635328e+01
  3     59.0  2.8e-03  2.972304588560041e-05  9.999999999858203e+04  1.324791386868198e+01  8.203696783922651e+01
  4     35.4  4.4e-08  2.846606313276551e-05  1.000000000000000e+05  1.314972719955753e+01  6.341938349515824e+01
  5     35.4  1.7e-08  3.618269150912934e-05  1.000000000000000e+05  3.133303249424009e+01  4.698476717224572e+01
```

### Fits plot
![GUW1](https://raw.githubusercontent.com/DunstanBecht/lpa-output/6a8a310af33dfa4b833e4e82e90909c980c57c3f/tests/fits/10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output_analysis/fits_plot_GUW1/j1_177nm.svg)
![GUW2](https://raw.githubusercontent.com/DunstanBecht/lpa-output/6a8a310af33dfa4b833e4e82e90909c980c57c3f/tests/fits/10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output_analysis/fits_plot_GUW2/j1_094nm.svg)

# Abbreviations

Some abbreviations are used:

* **GUW1**: Groma Ungár Wilkens model
* **GUW2**: Groma Ungár Wilkens simplified model
* **W1**: Wilkens model
* **W2**: Wilkens simplified model
* **F0**: filter 0 (fits performed on the region of L for which A is positive)
* **F1**: filter 1 (fits performed on the region of L for which A is not oscillating)
* **F2**: filter 2 (fits performed on the region of L for which ln(A)/L^2 is linear)

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
