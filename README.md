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
![Output plot](https://raw.githubusercontent.com/DunstanBecht/lpa-output/5a1bf26f7dc0cef1639d769f9cf606cbc9e22171/tests/fits/10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output_analysis/output_plot.svg)

### Fits data
```
  0.9.10 # v: lpa-ouput version
5.00E+13 # d: dislocation density [m^-2]
    real # t: A(L) transformation
      KR # m: model function
      F2 # f: filter
# j Lmax[nm]    error          density[nm-2]     cut-off-radius[nm]
  1     15.0  9.8e-08  4.931549288551229e-05  6.607678158070707e+03
  1     20.0  9.3e-08  4.938883421705052e-05  6.539631112076370e+03
  1     25.0  8.4e-08  4.941491049741727e-05  6.515909541068970e+03
  1     30.0  9.2e-08  4.945746252328201e-05  6.477781851184817e+03
  1     35.0  8.7e-08  4.946979595567050e-05  6.466870006521509e+03
  1     40.0  8.4e-08  4.947948377763086e-05  6.458373859572396e+03
  1     45.0  8.1e-08  4.947640662722809e-05  6.461054460493163e+03
  1     50.0  8.5e-08  4.946038050006081e-05  6.474961311229919e+03
  1     55.0  9.8e-08  4.943548704913596e-05  6.496530818794436e+03
  1     60.0  1.3e-07  4.940070632894863e-05  6.526683598863994e+03
  1     65.0  1.7e-07  4.935561086086386e-05  6.565881451777638e+03
  1     70.0  2.2e-07  4.930026967620467e-05  6.614211441722557e+03
  1     75.0  2.9e-07  4.923532321997969e-05  6.671308949535432e+03
  1     80.0  3.7e-07  4.916156644546603e-05  6.736704404985521e+03
  2     15.0  4.4e-07  4.997989346224843e-05  1.252815309192732e+03
  2     20.0  5.2e-07  4.984628132167453e-05  1.273131101919866e+03
  2     25.0  7.0e-07  4.970749350973785e-05  1.294433392013414e+03
  2     30.0  1.0e-06  4.953824808412701e-05  1.320781191398869e+03
  3     15.0  1.3e-06  4.937734136490503e-05  5.708221010898709e+02
  4     15.0  3.9e-06  4.886626422113469e-05  3.350853311390957e+02
  5     15.0  1.2e-05  4.768064254813507e-05  2.412700062577946e+02
```

### Fits plot
![GUW](https://raw.githubusercontent.com/DunstanBecht/lpa-output/5a1bf26f7dc0cef1639d769f9cf606cbc9e22171/tests/fits/10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output_analysis/fits_plot_GUW/j1_090nm.svg)
![KR](https://raw.githubusercontent.com/DunstanBecht/lpa-output/5a1bf26f7dc0cef1639d769f9cf606cbc9e22171/tests/fits/10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output_analysis/fits_plot_KR/j1_080nm.svg)

# Abbreviations

Some abbreviations are used:

* **KR**: Krivoglaz Ryaboshapka model
* **W**: Wilkens model
* **GUW**: Groma Ungár Wilkens model
* **F0**: filter 0 (fits performed on the region of L for which A is positive)
* **F1**: filter 1 (fits performed on the region of L for which A is not oscillating)
* **F2**: filter 2 (fits performed on the region of L for which ln(A)/L^2 is linear)

# User guide

The directory `tests/` contains several examples of package module usage. To become familiar with the use of these modules you should go through these scripts in the following order:
* `test_models.py`
* `test_filters.py`
* `test_collect.py`
* `test_analyze.py`

In the sources the docstrings are carefully written and it is recommended to refer to the documentation with the `help()` python command to list the available functions, classes and parameters.

The installation from PyPI does not allow the modification of the code. To edit the package and contribute to the development use the following commands in your working directory.
```bash
pip uninstall lpa-output
git clone https://github.com/DunstanBecht/lpa-output.git
pip install -e lpa-output
cd lpa-output
git branch <name_of_your_new_branch>
```
