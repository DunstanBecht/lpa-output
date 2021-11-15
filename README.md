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
![Output plot](https://raw.githubusercontent.com/DunstanBecht/lpa-output/8b1f7791531426160d4d447cc74da2154d1ee301/tests/fits/10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output_analysis/output_plot.svg)

### Fits data
```
0.9.3 # v: lpa-ouput version
 real # t: A(L) transformation
 GUW1 # m: model function
   f1 # f: filter
# j Lmax[nm]    error          density[nm-2]     cut-off-radius[nm]            fluctuation                 R0[nm]
1     35.4  2.2e-05  5.046010895384647e-05  5.746480684984484e+03  0.000000000000000e+00  3.792077963560511e+00
1     47.2  5.5e-05  5.010523698573634e-05  5.966097508211400e+03  0.000000000000000e+00  5.264353229992364e+00
1     59.0  1.2e-04  4.963456740267008e-05  6.262649973210389e+03  0.000000000000000e+00  6.120550275854728e+03
1     70.8  2.7e-04  4.894072203944183e-05  6.722977129050416e+03  0.000000000000000e+00  9.013687172041084e+03
1     82.6  5.8e-04  4.786204535937146e-05  7.516482209924678e+03  0.000000000000000e+00  1.734922135425533e+03
1     94.4  1.0e-03  4.646886651167122e-05  8.723270031082595e+03  0.000000000000000e+00  1.558791063735578e+03
1    106.2  1.9e-04  4.083425754878637e-05  2.504794805909281e+04  3.984357008160724e+00  2.189627834412867e+02
1    118.0  2.7e-04  3.749439423621627e-05  4.861986897597290e+04  5.515590348235106e+00  2.353247289280539e+02
1    129.8  4.2e-04  3.444455536098196e-05  9.999999116155526e+04  7.435702275778597e+00  2.465240744437884e+02
1    141.6  7.1e-04  3.460154681099660e-05  9.999998030869415e+04  7.972362574717990e+00  2.416401288017078e+02
1    153.4  1.1e-03  3.478677102644651e-05  9.999999672801225e+04  8.553116700082276e+00  2.374100819526699e+02
1    165.2  1.6e-03  3.496926995560776e-05  9.999999860888859e+04  9.080339238122207e+00  2.343120706436531e+02
1    177.0  2.0e-03  3.510399059132397e-05  9.999999921568978e+04  9.440535858267157e+00  2.325680833362240e+02
2     35.4  5.2e-09  3.111017377981547e-05  9.999999975233839e+04  1.283725675591809e+01  1.147009858685141e+02
2     47.2  5.3e-04  3.068994824498319e-05  1.000000000000000e+05  9.427294166626867e+00  1.331398091377537e+02
2     59.0  4.7e-04  3.061274475416887e-05  1.000000000000000e+05  8.998151057410809e+00  1.357576627736544e+02
2     70.8  7.5e-04  3.081690077326794e-05  1.000000000000000e+05  9.903535725408801e+00  1.313493778912255e+02
2     82.6  1.9e-03  3.115052610643318e-05  1.000000000000000e+05  1.116210754744593e+01  1.270184243155379e+02
3     35.4  1.9e-08  2.898105478250269e-05  1.000000000000000e+05  9.656557617399429e+00  9.437426648011720e+01
3     47.2  7.7e-04  2.923568668191490e-05  1.000000000000000e+05  1.101988289287063e+01  8.978184037808876e+01
3     59.0  2.8e-03  2.961265055912415e-05  1.000000000000000e+05  1.260574318474640e+01  8.638695624645894e+01
4     35.4  2.0e-08  2.825603238961896e-05  1.000000000000000e+05  1.207503722215676e+01  6.753324779415533e+01
5     35.4  2.6e-08  3.520019399432050e-05  1.000000000000000e+05  2.990921693333900e+01  4.796256546207081e+01
```

### Fits plot
![GUW2 harmonic 1](https://raw.githubusercontent.com/DunstanBecht/lpa-output/8b1f7791531426160d4d447cc74da2154d1ee301/tests/fits/10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output_analysis/fits_plot_GUW2/j1_094nm.svg)
![GUW2 harmonic 2](https://raw.githubusercontent.com/DunstanBecht/lpa-output/8b1f7791531426160d4d447cc74da2154d1ee301/tests/fits/10_rho5e13m-2_square_3200nm_RDD_d5e-5nm-2_screw_S0_PBC1_output_analysis/fits_plot_GUW2/j2_047nm.svg)

# Abbreviations

Some abbreviations are used:

* **GUW1**: Groma Ungár Wilkens model
* **GUW2**: Groma Ungár Wilkens simplified model
* **W1**: Wilkens model
* **W2**: Wilkens simplified model

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
