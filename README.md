# CALI
A Bayesian method for the inter-calibration of spectra in reverberation mapping

reference: http://adsabs.harvard.edu/abs/2014ApJ...786L...6L

## Compiling:  
First install the third-party package GSL and LAPACKE and then adjust Makefile with your system's library path. Type 
``make`` will generate an execuutable file ``cali``.

## Usage:

To run a code, type a commend in a Linux terminal:

```Bash
./cali param.txt
```

* The code read the configurations from the file ``param.txt`` and do intercalibration by a Monte-Carlo Markov chain （MCMC） technique. ``param.txt`` specifies file names for continuum and line data, the steps of MCMC sampling and the built-in steps that are discarded.

* The inter-clibrated continum and emission line fluxes are output into ``cont.txt`` and ``line.txt``, respectively. 
The factors for inter-calibration are output into ``factor.txt``. As a by-product, the code also reconstructs light curves, output into ``cont_recon.txt`` and ``line_recon.txt``, respectively.


## Author:
Yan-Rong Li,  liyanrong@ihep.ac.cn

Please contact me if you have any problem.
