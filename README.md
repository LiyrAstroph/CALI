# CALI
A Bayesian method for the inter-calibration of spectra in reverberation mapping

reference: http://adsabs.harvard.edu/abs/2014ApJ...786L...6L

## Compiling:  
First install the third-party package GSL and LAPACKE and then adjust Makefile with your system's library path. Type 
``make`` will generate a execuutable file ``cali``.

## Usage:

To run a code, type a commend in a Linux terminal:

```Bash
./cali
```

* The code read the data from the file "``ngc5548_year1.txt``" and do intercalibration by a Monte-Carlo Markov chain （MCMC） technique. 
The inter-clibrated continum and emission line fluxes are outputed into "``opt.txt``" and "``hb.txt``", respectively. 
The factors for inter-calibration are outputed into "``factor.txt``".

* The MCMC run with 100,000 steps and use the first 50,000 as built-in steps. 

* If you want to input a different data file, change the line 22 in main.c.

* If you want run with different MCMC steps, change the lines 11 and 12 in calibrations.c
