xpy3-tools is a collection of Python scripts for NMR data processing, designed to be suitable for use with Bruker Topspin software. As such, they rely on the Bruker Python API available since Topspin 4.1.3. The scripts primarily target the author's practice, but hopefully might be of interest to other NMR practitioners. The present version (v1.5) of xpy3-tools comprises 14 scripts, which number is meant to be extended in newer versions, as time goes on and the to-do list tend to grow.

Since Topspin 4.2, the scripts can be run directly from the Topspin command line by typing in a command xpy3 script_name. In this respect, they are similar to conventional Jython scripts and Bruker AU programs. However, they can be run as usual Python programs from a console or from a dedicated IDE such as Spyder or VSCode, provided that the relevant dataset is open in a Topspin window and the embedded web service is running. For detailed description of Topspin Python interface see 
https://www.bruker.com/en/products-and-solutions/mr/nmr-software/topspin/topspin-python-interface.html. 

xpy3-tools do not modify in place raw data from fid, ser, acqu, etc. files. If the raw data (time-domain NMR signals) are programmed to change, the whole data set will be copied under a new EXPNO thus keeping the original data intact. Otherwise the data are processed in-place, files with processed data being saved under the EXPNO from whence the script has been called.

## List of scripts: 
* baseline/baseline2 : baseline correction for a single spectrum or a series of spectra
* cpmg_fp : whole-echo processing for CPMG signal readouts
* cpmg_recursive : T2 relaxometry with mono-exponential recursive model
* cpmg_relaxation : T2 relaxometry with mono- and stretched-exponential fit
* denoise/denoise2 : noise removal in 1D/2D time-domain signals
* icoshift : peak alignment over a series of spectra
* lp : forward linear prediction on 1D data
* phase_fid : zero-order phasing (rotation) of time-domain signals (1D only)
* rotsync/rotsync_plus : reconstruction of rotor-syncronized MAS spectra from arbitrarily sampled fid's
* sino_fid : signal-to-noise ratio calculator for time-domain signals (1D only)
* symmetrize : symmetrization of spinning sidebands in a MAS spectrum of half-integer nuclei
* t1t2 : basic relaxation analysis on pseudo-2D spectra 
* utils : contains common utility functions and classes

## Installation
1) You can use pip to install required packages (dependencies) in Python 3 environment. Open a terminal/command prompt, cd to the folder containing a Python 3.x interpreter and run

python3.x pip3.x install -r path/to/requirements.txt path/to/ts_remote_api*.whl path/to/bruker_nmr_api*.whl

The file requirements.txt is from this repository and the Bruker API packages (the .whl files) can be found in [TOPSPIN]/python/examples or, in case of Topspin 4.1, downloaded from bruker.com. 

2) Topspin 4.1.4+ comes with a preinstalled Python 3 environment in [TOPSPIN]/python/lib/python3.x/ which includes the Bruker Python API.

3) A special care must be taken of the module BrukerIO. The package JTutils that contains this module cannot be installed via pip. So copy brukerIO.py directly from the JTutils homepage at https://github.com/jtrebosc/JTutils/blob/master/CpyLib/ and paste it in the destination folder for xpy3-tools (see next paragraph).

4) To be able to run xpy3-tools from the Topspin command line, place the whole package in [TOPSPIN]/python/examples. Should you prefer another location for the scripts, extend a Topspin search path for Python 3 modules accordingly, as described in Section 3.2 in [this_package]/doc/tutorial_v1.4.pdf.

## Citing
If you find xpy3-tools useful in your research please cite the package as:

https://github.com/olegvpetrov/xpy3-tools

## Contact
Please contact oleg.petrov@matfyz.cuni.cz with any quesitons, feedback, or suggestions.

## License
MIT



