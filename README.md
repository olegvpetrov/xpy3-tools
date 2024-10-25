xpy3-tools is a collection of Python scripts for NMR data processing intended to be used as command line tools in Bruker Topspin software. The scripts have primarily targeted the author's needs and the scope of his own practice, but hopefully might be of interest to other NMR practitioners. The present version (v1.0) comprises 12 scripts, which number is supposed to be extended in newer versions, as time goes on and the tasks tend to grow. 

The scripts are intended to be run by typing the command xpy3 script_name in the Topspin command line. In this respect they are similar to conventional Jython scripts and Bruker AU programs. They can also be run as regular Python programs from a console or within dedicated IDE such as Spider or VS Code, provided that the relevant dataset is loaded in Topspin and the embedded web service is running. For detailed description of Topspin Python interface see 
https://www.bruker.com/en/products-and-solutions/mr/nmr-software/topspin/topspin-python-interface.html. 

xpy3-tools do not modify in place raw data from fid, ser, acqu, etc. files. If the raw data (time-domain signals) are programmed to change, the whole data set is copied under a new EXPNO thus keeping the original data intact. Otherwise the data are processed in-place, files with processed data being saved under the EXPNO from whence a script was called.

## List of scripts: 
* baseline/baseline2 : baseline correction for a single spectrum or a series of spectra
* cpmg_fp : whole-echo processing for CPMG signal readouts
* cpmg_recursive : T2 relaxometry with mono-exponential recursive model
* cpmg_relaxation : T2 relaxometry with mono- and stretched-exponential fit
* denoise/denoise2 : noise reduction in 1D or 2D time-domain signals
* icoshift : peak alignment in stacks of 1D spectra
* phase_fid : zero-order phasing (rotation) of time-domain signals
* rotsync/rotsync_plus : reconstruction of rotor-syncronized MAS spectra from arbitrarily sampled fids
* symmetrize : symmetrization of spinning sidebands in 1D MAS spectra of half-integer nuclei
* utils : contains common utility functions and classes

## Installation
1) You can use pip to install required packages (dependencies) in Python 3 environment. Open a terminal/command prompt, cd to the folder containing a Python 3.x interpreter and run

python3.x pip3.x install -r path/to/requirements.txt path/to/ts_remote_api*.whl path/to/bruker_nmr_api*.whl

The file requirements.txt is from this repository and the Bruker API packages (the .whl files) can be found in [TOPSPIN]/python/examples or, in case of Topspin 4.1, downloaded from bruker.com. 

2) Topspin <=4.2 comes with a preinstalled Python 3 environment in [TOPSPIN]/python/lib/python3.x/ which includes the Bruker API packages.

3) A special care must be taken of the module BrukerIO, for the package JTutils that contains this module cannot be installed via pip. Copy brukerIO.py directly from the JTutils homepage at https://github.com/jtrebosc/JTutils/blob/master/CpyLib/ and paste it in the destination folder for xpy3-tools scripts (see next paragraph).

4) To run xpy3-tools from the Topspin command line, you should place the scripts in [TOPSPIN]/python/examples. Should you prefer another location for the scripts, you need to extend a Topspin serarch path for Python 3 modules accordingly, as described in tutorial.pdf.

## Citing
If you find xpy3-tools useful in your research please cite the package as:

O.V. Petrov, xpy3-tools: Topspin add-ons for NMR data processing, XXX (2024) xx, xxx-xxx.




