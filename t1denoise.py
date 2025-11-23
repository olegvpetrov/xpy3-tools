"""
Auto-removal of t1 noise in 2D spectra.

Uses the noise recognition method std_distribution [1] - same as in baseline correction 
scripts (baseline, baseline2, baseplane) - to select noise regions in data columns. Then 
it seeks a column with the smallest standard deviation (std) in these noise regions and 
uses this std as a scaling factor to level down noise regions in other columns.

There are two options provided at the command line, which are '-z' (or '--zero') and 
'-s' (or '--spline'). The former means that data points within the noise regions are 
set to zero instead of being leveled down. The latter means that the points are spline-
approximated by using the Whittaker type interpolation [2], with subsequent substraction 
of the approximations from the data assigned to noise. 

It is recommended to baseline correct the spectra beforehand (e.g. with baseplane).

[1] D. Erb (2022) pybaselines: A Python library of algorithms for the baseline 
    correction of experimental data. https://doi.org/10.5281/zenodo.5608581
[2] Whittaker core functionality used in the 'modape' package. 
    https://github.com/WFP-VAM/vam.whittaker

Usage requires Python 3.8+ environment

 - from shell (terminal):  python3 <script-name [option]>
 - from TopSpin command line: xpy3 <script-name [option]>

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of Nov 20, 2025
 
"""
import argparse
import brukerIO
import numpy as np
import sys
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from pybaselines.classification import std_distribution
from pybaselines.utils import optimize_window
from vam.whittaker import ws2d, ws2doptv
from utils import NMRDataSetAttributes, ProgressbarThread

import matplotlib
import matplotlib.pyplot as plt

def get_baseline(y, mask, l=None): 
    ''' Baseline point interpolation with the Whittaker smoother.'''
    w = np.array(mask, dtype='d')  
    if l is None:  
       # run version with V-curve optimization of lambda
        lrange = array.array('d',np.linspace(2,13,111))
        z, l = ws2doptv(y, w, lrange)
        return z, np.log10(l)
    else: 
       # run basic version
        z = ws2d(y, 10**l, w)
        return z
     
# create TopSpin interface object
top = Topspin()
dp = top.getDataProvider()

# initialize NMRDataSet object with current dataset in Topspin
hsqc = dp.getCurrentDataset()
if not hsqc or hsqc.getDimension() != 2: 
    top.showError("This application requires a 2D data set.")
    sys.exit(0)
    
# read 1D complex processed bruker data: 
a = NMRDataSetAttributes(hsqc)
if not a.is_processed():
    top.showError("Processed data not found.")
    sys.exit(0)

dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])
pdata = dta.readspect2d()

pb = ProgressbarThread(title='t1denoise in progress ') 

# recognize noise regions in pdata columns:
imax = np.argmax(np.linalg.norm(pdata, axis=0))
if dta.readprocpar('PH_mod') == 2 or dta.readprocpar('PH_mod', dimension=2) == 2:
    hw = optimize_window(pdata[:,imax])
else:    
    hw = max(optimize_window(pdata[:,imax]), optimize_window(-pdata[:,imax]))
mask = np.zeros(pdata.shape, dtype=bool)
for i in range(pdata.shape[1]):
    mask[:,i] = std_distribution(pdata[:,i], num_std=3., half_window=hw//4)[1]['mask']    

# reduce recognized noise according to args provided at the command line
parser = argparse.ArgumentParser()   
parser.add_argument("-z", "--zero", action="store_true")
parser.add_argument("-s", "--spline", action="store_true")
args = parser.parse_args() 

if args.zero:    
    for i in range(pdata.shape[1]):
        pdata[mask[:,i],i] = 0.

elif args.spline:
    baseline = np.zeros(pdata.shape)
    for i in range(pdata.shape[1]):
        baseline[:,i] = get_baseline(pdata[:,i], mask[:,i], l=-1.)
        pdata[mask[:,i],i] -= baseline[mask[:,i],i]  
        
else: # no option provided -> default noise reduction ('levellng down')
    noise = pdata * mask
    noise = np.where(noise == 0., np.nan, noise)

    noise_std = np.nanstd(noise, axis=0)
    scale = noise_std / np.min(noise_std)

    for i in range(pdata.shape[1]):
        pdata[mask[:,i],i] /= scale[i]

# save pdata to a 2rr file, refresh Topspin window:
dta.writespect2d(pdata)
top.getDisplay().show(hsqc, newWindow=False) 

pb.close()

