"""
Auto-correction of baseline in 2D spectra, applied to both axes.

(For interactive baseline correction in arrayed 1D spectra, use baseline2.) 

Exploits two external algorithms: std_distribution [1] for a baseline 
classification part and ws2d [2] for a baseline interpolation part. 

[1] D. Erb (2022) pybaselines: A Python library of algorithms for the baseline 
    correction of experimental data. https://doi.org/10.5281/zenodo.5608581
[2] Whittaker core functionality used in the 'modape' package. 
    https://github.com/WFP-VAM/vam.whittaker 
    
Usage requires Python 3.8+ environment

 - from shell (terminal):  python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of Nov 17, 2025
 
"""
import brukerIO
import numpy as np
import sys
import warnings
warnings.simplefilter("ignore", UserWarning)

from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from pybaselines.classification import std_distribution
from pybaselines.utils import optimize_window
from vam.whittaker import ws2d, ws2doptv
from utils import NMRDataSetAttributes, ProgressbarThread


def get_baseline(y, mask, l=None): 
    ''' Baseline point interpolation using the Whittaker smoother.'''
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

pb = ProgressbarThread(title='Baseplane correction in progress ') 

# row-wise identification of baseline points:
mask = np.zeros(pdata.shape, dtype=bool)
imax = np.argmax(np.linalg.norm(pdata, axis=-1))
if dta.readprocpar('PH_mod') == 2:
    hw = optimize_window(pdata[imax])
else:
    hw = max(optimize_window(pdata[imax]), optimize_window(-pdata[imax]))
for i in range(pdata.shape[0]):
    mask[i] = std_distribution(pdata[i], num_std=2., half_window=hw//4)[1]['mask']    

# row-wise interpolation and subtraction:
_, lrow = get_baseline(pdata[imax], mask[imax])
for i in range(pdata.shape[0]):
    pdata[i] -= get_baseline(pdata[i], mask[i], l=lrow)
    
# column-wise identification of baseline points:
imax = np.argmax(np.linalg.norm(pdata, axis=0))
if dta.readprocpar('PH_mod', dimension=2) == 2:
    hw = optimize_window(pdata[:,imax])
else:    
    hw = max(optimize_window(pdata[:,imax]), optimize_window(-pdata[:,imax]))
for i in range(pdata.shape[1]):
    mask[:,i] = std_distribution(pdata[:,i], num_std=2., half_window=hw//4)[1]['mask']    

# column-wise interpolation and subtraction:
_, lcol = get_baseline(pdata[:,imax], mask[:,imax])
for i in range(pdata.shape[1]):
    pdata[:,i] -= get_baseline(pdata[:,i], mask[:,i], l=lcol)

# save pdata to a 2rr file, refresh Topspin window:
dta.writespect2d(pdata)
top.getDisplay().show(hsqc, newWindow=False) 

pb.close()

