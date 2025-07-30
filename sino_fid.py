"""
Signal-to-noise ratio (SNR) calculator for time-domain signals (1D).
Applicable to truncated signals without defined pure noise regions.

SNR is measured as the ratio of the highest intensity in the signal - either fid or echo - 
to 2 * standard deviation of the noise, according to the Topspin's definition of SNR for 
1D spectra. 

The noise component can be readly extracted from an fid (or an echo) by means of principal 
component analysis (PCA), due to a high correlation between neighboring fid points and an 
even distribution of the noise along the fid.   

Thus, the noise std is defined as a scalar sigma in X = X0 + sigma*Z, where X0 is noiseless 
data and Z the noise matrix, the matrix X being constructed from 1D signal as a Toeplitz 
matrix. After SVD of X, a Gavish and Donoho's formula of the optimum threshold for singular 
values is applied [http://arxiv.org/abs/1305.5870] to extract sigma. 

Returns a message with the SNR value reported.
    
Usage requires Python 3.8+ environment

 - from shell (terminal):  python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of July 29, 2025 
 
"""
import brukerIO
import numpy as np
import sys

from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from scipy.linalg import toeplitz
from utils import NMRDataSetAttributes, PCA  

# get dataset displayed on the screen as a NMRDataSet object
top = Topspin()
dp = top.getDataProvider()
proton = dp.getCurrentDataset()

if not proton or proton.getDimension() != 1: 
    top.showError("Please open a 1D data set ")
    sys.exit(0)

# initialize a brukerIO's dataset object
a = NMRDataSetAttributes(proton)
dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])

# read a bruker fid file in a complex numpy array
fid = dta.readfidc(rmGRPDLY=True)

# normalize fid
idx = np.argmax(np.abs(fid))
fid /= fid[idx]

# undersample fid when too long  
skip = (len(fid) / 1024).__ceil__() 
fid = fid[::skip]

# construct a Toeplitz matrix out of fid
nrows = len(fid) // 2	# almost square matrix
X = toeplitz(fid[nrows-1::-1], fid[nrows-1::1])

# estimate sigma
pca = PCA(X)
sigma = pca.get_noise_std()

# measure and report SNR
sino = 1. / (sigma * 2)
top.showMessage(f"SINO (FID): {sino:.2f}")

