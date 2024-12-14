"""
PCA-based denoising in 1D: a Toeplitz matrix approach (aka Cadzow’s denoising)

Usage requires Python 3.8+ environment

 - from shell (terminal):     python3 denoise [--offset n]
 - from TopSpin command line: xpy3 denoise [--offset n] 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
first version on October 8, 2023

last modification on December 2, 2024:
 -- data converted to analog type prior to denoising
 -- optional argument --offset used for threshold adjustment
 -- denoised signal decorated with 5% of original noise
 
"""
import argparse
import brukerIO
import numpy as np
import sys

from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from scipy.linalg import toeplitz
from utils import ProgressbarThread, NMRDataSetAttributes, PCA
        
def toeplitz_vector(mat):
    """
    Convert Toeplitz matrix to one-dimensional data

    """
    row, col = mat.shape
    points = row+col-1
    data = np.zeros(points, dtype=mat.dtype)
    for i in range (0, points):
        data[i] = np.mean(np.diag(mat[:,:],i-row+1))
    return data        

# define external arguments
parser = argparse.ArgumentParser()
parser.add_argument("--offset", type=int, default=0)
args = parser.parse_args()

# get dataset displayed on the screen as a NMRDataSet object
top = Topspin()
dp = top.getDataProvider()
proton = dp.getCurrentDataset()

if not proton or proton.getDimension() != 1: 
    top.showError("Please open a 1D data set ")
    sys.exit(0)

a = NMRDataSetAttributes(proton)
new_expno = a.get_max_expno() + 1

if int(proton.getPar('status DSPFVS')) != -1:
    # convert to analog data and copy to a new EXPNO
    proton.launch('convdta '+ str(new_expno))
    
else:    
    # copy to a new EXPNO
    proton.launch('wra '+ str(new_expno))
    proton.launch('re '+ str(new_expno))

proton = dp.getCurrentDataset()

# read data files using brukerIO methods
dta = brukerIO.dataset([a.name, str(new_expno), a.procno, a.dir])
data = dta.readfidc(rmGRPDLY=False)

TDeff = int(proton.getPar('TDeff'))
if TDeff >= 4:
    data = data[:TDeff//2]

# construct a Toeplitz matrix
row = data.size // 2	# almost square matrix
mat = toeplitz(data[row-1::-1], data[row-1::1])

# start denoising...
pb = ProgressbarThread(title='Denoising in progress')
pca = PCA(mat)

# approximate mat using npc principal components
npc = pca.gavish_donoho_indicator() + args.offset
mat_approx = pca.approximate(n=npc)

# convert Toeplitz matrix to 1d array
data_approx = toeplitz_vector(mat_approx)
data_approx = .95*data_approx + .05*data  # decorate with 5% of original noise

if TDeff >= 4:
    TD = int(proton.getPar('TD'))
    data_approx = np.r_[data_approx, np.zeros(((TD-TDeff)//2,), dtype=complex)]
      
# write denoised data back to bruker files
dta = brukerIO.dataset([a.name, str(new_expno), a.procno, a.dir])
dta.writefidc(data_approx)

# process denoised data
proton.launch('efp') 
pb.close()
