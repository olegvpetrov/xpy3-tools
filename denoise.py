"""
PCA-based denoising in 1D: a Toeplitz matrix approach (aka Cadzow’s denoising)

Usage requires Python 3.8+ environment

 - from shell (terminal):     python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of October 8, 2023 
 
"""
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

# get dataset displayed on the screen as a NMRDataSet object
top = Topspin()
dp = top.getDataProvider()
proton = dp.getCurrentDataset()

if not proton or proton.getDimension() != 1: 
    top.showError("Please open a 1D data set ")
    sys.exit(0)

# read fid as numpy array
a = NMRDataSetAttributes(proton)
dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])
data = dta.readfidc(rmGRPDLY=False)

TDeff = int(proton.getPar('TDeff'))
if TDeff >= 4:
    data = data[:TDeff//2]

# construct a Toeplitz matrix
row = data.size // 2                  # almost square matrix
mat = toeplitz(data[row-1::-1], data[row-1::1])

# low-rank data approximation
pb = ProgressbarThread(title='Denoising in progress')

pca = PCA(mat)
npc = pca.gavish_donoho_indicator()
#npc = pca.malinowski_indicator()

mat_approx = pca.approximate(n=npc)
data_approx = toeplitz_vector(mat_approx)

if TDeff >= 4:
    TD = int(proton.getPar('TD'))
    zf = (TD - TDeff)//2
    data_approx = np.r_[data_approx, np.zeros((zf,), dtype=complex)]

# switch to a new EXPNO
new_expno = a.get_max_expno() + 1
proton.launch('wra '+ str(new_expno))
proton.launch('re '+ str(new_expno))
proton = dp.getCurrentDataset()

# write denoised data to new EXPNO 
dta = brukerIO.dataset([a.name, str(new_expno), a.procno, a.dir])
dta.writefidc(data_approx)

# process denoised data
proton.launch('efp') 
pb.close()
