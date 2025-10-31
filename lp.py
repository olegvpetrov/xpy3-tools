'''
Forward linear prediction (LP) on complex 1D data.

It starts off with PCA on a data matrix composed of lagged sections of the FID
(a.k.a. trajectory matrix), which is chosen to be of a Toeplitz format. Then it 
proceeds to iteractive projection of the PC's on a subset of zero-padded vectors 
of the matrix with a subsequent diagonal averaging. Being iterative, it takes 
longer than convential (recursive) LP procedures, so be patient.

Usage requires Python 3.8+ environment

 - from shell (terminal):  python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 
 
written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of Oct 28, 2025

Copyright 2025 Oleg Petrov 
SPDX-License-Identifier: MIT

'''
import argparse
import brukerIO
import numpy as np
import sys
import time
import utils

from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from scipy.linalg import toeplitz


def toeplitz_vector(mat):
    """
    Convert Toeplitz matrix to 1D data
    
    Credit to Guillaume Laurent & Pierre-Aymeric Gilles (denoise_nmr.py)

    """
    row, col = mat.shape
    points = row + col - 1
    data = np.zeros(points, dtype=mat.dtype)
    
    for i in range (0, points):
        data[i] = np.mean(np.diag(mat[:,:], i-row+1))
        
    return data 

    
def lpf(fid, lpbin, lpblk=0.7):
    '''   
    A forward LP procedure.
    
    Consists of iterative projection of the FID's PCs on a lower triagonal part of
    a Toeplitz matrix made from a zero-padded FID followed by diagonal averaging
    over the whole matrix. 
    
    Parameters
    ----------
    fid : array_like
        Input data (real or complex)
        
    lpbin : int
        Number of output points for LP; has similar effect to LPBIN from TopSpin's
        processing parameters.
        
    lpblk : float : {0.5,..., 0.99}
        Number of simultaneously predicted points, expressed as a fraction of
        len(fid). This has an effect similar to TopSpin's NCOEF.
         
    Returns 
    ----------
    fid_pred : ndarray
        Concatenation of input data with predicted values, of total size lpbin
    
    '''    
    N = len(fid)   
    L = int(N * (1 - lpblk))	# number of rows of a Toeplitz matrix
    K = N - L + 1       # number of columns of a Toeplitz matrix
    
    X0 = toeplitz(fid[L-1::-1], fid[L-1:])       
    pca = utils.PCA(X0)

    npc = pca.gavish_donoho_indicator()
    VVt = pca.pcs[:npc].T.conj() @ pca.pcs[:npc]

    while True:
    
        fid_pred = np.r_[fid, np.zeros(K-1)]
        norm = 0. 
    
        for j in range(30):

            X = toeplitz(fid_pred[N-1::-1], fid_pred[N-1:])
            X[:-L] = X[:-L] @ VVt
            
            fid_pred = toeplitz_vector(X)
            fid_pred[:N] = fid
                    
           # use difference between two consecutive norms as stopping criterion
            new_norm = np.linalg.norm(fid_pred[N:])
            if abs(new_norm - norm) / new_norm < 5e-3:
                break
            norm = new_norm
    
       # check if the required number of output points has been delivered 
        if len(fid_pred) >= lpbin:
            break
            
        fid = fid_pred
        N = len(fid)
        L = N - K + 1           
    return fid_pred[:lpbin]


def main():
   
    top = Topspin()
    dp = top.getDataProvider()
    
    proton = dp.getCurrentDataset()
    if not proton or proton.getDimension() != 1: 
        top.showError("This application requires a 1D data set.")
        sys.exit(0)

   # provide necessary LP parameters         
    ME_mod = proton.getPar('ME_mod')
    if ME_mod in ('0', '3', '4'):
        top.showError("Please specify ME_mod for forward LP.")
        sys.exit(0)
        
    LPBIN = int(proton.getPar('LPBIN'))
    if LPBIN == 0:
        top.showError("Please specify a non-zero LPBIN.")
        sys.exit(0) 
   
    parser = argparse.ArgumentParser()   
    parser.add_argument("--lpblk", type=float, default=0.7)
    args = parser.parse_args()
        
   # convert data set to AMX format under new_expno
    a = utils.NMRDataSetAttributes(proton)
    new_expno = a.get_max_expno() + 1  
           
    proton.launch('convdta '+ str(new_expno))
    proton = dp.getCurrentDataset()
       
   # read fid into a complex numpy array 
    dta = brukerIO.dataset([a.name, str(new_expno), a.procno, a.dir])    
    fid = dta.readfidc()
        
   # pre-processing steps
    if ME_mod in ('5', '6'): # mirror image foward LP
        ph0 = np.angle( fid[0] )
        fid *= np.exp(-1j * ph0)
        
        data = np.r_[fid[1:][::-1].conj(), fid]
        lpbin = LPBIN//2 + len(fid)
            
    else:
        data = fid
        lpbin = LPBIN//2 
         
   # run the LP procedure...
    pb = utils.ProgressbarThread(title='Forward LP in progress...') 
    time.sleep(2) 
    fid_pred = lpf(data, lpbin, lpblk=args.lpblk) 

   # post-processing steps         
    if ME_mod in ('5', '6'):
        fid_pred = fid_pred[len(fid)-1:]
        fid_pred *= np.exp(1j * ph0)
    fid_pred[:len(fid)] = fid  
    
   # write fid_pred to disk and update corresponding sample parameters
    dta.writefidc(fid_pred)    
    proton.launch('1 TD ' + str(len(fid_pred) * 2))
    
   # display fid_pred under SPECTRUM, for easy comparison with original data
    FT_mod = proton.getPar('FT_mod')
    proton.setPar('FT_mod', 'no') 

    SI = proton.getPar('SI')
    proton.setPar('SI', str(2**int(np.log2(LPBIN//2).__ceil__())))

    proton.setPar('ME_mod', 'no')
    proton.launch('trf')

    proton.setPar('FT_mod', FT_mod) 
    proton.setPar('SI', SI)
    
    pb.close()  
    
if __name__ == "__main__":
    main()

