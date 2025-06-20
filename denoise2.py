"""
PCA based denoising with the aid of non-uniform sampling (NUS)

For detail of the method see:
Petrov, O.V. The use of self-adaptive principal components in PCA-based denoising, JMR (2024) submitted

Usage requires Python 3.8+ environment

 - from shell (terminal):  python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of August 21, 2024 
 
"""
import argparse
import brukerIO
import numpy as np
import shutil
import statistics
import sys
import utils

from bruker.api.topspin import Topspin
from bruker.data.nmr import *

def sine_bell(dim, ssb=1, convex=False):    
    '''
    Return a window with a sine bell shape.
    
    Parameters
    ----------
    dim : int
        Number of points in the output window.
    ssb : int
        Sine-bell shift, analogous to Topspin's SSB. For ssb = 2 the window's 
        maximum is its first element. For ssb > 2 the maximum is gradually shifted
        towards the center. For ssb < 2 the maximum stays at the center. 
    convex : bool, optional
        When True, generates an upside (convex) window, for use in sampling schedules.
        
    Returns
    -------
    w : ndarray
        The window, with the maximum value normalized to 1
    '''
    assert 0 <= ssb <= dim
        
    if ssb < 2:
        w = np.sin(np.pi * (np.arange(dim) + 0.5) / dim)
    else:
        w = np.sin(np.pi / ssb + np.pi* (ssb - 1.) / ssb * (np.arange(dim) + 0.5) / dim)   
    
    if convex:    
        return (1. - w)
    else:    
        return w

def sampling_schedule(z, nus, weighted=False):
    '''
    Poisson-gap sampling with optional sinusoidal weights of the gaps
    (see Hyberts et al. JACS 132 (2010) 2145–2147; 10.1021/ja908004w)

    Parameters:  
    ----------
    z : total number of data points
    nus : fraction of sampled data points (a.k.a. 'NUS amount' in Topspin)
    weighted : when True, modulates the gaps between samples with a sine bell function
                 
    Returns:     
    ----------
    v : random numbers within [0:z-1] with Poisson distribution of gaps
    '''
    if weighted:
#        w = sine_bell(z, ssb=2)[::-1]	# same as in [Hyberts_2010]
        w = sine_bell(z, convex=True)	# complementary to [Hyberts_2010]
       
    v = np.zeros((z,), dtype=int) 
    p = nus * z 
    adj = float(z) / float(p) - 1 # initial guess of adjustment
    while True:
        i, n = 0, 0
        while i < z:    
            v[n] = i
            i += 1
            if weighted:
                k = np.random.poisson(adj * w[i-1])
            else:
                k = np.random.poisson(adj)                
            i += k
            n += 1    
        if n == int(p): break 
        if n > p: adj *= 1.02
        if n < p: adj /= 1.02 
              
    return v[:n]

# Special case of k-space sampling - with emphasis to central rows
def sampling_schedule_mri(z, nus, weighted=False):
    v1 = sampling_schedule(z//2+1, nus, weighted=weighted)
    v2 = sampling_schedule(z//2, nus, weighted=weighted) 
    return np.r_[-v1[::-1], v2[1:]]+z//2

def denoise(data, nus=None, samples=None, runs=10, weighted=True, mri=False):
    '''
    Driver for _denoise() function. Calls _denoise() twice to obtain attenuated
    versions of PCs of the original data matrix and its transpose. Thus obtained 
    PCs are transformed into 2-D PC layers for use in data reconstruction.
     
    Parameters:
    ----------
    data : dim1 x dim2 array of complex data

    nus : (optional) Fraction of dim1 to sample. When given, must be in [0, 1]

    samples : (optional) Number of inner loops where data samples are generated 
              and subjected to SVD.
            
    runs : (optional) Number of outer loops where reconstructed PCs are averaged by
           simple adding the output of _denoise()
 
    Returns: 
    ----------
    approximated (denoised) data, of the original shape dim1 x dim2
    ''' 
    if nus is None:
        nus_rng = np.linspace(.1, 1., num=24)
        pb = utils.ProgressbarThread(title='Estimating NUS amount... ', maximum=len(nus_rng))

        pca_full = utils.PCA(data)
        npc = min(10, nus_rng[0] * data.shape[0])
        
        csm = []    # cosine similarity   
        for i, nus in enumerate(nus_rng):
            d = np.zeros(npc)
            for m in range(10): 
                if mri:
                    ws = sampling_schedule_mri(len(data), nus, weighted=weighted)
                else:
                    ws = sampling_schedule(len(data), nus, weighted=weighted)    
                pca = utils.PCA(data[ws])
                           
                # calculate the cosine similarity between PC vectors
                A = np.abs(pca_full.pcs[:npc])
                B = np.abs(pca.pcs[:npc])
                d += np.diag(A.dot(B.T))
                
            if pb.stop_event.is_set():
                shutil.rmtree(f'{a.dir}/{a.name}/{str(new_expno)}/')
                sys.exit()
            pb.update(i + 1)
            csm.append(d / 10)

        pb.close()

       # linear fit of csm vs. log(nus)
        y = csm[-1] - csm
        x = np.log(nus_rng)
        x-= x[-1]

        # find where mean csm drops below 0.99, watch the bounds:
        ymean = np.mean(np.outer(x, x.dot(y) / x.dot(x)), axis=-1)
        nus = min(max(0.5, nus_rng[(ymean > 0.01).sum()]), 0.9)
        
    else:
        nus = float(nus)
        assert 0. < nus < 1.
   
    pb = utils.ProgressbarThread(title='De-noising in progress', maximum=runs*2)
            
    pcs1 = _denoise(data, nus, samples, weighted, mri)
    pb.update(1) 
     
    for i in range(1, runs):
        pcs1 += _denoise(data, nus, samples, weighted, mri)
        
        if pb.stop_event.is_set():
            shutil.rmtree(f'{a.dir}/{a.name}/{str(new_expno)}/')
            sys.exit()               
        pb.update(i+1)           

    data = data.T
    
    pcs2 = _denoise(data, nus, samples, weighted, mri)
    pb.update(runs + 1)  
    
    for i in range(1, runs):
        pcs2 += _denoise(data, nus, samples, weighted, mri)
        
        if pb.stop_event.is_set():
            shutil.rmtree(f'{a.dir}/{a.name}/{str(new_expno)}/')
            sys.exit()
        pb.update(runs + i + 2) 

    pcs1 /= runs      
    pcs2 /= runs
        
    # calculate outer products between vectors stored in pcs1 and pcs2
    npc = min(len(pcs1), len(pcs2))
    pc_layers = pcs1[:npc,:,None] @ pcs2[:npc,None,:]

    # calculate Frobenius inner products between each kernel and data
    fscores = np.tensordot(pc_layers, data.conj())
    
    # filter data
    data_approx = np.tensordot(fscores.conj(), pc_layers, axes=(0, 0))
    
    pb.close()
    return data_approx.T

def _denoise(data, nus, samples, weighted, mri):
    ''' 
    Repeat data sampling (samples of size dim1 x nus) and follow by PCA to get 
    a series of dim1 x nus matrices of psudo-random PCs. The psudo-random PCs 
    of same rank are used to approximate a full-data PC of that rank. Thus
    approximated full-data PCs ('denoised PCs') are sent back to denoise() to 
    proceed with data approximation. 
    '''
    
    # PCA on the whole of data
    pca_full = utils.PCA(data)

    # allocate 3D array to store samples' PCs
    pcs = np.zeros((samples, *pca_full.pcs.shape), dtype=complex)
    
    # populte its 1st layer with full-data PCs
    pcs[0] = pca_full.pcs        

    # populate the rest of its layers with samples' PCs 
    for i in range(1, samples):
        if mri:
            ws = sampling_schedule_mri(data.shape[0], nus, weighted=weighted) 
        else:  
            ws = sampling_schedule(data.shape[0], nus, weighted=weighted)
         
        pca = utils.PCA(data[ws])
        pcs[i,:pca.npc] = pca.pcs

    # n-rank approximation of full-data PCs
    n = np.zeros(pca.npc, dtype=int)
    for i in range(len(n)):
        pca = utils.PCA(pcs[:,i])  
              
        # estimate svd threshold
        if statistics.mode(n) == 0:
            n[i] = pca.gavish_donoho_indicator() 
        else:    
            n[i] = n[i-1]

        # approximate              
        pcs[0,i] = pca.approximate(n=n[i])[0]   

    return pcs[0,:len(ws)]

# initialize NMRDataSet object with the dataset in TopSpin window 
top = Topspin()
dp = top.getDataProvider()
hsqc = dp.getCurrentDataset()

if not hsqc or hsqc.getDimension() != 2: 
    top.showError("Please open a 2D data set ")
    sys.exit(0)

# switch to a new EXPNO
a = utils.NMRDataSetAttributes(hsqc)
new_expno = a.get_max_expno() + 1
hsqc.launch('wra '+ str(new_expno))

# read a ser file as numpy 2D array
ds = brukerIO.dataset([a.name, str(new_expno), a.procno, a.dir])
data = ds.readserc(rmGRPDLY=False)

# parse external arguments if passed         
parser = argparse.ArgumentParser()    
parser.add_argument("--nus", type=float, default=None)
parser.add_argument("--samples", type=int, default=20)

# denoise data, unweight data
args = parser.parse_args()
data_denoised = denoise(data, nus=args.nus, samples=args.samples, weighted=True)

# write denoised data to a ser file  
data_denoised = np.dstack((data_denoised.real, data_denoised.imag))
data_denoised = np.array([x.ravel() for x in data_denoised])  
ds.writeser(data_denoised) 

hsqc.launch('re '+ str(new_expno))
hsqc = dp.getCurrentDataset()

