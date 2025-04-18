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

def sampling_schedule(z, nus, weighted=True):
    '''
    Poisson-gap sampling with optional sinusoidal weights of the gaps
    (see Hyberts et al. JACS 132 (2010) 2145–2147; 10.1021/ja908004w)

    Parameters:  z : total number of points
                 nus : fraction of sampled points (a.k.a. nus amount)
                 weighted : sinusoidal weighting of the gaps between samples
                 
    Return:      v : random numbers within [0:z-1] with Poisson-distributed gaps
    '''
    v = np.zeros((z,), dtype=int) 
    p = nus * z 
    adj = float(z) / float(p) - 1 # initial guess of adjustment
    while True:
        i, n = 0, 0
        while i < z:    
            v[n] = i
            i += 1
            if weighted:
                k = np.random.poisson(adj * (1.- np.sin((i+0.5)/(z+1)*np.pi)))
            else:
                k = np.random.poisson(adj)                
            i += k
            n += 1    
        if n == int(p): break 
        if n > p: adj *= 1.02
        if n < p: adj /= 1.02
    return v[:n]

# A special case of k-space sampling, to accent central points
def sampling_schedule_mri(z, nus, weighted=True):
    v1 = sampling_schedule(z//2+1, nus, weighted=weighted)
    v2 = sampling_schedule(z//2, nus, weighted=weighted) 
    return np.r_[-v1[::-1], v2[1:]]+z//2

def denoise(data, nus=None, samples=None, runs=10, weighted=True, mri=False):
    '''
    A driver for _denoise() function. Calls _denoise() twice to obtain attenuated
    versions of PCs of the data matrix as well as its transpose. The PCs thus
    obtained are transformed into 2D 'PC layers' (kernels) to reconstruct the data.
     
    Input:
    data : dim1 x dim2 array of complex data

    nus : (optional) The fraction of dim1 to sample, in [0, 1]

    samples : (optional) The number of inner loops for data sampling and PCA
            
    runs : (optional) The number of outer loops for averaging the reconstructed PC's
 
    Output: approximated (denoised) data, of the original dim1 x dim2 shape
    ''' 
    if nus is None:
    # estimate nus amount     
        pca_full = utils.PCA(data)
        npc = min(5, pca_full.gavish_donoho_indicator())

        nus_rng = np.geomspace(.25, 1., num=16)
        csm = []	# cosine similarity

        pb = utils.ProgressbarThread(title='Estimating NUS amount... ', maximum=len(nus_rng))

        for i, nus in enumerate(nus_rng):
            d = np.zeros(npc)
            for m in range(10): 
                ws = sampling_schedule(len(data), nus)
                pca = utils.PCA(data[ws])
        
                # calculate the cosine similarity between PC vectors
                A = np.abs(pca_full.pcs[:npc])
                B = np.abs(pca.pcs[:npc])
                d += np.diag(A.dot(B.T))
                
            if pb.stop_event.is_set():
                shutil.rmtree(f'{a.dir}/{a.name}/{str(new_expno)}/')
                sys.exit()
            pb.update(i+1)
            csm.append(d/10)

        pb.close()

        # compute slopes of csm vs nus:
        y = (csm[-1]- csm)[::-1]
        x = nus_rng - nus_rng[0]
        slopes = x.dot(y) / x.dot(x)

        # find where csm drops below 0.9:
        ymean = np.mean(np.outer(x, slopes), axis=-1)
        nus = nus_rng[(ymean > .01).sum()]
        
    else:
        nus = float(nus)
        assert 0. < nus < 1.
        
    pb = utils.ProgressbarThread(title='Denoising in progress', maximum=runs*2)
            
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
    pb.update(runs+1)  
    
    for i in range(1, runs):
        pcs2 += _denoise(data, nus, samples, weighted, mri)
        
        if pb.stop_event.is_set():
            shutil.rmtree(f'{a.dir}/{a.name}/{str(new_expno)}/')
            sys.exit()
        pb.update(runs+i+2) 

    pcs1 /= runs      
    pcs2 /= runs
        
    # calculate outer products between vectors stored in pcs1 and pcs2
    npc = min(len(pcs1), len(pcs2))
    kernels = pcs1[:npc,:,None] @ pcs2[:npc,None,:]

    # calculate Frobenius inner products between each kernel and data
    fscores = np.tensordot(kernels, data.conj())
    
    # filter data
    data_approx = np.tensordot(fscores.conj(), kernels, axes=(0, 0))
    
    pb.close()
    return data_approx.T

def _denoise(data, nus, samples, weighted, mri):
    ''' 
    Repeated data sampling, of size dim1 x nus, followed by PCA, resulting in 
    a series of dim1 x nus matrices of psudo-random PCs. The psudo-random PCs 
    of same rank are used to approximate a full-data PC of that rank. The full-
    data PCs thus approximated ('denoised') are sent to the function denoise() 
    to proceed with approximation of data themselves. 
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

# denoise data
args = parser.parse_args()
data_denoised = denoise(data, nus=args.nus, samples=args.samples, runs=10)

# write denoised data to a ser file  
data_denoised = np.dstack((data_denoised.real, data_denoised.imag))
data_denoised = np.array([x.ravel() for x in data_denoised])  
ds.writeser(data_denoised) 

hsqc.launch('re '+ str(new_expno))
hsqc = dp.getCurrentDataset()

