"""
This script transforms CPMG FID into a cumulative sum of the individual echoes,
each echo being represented by either its peak intensity or integral. 
The cumulative sum as a function of echo mumber is (optionally) fitted with 
a mono-exponential relaxation model. The results are saved to a text file 
<expno>/pdata/1/ct1t2.001.

Suitable for 1D qcpmg and similar pulse programs with continuous sampling provided 
that individual echoes are separated by zeroes. Having zeroes is necessary because 
the script splits the echo train into individual echoes based on zero positions. 
Enclosed is a pulse program op_qcpmg that outputs a cpmg signal with zeroes between 
individual echoes - just on this purpose.

Fitting to a cumulative sum of the echoes seems more robust than fitting to the 
echoes themselves, as summation decreases point-to-point scatter [J. Phys. Chem. 
Lett. 2016, 7, 1249-1253]. The used formula of fitting applies only to periodic 
pulse sequences where each signal of a series can be expressed as a function of 
preceding signals (hence the word 'recursive' in the script's name). 

Usage requires Python 3.8+ environment

 - from shell (terminal):     python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of June 30, 2023 

"""
import brukerIO
import numpy as np
import sys

import tkinter as tk
import tkinter.ttk as ttk

from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from matplotlib.offsetbox import AnchoredText
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from pathlib import Path
from pybaselines.classification import _find_peak_segments
from scipy.optimize import leastsq
from statistics import mode
from utils import NMRDataSetAttributes, PCA

# The function to be called on changing the radio variable
def select_signal(): 
    global measure

    if var.get() != measure: 
        measure = var.get()        

        # clear the results of fitting, if any
        if annotation.txt.get_text() != '':
            annotation.txt.set_text('')
            figure.axes[1].lines[1].set_ydata([None])
            figure.axes[1].legend([figure.axes[1].lines[0]], ['echo signal'])        

        # replace signals[measure] for a selected signal and update the plots           
        figure.axes[0].lines[2].set_ydata(signals[measure])
        figure.axes[1].lines[0].set_ydata(np.cumsum(signals[measure]))            

        canvas.draw()

# The function to be called anytime 'Fit Me!' is pressed
def fitfunc():
    # get best-fit parameters
    p = _fitfunc(np.cumsum(signals[measure]), te)
    
    # plot best-fit function
    figure.axes[1].lines[1].set_xdata(xdata)
    figure.axes[1].lines[1].set_ydata(func(p, te))         
    figure.axes[1].lines[1].set_linestyle('solid')

    annotation.txt.set_text(f'\n\
        Model: y = A1*exp(-x/t1) + y0\n\n\
        y0          {p[2]:.1f}\n\
        A1          {p[0]:.1f}\n\
        t1           {1e3/p[1]:0.2f} ms')

    figure.axes[1].legend().texts[1].set_text('fit')
    canvas.draw()

def _fitfunc(y, te) -> list: 
    '''
    Fit func to data (x, y) using leastsq() from a SciPy optimization library

    :return: best-fit parameters for func
    '''
    # initialize function parameters
    try: 
       # solve Az = b
        A = np.array((np.ones(x.size//2), x[0:x.size//2]))    
        b = np.log(abs(np.diff(y, prepend=0)[0:x.size//2]))     
        z = np.linalg.lstsq(A.T, b, rcond=None)
        amplitude = np.exp(z[0][0])
        rate = -z[0][1]
    except:
        amplitude = np.abs(y[0])
        rate = 1 / (x[-1]-x[0])
    offset = y[0]
    p0 = [amplitude, rate, offset] 
  

    # run least-squares optimization
    plsq = leastsq(_residual, p0, args=(y, te))
    return plsq[0] # best-fit parameters

def _residual(p, y, te):
    return y - func(p, te)

# The function to fit:
def func(p, te):
    w = np.arange(1, len(xdata) + 1)
    z = np.exp(-p[1] * te)
    return p[0] * (2 * z**(w + 1) - z - 1) / (z - 1) + p[2] * w

# The function to be called anytime 'Save & Quit' is pressed
def save():

    header = f'tau       {measure}'
    data = np.column_stack((xdata, signals[measure]))

    path = proton.getIdentifier() + '/ct1t2.001'
    print(path)
    np.savetxt(path, data, fmt='%.4e', header=header, comments='', newline='\r\n') 

    root.destroy()  

# create TopSpin interface object
top = Topspin()
dp = top.getDataProvider()

# get the dataset visible on the screen as a NMRDataSet object
proton = dp.getCurrentDataset()
if not proton or proton.getDimension() > 1: 
    top.showError("Please switch to a 1D dataset window ")
    sys.exit(0)

# instantiate a brukerIO dataset object for read / write functionality:
a = NMRDataSetAttributes(proton)
dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])

# read fid:  
y = dta.readfidc(rmGRPDLY=False, applyNC=False)
x = dta.getxtime()[::2]

# find echo segments:
starts, ends = _find_peak_segments(y == 0)
period  = mode(ends - starts)

if len(starts) == 1:
    top.showError("Cannot find zero-separated segments in the signal:\nDivision into echoes aborted.")
    sys.exit(0)

# deal with missing echoes, if any:    
idx = (np.diff(starts) > 1.5*period).nonzero()[0] + 1
starts = np.insert(starts, idx, starts[idx-1] + period)

y2 = [y[i:j+1] for i,j in zip(starts, starts + period)]
x2 = [x[i:j+1] for i,j in zip(starts, starts + period)]  
if len(y2[-1]) != len(y2[0]):
    y2 = y2[:-1]
    x2 = x2[:-1]

y2 = np.asarray(y2)
x2 = np.asarray(x2)
xdata = np.mean(x2, axis=1)
te = np.mean(np.diff(xdata))

# evaluate echo signal:
pca = PCA( y2 )
npc = pca.malinowski_indicator()

midpoints = pca.evaluate_at_index(y2.shape[1]//2, n=npc)
integrals = pca.evaluate(n=npc)
integrals *= midpoints[0] / integrals[0]

signals = {'intensity': midpoints.real, 'area': integrals.real}
measure = 'intensity' 

# prepare matplotlib drawables :   
figure = Figure(figsize=(1.25*6.4, 2.2*4.8))
ax1 = figure.add_subplot(211)
ax2 = figure.add_subplot(212)
figure.subplots_adjust(left=0.1, bottom=0.05, right=0.95, 
                       top=0.95, wspace=0.4, hspace=0.14)                                     
                    
ax1.grid(alpha = 0.25)
ax2.grid(alpha = 0.25)                    

sl = np.log10(ax1.get_ylim()[1]).__floor__()
ax1.ticklabel_format(axis='y', scilimits=(sl,sl))

ax2.set_xlabel('s')
ax2.xaxis.set_label_coords(1.02, -0.03)  

# draw it (on a tkinter canvas)
ax1.plot(x, y.real, label='Re')  
ax1.plot(x, y.imag, label='Im')
ax1.plot(xdata, signals[measure], '*', label='echo signal') 

ax2.semilogy(xdata, np.cumsum(signals[measure]), '*', c='g', label='echo signal')  
ax2.semilogy([], [], linestyle='none', c='r', label=' ')  

ax1.legend(frameon=False) 
ax2.legend(frameon=False)    

annotation = AnchoredText('', loc=5,  frameon=False)
ax2.add_artist(annotation) 

# raise a tkinter window, await user's commands
root = tk.Tk()
root.title(a.get_ts_title())

controls = tk.Frame(root)
var = tk.StringVar(value=measure)        
tk.Label(controls, text='Measure of signal:').grid(row=0, column=0, padx=(80,0), pady=(10,0), sticky='we')
tk.Radiobutton(controls, text='Intensity', variable=var, value='intensity', command=select_signal).grid(row=1, column=0, padx=(70,0), sticky='w')
tk.Radiobutton(controls, text='Area', variable=var, value='area', command=select_signal).grid(row=2, column=0, padx=(70,0), sticky='w')
ttk.Button(controls, text='Fit Me!', command=fitfunc).grid(row=0, column=2, padx=(380,0), pady=(10,0), sticky='w')
ttk.Button(controls, text='Save & Quit', command=save).grid(row=0, column=3, padx=5, pady=(10,0), sticky='w')
      
canvas = FigureCanvasTkAgg(figure, root)
canvas_widget = canvas.get_tk_widget()

NavigationToolbar2Tk(canvas, root)        
controls.pack(side='bottom', anchor='w')
canvas_widget.pack(side='top', fill='both', expand=True) 

root.mainloop()
