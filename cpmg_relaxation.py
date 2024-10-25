"""
This script splits a CPMG FID into individual echo segments and evaluates the echoes
using either peak intensity or integral. The resulting relaxation curve are fitted 
with mono or stretched exponential decay functions. The curve is saved to a text file 
<expno>/pdata/1/ct1t2.001.

Suitable for qcpmg type pulse programs with continuous sampling provided 
that individual echoes are separated by zero data points, as the script splits
the train into individual echoes by looking for zero points in a data array. 
Enclosed is a pulse program op_qcpmg that ensures the signal the zero-point 
separators - consider using it instead of qcpmg.

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
from math import gamma
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from pybaselines.classification import _find_peak_segments
from scipy.optimize import leastsq
from statistics import mode
from utils import NMRDataSetAttributes, PCA

# The function called on changing a radio-button variable var1
def select_signal(): 
    global measure

    if var1.get() != measure: 
        measure = var1.get()

        # clear the results of fitting, if any
        if annotation.txt.get_text() != '':
            annotation.txt.set_text('')
            figure.axes[1].lines[1].set_ydata([None])
            figure.axes[1].legend([figure.axes[1].lines[0]], ['echo signal'])

        # replace ydata for selected signal and update the plots
        figure.axes[0].lines[2].set_ydata(signals[measure])
        figure.axes[1].lines[0].set_ydata(signals[measure])

        canvas.draw()

# The function to be called anytime the 'Fit Me!' button is pressed
def fitfunc():    
    # get best-fit parameters
    p = _fitfunc(xdata, signals[measure])

    # plot best-fit function
    figure.axes[1].lines[1].set_xdata(xdata)
    figure.axes[1].lines[1].set_ydata(func(p, xdata))         
    figure.axes[1].lines[1].set_linestyle('solid')

    if len(p) == 3: 
        annotation.txt.set_text(f'\n\
            Model: y = A1*exp(-x/t1) + y0\n\n\
            y0          {p[2]:.1f}\n\
            A1          {p[0]:.1f}\n\
            t1           {1e3/p[1]:0.2f} ms')
    elif len(p) == 4: 
        annotation.txt.set_text(f'\n\
            Model: y = A1*exp(-(x/t1)^beta) + y0\n\n\
                         y0           {p[2]:.1f}\n\
                         A1           {p[0]:.1f}\n\
                         beta        {p[3]:.3f}\n\
                       <t1>        {1e3*(gamma(1/p[3])/p[1])/p[3]:0.2f} ms')

    figure.axes[1].legend().texts[1].set_text('fit')
    canvas.draw()

def _fitfunc(x, y) -> list: 
    '''
    Fit function to data using leastsq() method from a SciPy optimization library
     
    :return: best-fit parameters for func
    '''
    # initialize function parameters
    try: 
       # solve Az = b
        A = np.array((np.ones(x.size//2), x[0:x.size//2]))    
        b = np.log(np.abs(y[0:x.size//2]))     
        z = np.linalg.lstsq(A.T, b, rcond=None)
        amplitude = np.exp(z[0][0])
        rate = -z[0][1]
    except:
        amplitude = np.abs(y[0])
        rate = 1/(x[-1]-x[0])
    offset = np.min(y)
    beta = 1.0

    if var2.get() == 'mono':
        p0 = [amplitude, rate, offset]   
    elif var2.get() == 'stretched':
        p0 = [amplitude, rate, offset, beta] 

    # run least-squares optimization
    plsq = leastsq(_residual, p0, args=(x, y))

    return plsq[0] # best-fit parameters    

def _residual(p, x, y):
    return y - func(p, x)

def func(p, x):
    if len(p) == 3:		# mono-exponential function
        return p[0] * np.exp(-p[1]*x) + p[2]
    elif len(p) == 4:		# stretched exponential function
        return p[0] * np.exp(-(p[1]*x) ** p[3]) + p[2]

# The function to be called anytime 'Save & Quit' is pressed
def save():

    header = f'tau      {measure}'
    data = np.column_stack((xdata, signals[measure]))

    path = proton.getIdentifier() + '/ct1t2.001'
    np.savetxt(path, data, fmt='%.4e', header=header, newline='\r\n') 

    root.destroy()

# create TopSpin interface object
top = Topspin()
dp = top.getDataProvider()

# get the dataset visible on the screen as a NMRDataSet object
proton = dp.getCurrentDataset()
if not proton or proton.getDimension() > 1: 
    top.showError("Please switch to a 1D data set ")
    sys.exit(0)

# instantiate a brukerIO dataset object
a = NMRDataSetAttributes(proton)
dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])

# read cpmg fid  
y = dta.readfidc(rmGRPDLY=False, applyNC=False)
x = dta.getxtime()[::2]

# find echo segments
starts, ends = _find_peak_segments(y == 0)
period  = mode(ends - starts)

if len(starts) == 1:
    top.showError("Cannot find zero-separated segments in the signal:\nDivision into echoes aborted.")
    sys.exit(0)

# deal with missing echoes, if any    
idx = (np.diff(starts) > 1.5*period).nonzero()[0] + 1
starts = np.insert(starts, idx, starts[idx-1] + period)

y2 = [y[i:j+1] for i,j in zip(starts, starts + period)]
x2 = [x[i:j+1] for i,j in zip(starts, starts + period)]
if len(y2[-1]) != len(y2[0]):
    y2 = y2[:-1]
    x2 = x2[:-1]

y2 = np.asarray(y2)
x2 = np.asarray(x2)

# evaluate echo signal
pca = PCA( y2 )
npc = pca.malinowski_indicator()

midpoints = pca.evaluate_at_index(y2.shape[1]//2, n=npc)
integrals = pca.evaluate(n=npc)
integrals *= midpoints[0] / integrals[0]

signals = {'intensity': midpoints.real, 'area': integrals.real}
measure = 'intensity'

xdata = np.mean(x2, axis=1)
ydata = signals[measure]

# initialize matplotlib drawables
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

# draw (on a tkinter canvas)
ax1.plot(x, y.real, label='Re')  
ax1.plot(x, y.imag, label='Im')
ax1.plot(xdata, ydata, '*', label='echo signal') 

ax2.semilogy(xdata, ydata, '*', c='g', label='echo signal')  
ax2.semilogy([], [], linestyle='none', c='r', label=' ')  

ax1.legend(frameon=False) 
ax2.legend(frameon=False)    

annotation = AnchoredText('', loc=5,  frameon=False)
ax2.add_artist(annotation) 

# raise a window, await user's commands
root = tk.Tk()
root.wm_title(a.get_ts_title())

controls = tk.Frame(root)
var1 = tk.StringVar(value=measure)  
var2 = tk.StringVar(value='mono')      

tk.Label(controls, text='Measure of signal:').grid(row=0, column=0, padx=(50,0), pady=(10,0), sticky='we')
tk.Label(controls, text='Exponential decay model:').grid(row=0, column=1, padx=(50,0), pady=(10,0), sticky='we')

tk.Radiobutton(controls, text='Intensity', variable=var1, value='intensity', command=select_signal).grid(row=1, column=0, padx=(41,0), sticky='w')
tk.Radiobutton(controls, text='Area', variable=var1, value='area', command=select_signal).grid(row=2, column=0, padx=(41,0), sticky='w')
tk.Radiobutton(controls, text='Mono', variable=var2, value='mono').grid(row=1, column=1, padx=(41,0), sticky='w')        
tk.Radiobutton(controls, text='Stretched', variable=var2, value='stretched').grid(row=2, column=1, padx=(41,0), sticky='w')

ttk.Button(controls, text='Fit Me!', command=fitfunc).grid(row=0, column=2, padx=(190,0), pady=(10,0), sticky='w')
ttk.Button(controls, text='Save & Quit', command=save).grid(row=0, column=3, padx=5, pady=(10,0), sticky='w')

canvas = FigureCanvasTkAgg(figure, root)
canvas_widget = canvas.get_tk_widget()

NavigationToolbar2Tk(canvas, root)
controls.pack(side='bottom', anchor='w')
canvas_widget.pack(side='top', fill='both', expand=True) 

root.mainloop()
