"""
This script implements interactive baseline correction in 1D. 

Makes use of two external algorithms: std_distribution [1] for the baseline 
classification part and ws2d [2] for the baseline interpolation. 

Includes GUI for interactive optimization of following parameters: 
 - half_window: the number of data points involved in rolling std calculation;
 - num_std: deviation threshold for ascribing data points to baseline;
 - lambda: interpolation smoothing parameter.

[1] Erb, D. (2022) pybaselines: A Python library of algorithms for the baseline 
    correction of experimental data. https://doi.org/10.5281/zenodo.5608581
[2] Whittaker core functionality used in the modape package. 
    https://github.com/WFP-VAM/vam.whittaker 
    
Usage requires Python 3.8+ environment

 - from shell (terminal):  python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
first version on March 10, 2023 

last modification on December 4, 2024:
 -- applies to both 1r and 1i files only if their modification times coincide 
 -- applies even if 1i file is missing (e.g. to magnitute spectra)
 
"""
import brukerIO
import numpy as np
import sys

import tkinter as tk
import tkinter.ttk as ttk

from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from pathlib import Path
from pybaselines.classification import std_distribution
from pybaselines.utils import optimize_window
from scipy.ndimage import zoom
from scipy.signal import hilbert
from utils import NMRDataSetAttributes, ScrollEventHandler
from vam.whittaker import ws2d, ws2doptv

def get_mask(y, ns, hw): 
     ''' Baseline classification function'''
     mask = std_distribution(y, num_std=ns, half_window=hw)[1]['mask'] 
     return mask

def get_baseline(y, mask, l=None): 
    ''' Baseline interpolation function'''
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

def redraw():    
    figure.axes[0].lines[0].set_xdata(x[mask])
    figure.axes[0].lines[0].set_ydata(y[mask])
    figure.axes[0].lines[2].set_ydata(baseline)
    figure.axes[1].lines[0].set_ydata(y - baseline)
    canvas.draw_idle()  

def update(val):
    global mask, baseline
    mask = get_mask(y, ns.get(), hw.get())
    baseline = get_baseline(y, mask, l=l.get()) 
    redraw()
        
def reset():
    hw.set(hw0) 
    ns.set(f'{ns0:0.2f}')
    l.set(f'{l0:0.2f}') 
    update(1)
    canvas.draw_idle()  
        
def ok():
    global baseline
    if len(pdata) > len(y):
        # zoom baseline to match full-size pdata
        baseline = zoom(baseline, len(pdata)/len(y)) 
               
    # write baseline-corrected pdata in 1r and 1i files
    if not Path(proton.getIdentifier() + '/1i').exists():
        dta.writespect1dri(pdata - baseline, -np.imag(hilbert(baseline))) 
        Path.unlink((proton.getIdentifier() + '/1i')) 
    elif int(Path(proton.getIdentifier() + '/1r').stat().st_mtime) == int(Path(proton.getIdentifier() + '/1i').stat().st_mtime):
        dta.writespect1dri(pdata - baseline, dta.readspect1d("1i") + np.imag(hilbert(baseline)))
    else:
        dta.writespect1dri(pdata - baseline, dta.readspect1d("1i"))
        
    # refresh TopSpin window
    top.getDisplay().show(proton, newWindow=False) 
    root.destroy()
        
# create TopSpin interface object
top = Topspin()
dp = top.getDataProvider()

# get the dataset visible on the screen as a NMRDataSet object
proton = dp.getCurrentDataset()
if not proton or proton.getDimension() > 1: 
    top.showError("Please open a 1D data set ")
    sys.exit(0)
       
# read 1D complex processed bruker data: 
a = NMRDataSetAttributes(proton)
if not a.is_processed():
    top.showError("Processed data not found.")
    sys.exit(0)
    
dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])
#pdata = dta.readspect1dri()
pdata = dta.readspect1d("1r")

# downsize data for smoother graphics
skip = (len(pdata)/12288).__ceil__()
x = dta.getprocxppm()[::skip]

#if Path(proton.getIdentifier() + '/1i').exists()):
y = pdata[::skip]

# start baseline fitting... 
hw0 = max(optimize_window(y), optimize_window(-y))
hw0 = (hw0 / 4).__ceil__()
ns0 = 2.0

# classify baseline points 
mask = get_mask(y, ns0, hw0)

# approximate baseline with a Whittaker smoother
baseline, l0 = get_baseline(y, mask)

# plot the result
figure = Figure(figsize=(8.2, 9.2))
ax1 = figure.add_subplot(211)
ax2 = figure.add_subplot(212, sharey=ax1)
figure.subplots_adjust(left=0.1, bottom=0.05, right=0.95, top=0.95, wspace=0.4, hspace=0.17)

ax1.plot(x[mask], y[mask], 'o', c='orange', ms=4, label='baseline mask')
ax1.plot(x, y)
ax1.plot(x, baseline, c='orangered', label='approximation')
ax2.plot(x, y - baseline)

ax1.legend(frameon=False)
ax1.invert_xaxis()
ax2.invert_xaxis()

ymax = ax1.get_ylim()[1]
ax1.set_ylim(-ymax*0.5, ymax*0.5)
ax2.set_ylim(-ymax*0.5, ymax*0.5)

handler = ScrollEventHandler(ax1)
figure.canvas.mpl_connect('scroll_event', handler.on_scroll)

# raise the window, await user commands
root = tk.Tk()
root.title(a.get_ts_title())
       
hw = tk.IntVar(value=hw0)
ns = tk.DoubleVar(value=f'{ns0:0.2f}')
l = tk.DoubleVar(value=f'{l0:0.2f}')

controls = tk.Frame(root)
tk.Label(controls, text="half_window").grid(row=0, column=0, padx=(25,10), sticky='e')
tk.Label(controls, text="num_std").grid(row=1, column=0, padx=10, sticky='e')
tk.Label(controls, text="lambda").grid(row=2, column=0, padx=10, sticky='e')
tk.Label(controls, textvariable=hw, width=5).grid(row=0, column=2, sticky='w')
tk.Label(controls, textvariable=ns, width=5).grid(row=1, column=2, sticky='w')
tk.Label(controls, textvariable=l, width=5).grid(row=2, column=2)
               
tk.Scale(controls, from_=max(hw0//10, 3), to=hw0*4, resolution=1, variable=hw, 
    orient='h', length=500, showvalue=0, command=update).grid(row=0, column=1, pady=5)
tk.Scale(controls, from_=max(ns0//2, 0.5), to=max(ns0*1.5, 3.0), resolution=.01, variable=ns, 
    orient='h', length=500, showvalue=0, command=update).grid(row=1, column=1, pady=5)
tk.Scale(controls, from_=max(l0-2.,0.1), to=l0+3., resolution=.01, variable=l, 
    orient='h', length=500, showvalue=0, command=update).grid(row=2, column=1, pady=5)
        
ttk.Button(controls, text="Reset", command=reset).grid(row=0, column=3, padx=30, pady=5, sticky='e')
ttk.Button(controls, text="Save", command=ok).grid(row=1, column=3, padx=30, sticky='e')
                
canvas = FigureCanvasTkAgg(figure, master=root)
NavigationToolbar2Tk(canvas, root)
controls.pack(side='bottom', anchor='w')
canvas.get_tk_widget().pack(fill='both', expand=True) 

root.mainloop()
