"""
This script picks up MAS rotory echoes from an arbitrarily sampled fid, save them as a new 
(undersampled) fid and processes with efp. The resulting spectrum is analogous to one from 
an actual rotor-synchronized experiment. It has no spinning sidebands and the centerband 
consists of central and satellite transitions on a par [see J. Magn. Reson. 177 (2005) 44-
55)], which makes it useful for quantitation by integration of species with distinct 
chemical shifts (applies to spins-3/2, 5/2, etc.). The spectrum is displayed in a new 
Topspin window, with EXPNO set next to the maximum EXPNO for a given NAME.

Usage requires Python 3.8+ environment

 - from shell (terminal):     python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of May 8, 2023 
 
"""
import brukerIO
import numpy as np
import shutil
import sys

import tkinter as tk
import tkinter.ttk as ttk

from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from pathlib import Path
from scipy import ndimage, signal
from statistics import mode
from utils import NMRDataSetAttributes, ScrollEventHandler, PCA

def get_cb_range(top:Topspin, ds:NMRDataSet):
    '''
    Generates an integrate regions file intrng with the Topspin command 'int auto', 
    the region with highest integral value is assigned to the centerband region.
    '''    
    if not NMRDataSetAttributes(ds).is_processed():
        top.showError("Processed data not found.")
        sys.exit(0)

    pdatadir = ds.getIdentifier()
    path = pdatadir+'/intrng'
    
    # back up the current intrng content
    backup = pdatadir+'/intrng.backup'
    if Path(path).exists():
        shutil.copy(path, backup)
    
    # run auto-integrate to get integral regions
    ds.launch('int auto')
    intrng = np.loadtxt(path, skiprows=1)
    
    # calculate the integral value for each region
    intval = np.zeros(len(intrng))
    for i, r in enumerate(intrng):
        pdata = ds.getSpecDataPoints(physRange=[PhysicalRange(r[0], r[-1])])
        intval[i] = np.sum(pdata['dataPoints'])

    # select the  highest integral region  
    maxrng = intrng[np.argmax(intval)]
    start = ds.getIndexFromPhysical(maxrng[0])
    stop = ds.getIndexFromPhysical(maxrng[-1])

    # restore the backup integrate regions
    if Path(backup).exists():
        shutil.move(backup, path)
    top.getDisplay().show(ds, newWindow=False) 
   
    return np.arange(start, stop)
            
def get_sb_period(data, cb) -> int:
    '''
    Measure spinning sideband period by means of auto-correlation method
    
    :data: input spectrum with sidebands
    :cb: centerband indices 
        
    :return: spinning sideband period (as a number of points)
    '''
    # zeroize the centerband 
    da_ta = np.copy(data)
    da_ta[cb] = 0.

    # compute auto-correlation function
    selfcor = np.correlate(da_ta, da_ta, mode='same')
    selfcor = selfcor[len(selfcor)//2:]

    # find most prominent peaks of auto-correlation function
    peaks = signal.find_peaks(selfcor)[0]
    prominences = signal.peak_prominences(selfcor, peaks)[0]
   
    main_peaks = [0]
    while peaks.any():
        i = np.argmax(prominences)
        main_peaks = np.r_[main_peaks, peaks[i]]
        prominences = prominences[i+1:]
        peaks = peaks[i+1:]
        
    # assign the mode of the main peak interval to period
    period = mode( np.diff(main_peaks) )

    return int(period)             
        
def redraw(): 
    figure.axes[0].lines[1].set_xdata(x[samples])
    figure.axes[0].lines[1].set_ydata(abs(fid)[samples])
    canvas.draw_idle()
    
def update_offset(val):
    global samples
    shift = samples[0] - samples_copy[0]
    samples += (int(val) - shift)
    samples = samples[samples < len(x)] 
    redraw()
   
def update_interval(val):
    global samples
    start = samples[0]    
    samples = np.arange(start, len(x), int(val)) 
    redraw()
        
def reset():
    global samples
    samples = samples_copy
    offset.set(0)
    interval.set(period)
    redraw() 
        
def ok():  
    global ok_flag
    ok_flag = True
    root.destroy()       
            
# create TopSpin interface object
top = Topspin()
dp = top.getDataProvider()

# get the dataset visible on the screen as a NMRDataSet object
proton = dp.getCurrentDataset()
if not proton or proton.getDimension() > 1: 
    top.showError("Please open a 1D data set ")
    sys.exit(0)
    
# get the same dataset as a brukerIO.dataset object
a = NMRDataSetAttributes(proton)
dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])
 
# read 1D processed data as a numpy array    
pdata = dta.readspect1d()

# define centerband limits (index range)
cb = get_cb_range(top, proton)

# estimate the rotor period from spining sidebands
period = get_sb_period(pdata, cb)
period = len(pdata) / period

# switch to a new EXPNO while removing the digital filter part from fid
new_expno = NMRDataSetAttributes(proton).get_max_expno() + 1
proton.launch('convdta '+ str(new_expno))
proton = dp.getCurrentDataset()

# refresh dataset attributes
a = NMRDataSetAttributes(proton)
dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])

# read fid as a numpy array
fid = dta.readfidc(rmGRPDLY=False)

# zoom in for better sampling (x10)
fid = ndimage.zoom(fid, zoom := 10)
period = np.rint(period * zoom).astype(int)

# estimate 1st echo's whereabouts
DW = proton.getPar('DW')
DE = proton.getPar('DE')
start = period - float(DE) / float(DW) * zoom

# pick fid peaks while bound to period
start = int(start)
peaks, _ = signal.find_peaks(abs(fid[start:]), distance=period*0.9)

# sample at constant step starting from the leftmost peak found 
samples = np.arange(peaks[0]+start, len(fid), period) 
samples_copy = np.copy(samples) 

# create a matplotlib figure 
figure = Figure(figsize=(8.2, 4.6))
ax = figure.add_subplot(111)
figure.subplots_adjust(left=0.08, bottom=0.1, right=0.95, top=0.95)

# plot fid and mark the found rotary echoes with asterisks               
x = np.arange(len(fid))
ax.plot(x, abs(fid), label='fid magnitude') 
ax.plot(x[samples], abs(fid)[samples], '*', label='rotary echo positions')

ax.legend(frameon=False)
ax.set_xlabel('points', loc='right')

handler = ScrollEventHandler(ax)
figure.canvas.mpl_connect('scroll_event', handler.on_scroll)
                  
# raise a window, await user commands
root = tk.Tk()
root.wm_title(a.get_ts_title())

controls = tk.Frame(root)
offset = tk.IntVar(value=0)
interval = tk.IntVar(value=period)
        
tk.Label(controls, text="Adjust for echo positions:").grid(row=0, column=1, pady=(5,0))
tk.Label(controls, text="offset").grid(row=1, column=0, padx=(25,10), sticky='e')
tk.Label(controls, textvariable=offset, width=5).grid(row=1, column=2, sticky='w')
tk.Label(controls, text="period").grid(row=2, column=0, padx=10, sticky='e')
tk.Label(controls, textvariable=interval, width=5).grid(row=2, column=2, sticky='w')  
        
tk.Scale(controls, variable=offset, from_=-samples[0]//10, to=samples[0]//10, 
    orient='h', length=540, showvalue=0, command=update_offset).grid(row=1, column=1, pady=5)
tk.Scale(controls, variable=interval, from_=int(period*0.95), to=int(period*1.05),  
    orient='h', length=540, showvalue=0, command=update_interval).grid(row=2, column=1, pady=5)

ttk.Button(controls, text="Reset", command=reset).grid(row=1, column=3, padx=30, pady=5, sticky='e')
ttk.Button(controls, text="OK", command=ok).grid(row=2, column=3, padx=30, sticky='e')
      
canvas = FigureCanvasTkAgg(figure, root)
canvas_widget = canvas.get_tk_widget()

NavigationToolbar2Tk(canvas, root)
controls.pack(side='bottom', anchor='w')
canvas_widget.pack(side='top', fill='both', expand=True) 

ok_flag = False
root.mainloop()

# on closing the window 
if ok_flag:

    # prepend fid's first point to the samples
    samples = np.r_[0, samples]
    
    # replace the original fid with the echo sample 
    dta.writefidc(fid[samples])
    
    # modify relevant status parameters as to suit undersampling
    SFO1 = float(proton.getPar('SFO1'))
    SF = float(proton.getPar('SF'))
    SR = float(proton.getPar('SR'))
    SW = float(proton.getPar('SW'))
    SW /= (period/zoom)
    
    # NMRDataSet.setPar() does not work with status parameters, using brukerIO instead
    dta.writeacqpar('TD', len(samples) * 2)
    dta.writeacqpar('SW', SW)
    
    # process the data
    proton.setPar('OFFSET', 1e6 * (SFO1/SF - 1.) + 0.5 * SW * SFO1/SF)
    proton.setPar('PHC1', 0.0)
    proton.setPar('SR', SR)
    proton.launch('efp')     
    
else:
    # delete the newly created EXPNO, switch to the old one
    proton.launch('close')
    shutil.rmtree(f'{a.dir}/{a.name}/{a.expno}/')
