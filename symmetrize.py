"""
This script equalizes opposite spinning sidebands in a quadupolar MAS spectrum 
while keeping the centerband intact. It may be useful for quadrupolar lineshape 
analysis when a first-order interaction model is being applied which is known 
to rely on symmetrical distribution of sideband intensity about the centerband. 

Usage requires Python 3.8+ environment

 - from shell (terminal):     python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of May 17, 2023 
 
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
from scipy.signal import find_peaks, peak_prominences
from statistics import mode
from utils import NMRDataSetAttributes, ScrollEventHandler 

def get_cb_range(top:Topspin, ds:NMRDataSet):
    '''
    Generates an integrate regions file intrng using the Topspin command 'int auto', 
    the region with highest integral value is to be assigned to the centerband region.
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
    Measure spinning sideband period through auto-correlation
        
    :data: input spectrum with sidebands
    :cb: centerband indices 
        
    :return: spinning sideband period, in points
    '''
    # zeroize the centerband 
    da_ta = np.copy(data)
    da_ta[cb] = 0.

    # compute auto-correlation 
    selfcor = np.correlate(da_ta, da_ta, mode='same')
    selfcor = selfcor[len(selfcor)//2:]

    # find most prominent peaks of auto-correlation 
    peaks = find_peaks(selfcor)[0]
    prominences = peak_prominences(selfcor, peaks)[0]
        
    main_peaks = [0]
    while peaks.any():
        i = np.argmax(prominences)
        main_peaks = np.r_[main_peaks, peaks[i]]
        prominences = prominences[i+1:]
        peaks = peaks[i+1:]
        
    # assign the mode of the main peak interval to period
    period = mode( np.diff(main_peaks) )
    return int(period)

def symmetrize(data, cb, period, mode='mean'):
    ''' 
    :data:   input spectrum with sidebands
    :cb:     centerband indices 
    :period: sideband period, in points
    :mode:   'mean', 'maximum' 
    
    :return: a copy of spectrum with opossite sidebands of equal intensity  
    
    In mode 'maximum' the opposite sidebands are compared one by one and that
    with a greater integral replaces the other. In mode 'mean' (default) it is
    their average that is placed on both sides of the centerband. In both cases
    the centerband is kept intact.  
    '''        
    data_sym = np.copy(data)
    i = 0
    while True:
        idx1 = range(cb[-1]-(i+2)*period, cb[0]-i*period)
        idx2 = range(cb[-1]+i*period, cb[0]+(i+2)*period)

        if idx1[0] <= 0 or idx2[-1] >= len(data):
            break
        
        if mode == 'mean':
            suma = 0.5*(data[idx1] + data[idx2])
            data_sym[idx1] = suma
            data_sym[idx2] = suma 
     
        elif mode == 'maximum':
            if np.sum(data[idx1]) > np.sum(data[idx2]):
                data_sym[idx1] = data[idx1]
                data_sym[idx2] = data[idx1]
            else:
                data_sym[idx1] = data[idx2]
                data_sym[idx2] = data[idx2] 
        i += 1                   
    return data_sym
    
def update(_):   
    global pdata_sym     
    cb = np.arange(int(scale_cb1.get()), int(scale_cb2.get()))
    pdata_sym = symmetrize(pdata, cb, period)
    figure.axes[0].lines[2].set_xdata(x[cb] )
    figure.axes[0].lines[2].set_ydata(pdata[cb])
    figure.axes[0].lines[3].set_xdata([x[cb[0]], x[cb[-1]]])
    figure.axes[0].lines[3].set_ydata([pdata[cb[0]], pdata[cb[-1]]])
    figure.axes[1].lines[0].set_ydata(pdata_sym) 
    canvas.draw_idle()
         
def reset():
    cb = cb_copy
    scale_cb1.set(cb[0])
    scale_cb2.set(cb[-1])
    canvas.draw_idle()  
        
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
    top.showError('Open a 1D data set ')
    sys.exit(0) 

# make it a brukerIO.dataset object for extended I/O capabilities :
a = NMRDataSetAttributes(proton)
dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])

# read processed data
pdata = dta.readspect1dri()[0]

# get centerband index range
cb = get_cb_range(top, proton)
cb_copy = np.copy(cb)

# measure sideband period
period = get_sb_period(pdata, cb)

# equalize oposite sideband intensities
pdata_sym = symmetrize(pdata, cb, period)

# checkpoint
figure = Figure(figsize=(8.2, 9.2))
ax1 = figure.add_subplot(211)
ax2 = figure.add_subplot(212, sharey=ax1)
figure.subplots_adjust(left=0.1, bottom=0.05, right=0.95, top=0.97, wspace=0.4, hspace=0.17)

# mark sideband positions with asterisk
w = period // 5
c = np.argmax(pdata)

ma = np.arange(len(pdata)) // period
ma = ma[w:] > ma[:-w]
ma = ma[:-(len(ma) % period)]
ma = np.roll(ma, c - period + w//2)

sb,_ = find_peaks(pdata[:len(ma)] * ma, distance=period*0.9) 
sb = np.delete(sb, np.where(sb==c))

# plot
xranges = proton.getSpecDataPoints()['physicalRanges']
x = np.linspace(xranges[-1]['start'], xranges[-1]['end'], len(pdata))

ax1.plot(x, pdata)
ax1.plot(x[sb], pdata[sb], '*', c='#ff7f0e', label='recognized sidebands')
ax1.plot(x[cb], pdata[cb], c='orangered', label='recognized centerband')
ax1.plot([x[cb[0]], x[cb[-1]]], [pdata[cb[0]], pdata[cb[-1]]], 's', c='orangered')
ax2.plot(x, pdata_sym, label='symmetrized spectrum')

ymax = ax1.get_ylim()[1]
for ax in figure.axes: 
    ax.set_xlabel(xranges[-1]['unit'])
    ax.xaxis.set_label_coords(-0.03,-.04) 
    ax.set_ylim(-ymax*.01, ymax*.1)
    ax.legend(frameon=False)
    ax.invert_xaxis()
    
handler = ScrollEventHandler(ax1)
figure.canvas.mpl_connect('scroll_event', handler.on_scroll)

# raise the window, await user commands
root = tk.Tk()
root.wm_title(a.get_ts_title())

controls = tk.Frame(root)
tk.Label(controls, text='Adjust centerband limits:').grid(row=0, column=1, pady=(10,0))
tk.Label(controls, text='left limit').grid(row=1, column=0, padx=(50,10), sticky='w')
tk.Label(controls, text="right limit").grid(row=2, column=0, padx=10, sticky='e')
        
scale_cb1 = ttk.Scale(controls, value=cb[0], from_=cb[0]+len(cb)//2-period, to=cb[0]+len(cb)//2-1, 
                 orient='h', length=550, command=update)        
scale_cb1.grid(row=1, column=1, pady=5)
scale_cb2 = ttk.Scale(controls, value=cb[-1], from_=cb[-1]-len(cb)//2+1, to=cb[-1]-len(cb)//2+period,
                 orient='h', length=550, command=update)
scale_cb2.grid(row=2, column=1, pady=5)
               
ttk.Button(controls, text="Reset", command=reset).grid(row=1, column=3, padx=20, pady=5, sticky='e')
ttk.Button(controls, text="OK", command=ok).grid(row=2, column=3, padx=20, sticky='e')
      
canvas = FigureCanvasTkAgg(figure, root)
canvas_widget = canvas.get_tk_widget()
NavigationToolbar2Tk(canvas, root)

# care for packing :
controls.pack(side='bottom', anchor='w')
canvas_widget.pack(side='top', fill='both', expand=True) 

ok_flag = False
root.mainloop()

# on closing the window:
if ok_flag:    
    # write symmetrized pdata to disk, refresh TopSpin window:
    dta.writespect1dri(pdata_sym, dta.readspect1dri()[1])
    top.getDisplay().show(proton, newWindow=False)
