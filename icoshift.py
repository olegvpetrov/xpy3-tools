"""
Peak alignment in a series of spectra with pyicoshift (github.com/sekro/pyicoshift)
 
----------------------- Comments to original pyicoshift: -----------------------

pyicoshift
from scratch python implementation of the
"icoshift - interval Correlation Optimized shifting" algorithm for peak alignment
Sebastian Krossa 08/2019
NTNU Trondheim
sebastian.krossa@ntnu.no

icoshift by Authors:
    Francesco Savorani - Department of Food Science
                         Quality & Technology - Spectroscopy and Chemometrics group
                         Faculty of Sciences
                         University of Copenhagen - Denmark
    email: frsa@life.ku.dk - www.models.life.ku.dk
    Giorgio Tomasi -     Department of Basic Science and Environment
                         Soil and Environmental Chemistry group
                         Faculty of Life Sciences
                         University of Copenhagen - Denmark
    email: giorgio.tomasi@ec.europa.eu - www.igm.life.ku.dk

--------------------------------------------------------------------------------
    
Usage requires Python 3.8+ environment

 - from shell (terminal):  python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of 1 June 2024 
 
"""
import brukerIO
import numpy as np
import re
import sys

import tkinter as tk
import tkinter.ttk as ttk

from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.offsetbox import AnchoredText
from matplotlib.widgets import SpanSelector
from pybaselines.utils import optimize_window
from pyicoshift import Icoshift
from time import sleep
from utils import NMRDataSetAttributes, ScrollEventHandler

# The function called on span selection
def onselect(xmin, xmax):
    
    indmin, indmax = np.searchsorted(ppm[::-1], (xmin, xmax))
    indmin = max(1, indmin)
    indmax = min(indmax, len(ppm)-1)
    
    # highlight selection
    if spans and align_mode.get() == 'whole': 
        clear_spans()        
    span = ax1.axes.axvspan(xmin=ppm[-indmin], xmax=ppm[-indmax], alpha=0.5, color='tab:blue', lw=0)
    spans.append(span) 
   
    # upgrade list of intervals     
    intervals.append((int(len(ppm)-indmax), int(len(ppm)-indmin))) 
    
    # check if any two intervals completely overlap
    intervals.sort()
    n = 1
    while n < len(intervals):
        if intervals[n][1] <= intervals[n-1][1]:
            del intervals[n]
        n += 1
        
    clear_btn.state(["!disabled"])
    canvas.draw()       

def clear_spans():  
    for span in spans:
        span.remove() 

    spans.clear()  
    intervals.clear() 
    
    clear_btn.state(["disabled"])     
    canvas.draw() 

def align():
    global icosh

    if proton.getDimension() == 1: 
        icosh.signals = signals
    else:
        icosh.signals = abs(signals)
    
    if target.get() == 'user_defined':
        icosh.target = icosh.signals[0]
    else:    
        icosh.target = target.get()  
        
    if intervals == []:
        icosh.inter = [(0, len(ppm)-1)]
        icosh.max_shift = optimize_window(signals[-1])
    else:    
        icosh.inter = intervals      
        icosh.max_shift = max([b - a for a, b in intervals])//2
    
    if align_mode.get() == 'whole':
        icosh.fill_mode = 'nan'
    else:
        icosh.fill_mode = 'adjacent'
        
    # run the actual icoshift
    icosh.run()
    
    # apply found interval shifts to whole signals, if 'whole' set
    if align_mode.get() == 'whole':    
        for i in range(len(signals)):
            signal = icosh.result[i]
            signal = signal[slice(*icosh.inter[0])]        
            nans = np.where(np.isnan(signal))[0]
            if nans.size:
                shift = nans[0]
                if shift == 0:
                    shift = nans[-1]
                else:
                    shift -= len(signal)
                icosh.result[i] = np.roll(icosh.signals[i], shift)

    for i, signal in enumerate(icosh.result):
        if proton.getDimension() == 2: 
            signal[mask[i]] *= -1
        ax2.lines[i].set_ydata(signal) 
        ax2.lines[i].set_xdata(ppm) 

    annotation.txt.set_text('')    
    canvas.draw_idle()
    
def save():
    global dta

    if proton.getDimension() == 1:
        for i, dta in enumerate(dta_list):
            dta.writespect1dri(icosh.result[i], dta.readspect1dri()[1])
        proton.launch('.ret')
        proton.launch('.md')    

    elif proton.getDimension() == 2:
        signals_full[:len(icosh.result)] = icosh.result
        dta.writespect2d(signals_full, name="2rr")        

    root.destroy()
            
# create TopSpin interface object
top = Topspin()
dp = top.getDataProvider()

# create Icoshift object
icosh = Icoshift()

# get the dataset visible on the screen as a NMRDataSet object
proton = dp.getCurrentDataset()
signals, intervals, spans = [], [], []

# read 1D datasets associated with proton via multiple display
if proton.getDimension() == 1:
    proton.launch('.ret')
    proton.launch('.md')

    with open(proton.getIdentifier()+'/assocs') as f:
        for line in f:
            if line.startswith('##$MULTDISP') and ("<>" not in line):
                path_list = line.rstrip().split(" ")[1][1:-1].split('|')
                path_list.insert(0, proton.getIdentifier())
                break

    dta_list = []
    for path in path_list:
        a = NMRDataSetAttributes(top.getDataProvider().getNMRData(path))
        dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])
        
        dta_list.append(dta)
        signals.append(dta.readspect1dri()[0])
        
    lens = [len(l) for l in signals]
    if all(x == lens[0] for x in lens):
        signals = np.array(signals)
    else:    
        top.showError("Spectra must have same dimension")
        sys.exit(0)
        
    ppm = dta_list[0].getprocxppm()

# read processed pseudo-2D data: 
if proton.getDimension() == 2:
    a = NMRDataSetAttributes(proton)
    if not a.is_processed():
        top.showError("Processed data not found.")
        sys.exit(0)
    
    dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])
    ppm = dta.getprocxppm() 

    # read 2D data and remove zero-valued rows if any
    signals_full = dta.readspect2d(name="2rr")
    signals = signals_full[np.any(signals_full, axis=1)]
    
    # prepare for using abs(signals) if negative peaks are present (as e.g. in t1ir)
    mask = signals != abs(signals)
    
# initialize matplotlib figure
fig = Figure(figsize=(1.25*6.4, 4.4*4.8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212, sharex=ax1, sharey=ax1)
fig.subplots_adjust(left=0.09, bottom=0.05, right=0.98, top=0.95, wspace=0.4, hspace=0.17) 

ax1.grid(alpha = 0.25)
ax2.grid(alpha = 0.25)         
ax1.invert_xaxis()

annotation = AnchoredText('', loc='center',  frameon=False)
ax2.add_artist(annotation) 
text = "Drag mouse over the top graph to select intervals."
annotation.txt.set_text(text)

for y in signals:
    ax1.plot(ppm, y)
    ax2.plot([], [])

# raise a window, await user's commands             
root = tk.Tk()
root.wm_title(NMRDataSetAttributes(proton).get_ts_title())

controls = tk.Frame(root)
target = tk.StringVar(value='maxcorr')  
align_mode = tk.StringVar(value='n_intervals')      

tk.Label(controls, text='Target:').grid(row=0, column=0, padx=(30,0), pady=(10,0), sticky='e')
tk.Label(controls, text='Align mode:').grid(row=1, column=0, padx=(30,0), pady=(0,0), sticky='e')

tk.Radiobutton(controls, text='average', variable=target, value='average').grid(row=0, column=1, padx=(0,0), pady=(10,0), sticky='w')
tk.Radiobutton(controls, text='median', variable=target, value='median').grid(row=0, column=2, padx=(0,0), pady=(10,0), sticky='w')
tk.Radiobutton(controls, text='max', variable=target, value='max').grid(row=0, column=3, padx=(0,0), pady=(10,0), sticky='w')
tk.Radiobutton(controls, text='average2', variable=target, value='average2').grid(row=0, column=4, padx=(0,0), pady=(10,0), sticky='w')
tk.Radiobutton(controls, text='maxcorr', variable=target, value='maxcorr').grid(row=0, column=5, padx=(0,0), pady=(10,0), sticky='w')
tk.Radiobutton(controls, text='user_defined', variable=target, value='user_defined').grid(row=0, column=6, padx=(0,0), pady=(10,0), sticky='w')

tk.Radiobutton(controls, text='n_intervals', variable=align_mode, value='n_intervals', command=clear_spans).grid(row=1, column=1, padx=(0,0), sticky='w')        
tk.Radiobutton(controls, text='whole', variable=align_mode, value='whole', command=clear_spans).grid(row=1, column=2, padx=(0,0), sticky='w')
            
clear_btn = ttk.Button(controls, text='Clear', command=clear_spans)
clear_btn.grid(row=0, column=6, padx=(170,0), pady=(10,0), sticky='w')
clear_btn.state(["disabled"])
ttk.Button(controls, text='Align', command=align).grid(row=1, column=6, padx=(170,0), pady=(3,0), sticky='w')
ttk.Button(controls, text='Save', command=save).grid(row=2, column=6, padx=(170,0), pady=(3,0), sticky='w')
      
canvas = FigureCanvasTkAgg(fig, root)
canvas_widget = canvas.get_tk_widget()
        
handler = ScrollEventHandler(ax1)
canvas.mpl_connect('scroll_event', handler.on_scroll)
      
NavigationToolbar2Tk(canvas, root)
controls.pack(side='bottom', anchor='w')
canvas_widget.pack(side='top', fill='both', expand=True) 
        
# use span selector
span_selector = SpanSelector(ax1, onselect, 'horizontal', useblit=True, props=dict(alpha=0.5, facecolor='tab:blue'))
canvas.mpl_connect('key_press_event', span_selector)

root.mainloop()
