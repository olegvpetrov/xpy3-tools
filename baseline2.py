"""
This script implements interactive baseline correction on a series of spectra. 

Makes use of two external algorithms: std_distribution [1] for its baseline 
classification part and ws2d [2] for baseline interpolation. 

Includes GUI for interactive optimization of the following parameters: 
 - half_window: the number of data points involved in rolling std calculation;
 - num_std: deviation threshold for assigning data points as baseline;
 - lambda: interpolation smoothing parameter.

[1] Erb, D. (2022) pybaselines: A Python library of algorithms for the baseline 
    correction of experimental data. https://doi.org/10.5281/zenodo.5608581
[2] Whittaker core functionality used in the modape package. 
    https://github.com/WFP-VAM/vam.whittaker 
    
Usage requires Python 3.8+ environment

 - from shell (terminal):  python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of Apr 20, 2024 
 
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
from pybaselines.classification import std_distribution
from pybaselines.utils import optimize_window
from scipy.ndimage import zoom
from scipy.signal import hilbert
from vam.whittaker import ws2d, ws2doptv
from utils import NMRDataSetAttributes, ScrollEventHandler

class BaseLinePlot():

    def __init__(self, x, y):
        self.x = x
        self.y = y

        self.hw = max(optimize_window(self.y), optimize_window(-self.y))
        self.hw = (self.hw / 1).__ceil__()
        self.ns = 2.0   
        self.l = None     
                        
        self.mask = self.get_mask()
        self.baseline, self.l = self.get_baseline() 
        
        # save initial values: 
        self.hw0 = self.hw
        self.ns0 = self.ns
        self.l0 = self.l
        self.mask0 = self.mask
        self.baseline0 = self.baseline
        
    def get_mask(self): 
        ''' Baseline classification function'''
        return std_distribution(self.y, num_std=self.ns, half_window=self.hw)[1]['mask']

    def get_baseline(self):
        ''' Baseline interpolation function'''
        w = np.array(self.mask, dtype='d')
        if self.l is None:
            # run version with V-curve optimization of lambda
            lrange = array.array('d',np.linspace(2,13,111))
            z, l = ws2doptv(self.y, w, lrange)
            return z, np.log10(l)
        else:
            # run basic version
            z = ws2d(self.y, 10**self.l, w)
            return z

class AppWindow():

    def __init__(self, slides, title=None):   
        self.root = tk.Tk()
        self.root.title(title)
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        
        self.slides = slides
        self.index = 0

        self.slide = slides[0]
        self.create_figure()
                
        self.create_ui()
        self.root.mainloop()  
        
    def create_figure(self):        
        self.figure = Figure(figsize=(8.2, 9.2))
        self.ax1 = self.figure.add_subplot(211)
        self.ax2 = self.figure.add_subplot(212, sharey=self.ax1)
        self.figure.subplots_adjust(left=0.1, bottom=0.05, right=0.98, 
                                    top=0.95, wspace=0.4, hspace=0.17)

        self.ax1.plot(self.slide.x[self.slide.mask], self.slide.y[self.slide.mask], 'o', c='orange', ms=4, label='baseline mask')
        self.ax1.plot(self.slide.x, self.slide.y)
        self.ax1.plot(self.slide.x, self.slide.baseline, c='orangered', label='approximation')
        self.ax2.plot(self.slide.x, self.slide.y - self.slide.baseline)

        self.ax1.legend(frameon=False, title=f'# {self.index+1}')
        self.ax1.invert_xaxis()
        self.ax2.invert_xaxis()

        ymax = max(np.abs(self.ax1.get_ylim())) 
        self.ax1.set_ylim(-ymax*0.5, ymax*0.5)
        self.ax2.set_ylim(-ymax*0.5, ymax*0.5)

        self.handler = ScrollEventHandler(self.ax1)
        self.figure.canvas.mpl_connect('scroll_event', self.handler.on_scroll)       
        
    def create_ui(self):
        self.hw = tk.IntVar(value=self.slide.hw)
        self.ns = tk.DoubleVar(value=f'{self.slide.ns:0.2f}')
        self.l = tk.DoubleVar(value=f'{self.slide.l:0.2f}')

        controls = tk.Frame(self.root)
        tk.Label(controls, text="half_window").grid(row=0, column=0, padx=(25,10), sticky='e')
        tk.Label(controls, text="num_std").grid(row=1, column=0, padx=10, sticky='e')
        tk.Label(controls, text="lambda").grid(row=2, column=0, padx=10, sticky='e')
        tk.Label(controls, textvariable=self.hw, width=5).grid(row=0, column=2, sticky='w')
        tk.Label(controls, textvariable=self.ns, width=5).grid(row=1, column=2, sticky='w')
        tk.Label(controls, textvariable=self.l, width=5).grid(row=2, column=2)
               
        tk.Scale(controls, from_=max(int(self.slide.hw*0.1), 3), to=self.slide.hw*4, resolution=1, variable=self.hw, 
            orient='h', length=460, showvalue=0, command=self.update_classification).grid(row=0, column=1, pady=5)
        tk.Scale(controls, from_=max(self.slide.ns*0.5, 0.5), to=max(self.slide.ns*1.5, 3.0), resolution=.01, variable=self.ns, 
            orient='h', length=460, showvalue=0, command=self.update_classification).grid(row=1, column=1, pady=5)
        tk.Scale(controls, from_=max(self.slide.l-2.,0.1), to=self.slide.l+3., resolution=.01, variable=self.l, 
            orient='h', length=460, showvalue=0, command=self.update_interpolation).grid(row=2, column=1, pady=5)
        
        ttk.Button(controls, text="Reset", command=self.reset).grid(row=0, column=3, padx=(15,3), pady=(5,3), sticky='e')
        ttk.Button(controls, text="<", command=self.previous).grid(row=1, column=3, padx=(15,3), pady=3, sticky='e') 

        ttk.Button(controls, text="Save", command=self.save).grid(row=0, column=4, padx=(3,15), pady=(5,3), sticky='e')
        ttk.Button(controls, text=">", command=self.next).grid(row=1, column=4, padx=(3, 15), sticky='e')
                
        self.canvas = FigureCanvasTkAgg(self.figure, self.root)
        canvas_widget = self.canvas.get_tk_widget()
        NavigationToolbar2Tk(self.canvas, self.root)
        controls.pack(side='bottom', anchor='w')
        canvas_widget.pack(fill='both', expand=True) 

    def update_classification(self, val):
        self.slide.hw = self.hw.get()
        self.slide.ns = self.ns.get()
        self.slide.mask = self.slide.get_mask()
        self.slide.baseline = self.slide.get_baseline() 
        self.redraw()
           
    def update_interpolation(self, val):
        self.slide.l = self.l.get()
        self.slide.baseline = self.slide.get_baseline() 
        self.redraw()
        
    def reset(self):        
        self.slide.hw = self.slide.hw0
        self.slide.ns = self.slide.ns0
        self.slide.l = self.slide.l0
        self.slide.mask = self.slide.mask0
        self.slide.baseline = self.slide.baseline0
        self.redraw()

    def redraw(self):
        self.hw.set(self.slide.hw) 
        self.ns.set(f'{self.slide.ns:0.2f}')
        self.l.set(f'{self.slide.l:0.2f}')
        self.ax1.lines[0].set_xdata(self.slide.x[self.slide.mask])
        self.ax1.lines[0].set_ydata(self.slide.y[self.slide.mask])
        self.ax1.lines[1].set_ydata(self.slide.y)
        self.ax1.lines[2].set_ydata(self.slide.baseline)
        self.ax2.lines[0].set_ydata(self.slide.y - self.slide.baseline)
        self.ax1.legend(title = f'# {self.index+1}')
        self.canvas.draw_idle()
        
    def next(self):
        if self.index == len(self.slides)-1:
            return 
        self.index += 1
        self.slide = self.slides[self.index]  
        self.redraw()    
       
    def previous(self):
        if self.index == 0:
            return 
        self.index -= 1
        self.slide = self.slides[self.index]  
        self.redraw()     
        
    def save(self):
        for i, slide in enumerate(self.slides):
            if len(pdata[i]) > len(slide.baseline):
                # zoom baseline to match full-size pdata
                slide.baseline = zoom(slide.baseline, len(pdata[i])/len(slide.baseline))

            pdata[i] -= slide.baseline 
            if hsqc.getPar('AQ_mod') == '1': 
                pdata_2ii[i] += np.imag(hilbert(slide.baseline))
         
        # write baseline-corrected spectra in a 2rr file
        dta.writespect2d(pdata)
        dta.writespect2d(pdata_2ii, name="2ii")
        self.close()

    def close(self):
        self.root.destroy()
        
# create TopSpin interface object
top = Topspin()
dp = top.getDataProvider()

# get the dataset visible on the screen as a NMRDataSet object
hsqc = dp.getCurrentDataset()
if not hsqc or hsqc.getDimension() != 2: 
    top.showError("Please open a 2D data set ")
    sys.exit(0)
    
# read 1D complex processed bruker data
a = NMRDataSetAttributes(hsqc)
if not a.is_processed():
    top.showError("Processed data not found.")
    sys.exit(0)

dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])
pdata = dta.readspect2d()
pdata_2ii = dta.readspect2d(name="2ii")

# downsize data for smoother graphics
skip = (len(pdata[0])/12288).__ceil__()
x = dta.getprocxppm()[::skip]
y = pdata[:,::skip] 

# start baseline correction... 
slides = []
for i in y:
    if i.any(): 
        slides.append(BaseLinePlot(x, i))
        
# show corrected spectra, write into a 2rr file
slideshow = AppWindow(slides, title=a.get_ts_title())

# refresh TopSpin window
top.getDisplay().show(hsqc, newWindow=False) 
