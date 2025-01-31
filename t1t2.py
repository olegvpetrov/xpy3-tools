"""
This script integrates a series of spectra measured at varying delays (pseudo 2D data) over an interactively selected region. The integral vs delay plot is analyzed against exponential relaxation models. The functionality is therefore similar to 
Topspin's T1T2 toolbox but with some differences:

1) Works with processed data only, not with FIDs.
2) Cannot plot the intensity of the biggest peak in the integration region, only the integral itself.
3) Includes a noise filter based on SVD threshold method (optional).
4) Subtracts baseline contribution from the integrals (optional).
5) Relaxation models are limited to mono-, two-, and stretched exponential.

The relaxation curve is saved to a text file <EXPNO>/pdata/1/ct1t2.001, along with 
a best-fit curve.

Usage requires Python 3.8+ environment

 - from shell (terminal):  python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of Jan 30, 2025 
 
"""
import brukerIO
import numpy as np
import sys

import tkinter as tk
import tkinter.ttk as ttk

from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from math import gamma
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.offsetbox import AnchoredText
from matplotlib.widgets import SpanSelector
from pybaselines.classification import std_distribution
from pybaselines.utils import optimize_window
from scipy.optimize import minimize
from vam.whittaker import ws2doptv
from utils import NMRDataSetAttributes, ScrollEventHandler, PCA, ProgressbarToplevel

def onselect(xmin, xmax):
    global intrng
   
    if xmax == xmin:
        return
    
    # set selection limits for integration
    indmin, indmax = np.searchsorted(ppm[::-1], (xmin, xmax))
    indmin, indmax = max(1, indmin), min(indmax, len(ppm)-1)
    intrng = (int(len(ppm)-indmax), int(len(ppm)-indmin))
    
    # draw a vertical span (rectangle) across ax1
    if ax1.patches: 
        ax1.patches.clear() # only one selection at a time
    span = ax1.axes.axvspan(xmin=ppm[-indmin], xmax=ppm[-indmax], alpha=0.5, color='tab:blue', lw=0)
    
    annotation1.txt.set_text('') 
    integrate()

def integrate():    
    if not intrng:
        return

    if noise_filter.get():
        pca = PCA( rows[:,slice(*intrng)])
        npc = pca.malinowski_indicator()         #npc = pca.gavish_donoho_indicator()
        integrals = pca.evaluate(n=npc).real
        
    else:
        integrals = np.sum(rows[:,slice(*intrng)].real, axis=-1)     
        
    if baseline.get():
        pb = ProgressbarToplevel(root, title='Baseline correction in progress', maximum=len(rows))

       # search for mask (baseline classification part)
        the_row = rows[np.argmax(np.sum(np.abs(rows), axis=-1))].real
        hw = max(optimize_window(the_row), optimize_window(-the_row)).__ceil__()
        mask = std_distribution(the_row, num_std=2.0, half_window=hw)[1]['mask']
    
       # given the mask, interpolate baseline for each row
        baselines = np.zeros(rows.shape)
        for i, row in enumerate(rows.real):
            baselines[i] = get_baseline(row, mask)
            pb.update(i+1) 
            
       # calculate total intensities at the mask positions of rows
        base_mask = np.sum(rows.real[:, mask], axis=-1) 
        
       # estimate baseline intensities in intrng 
        base_intrng = np.sum(baselines[:,slice(*intrng)], axis=-1)
        
       # fit intensity curve measured at mask positions to one for intrng
        def fun(p, c1, c2):
            y = c2 - (p[1] * c1 + p[0])
            return np.sum(y * y)
        
        res = minimize(fun, [0., 1.], args=(base_mask, base_intrng), method='powell')
        base_integrals = base_mask * res.x[1] + res.x[0]
        
       # subtruct thus estimated baseline contribution from total intensity at intrng
        integrals -= base_integrals 
        pb.close() 
        
    # clear the results of fitting, if any
    ax2.lines[1].set_ydata([None])    
    annotation2.txt.set_text('')
    
    # plot integrals    
    ax2.lines[0].set_xdata(xdata)
    ax2.lines[0].set_ydata(integrals)
    ax2.legend([ax2.lines[0]], ['integral'])
#    ax2.legend_ = None

    ax2.relim()
    ax2.autoscale()
    
    canvas.draw() 
    if vdlist_flag:
        fit_btn.configure(state = tk.NORMAL) 
    save_btn.configure(state = tk.NORMAL) 
       
def set_scale():
    if not intrng:
        return
        
    if semilogx.get():
        ax2.lines[0].set_xdata(xdata)
        ax2.lines[1].set_xdata(xdata)
        ax2.set_xscale('log')
    else:
        ax2.lines[0].set_xdata(xdata)
        ax2.lines[1].set_xdata(xdata)
        ax2.set_xscale('linear')   
        
    if semilogy.get():
        ax2.set_yscale('log')
    else:
        ax2.set_yscale('linear')  
        
    ax2.relim()
    ax2.autoscale()
    canvas.draw() 

def fitfunc():
    if not intrng:
        return
        
    # get best-fit parameters
    xdata = ax2.lines[0].get_xdata()
    ydata = ax2.lines[0].get_ydata()
    p = _fitfunc(xdata, ydata)
    
    # plot best-fit function
    ax2.lines[1].set_xdata(xdata)
    ax2.lines[1].set_ydata(func(p, xdata))

    if len(p) == 3: # Mono-exponential
        annotation2.txt.set_text(f'\n\
            Model: y = A1*exp(-x/t1) + y0\n\n\
            y0          {p[2]:.1f}\n\
            A1          {p[0]:.1f}\n\
            t1           {1e3/p[1]:0.2f} ms')
    elif len(p) == 4: # Stretched exponential
        annotation2.txt.set_text(f'\n\
            Model: y = A1*exp(-(x/t1)^beta) + y0\n\n\
                         y0           {p[2]:.1f}\n\
                         A1           {p[0]:.1f}\n\
                         beta        {p[3]:.3f}\n\
                       <t1>        {1e3*(gamma(1/p[3])/p[1])/p[3]:0.2f} ms')
    elif len(p) == 5: # Two-exponential
        annotation2.txt.set_text(f'\n\
            Model: y = A1*exp(-x/t1) + A2*exp(-x/t2) + y0\n\n\
                         y0          {p[2]:.1f}\n\
                         A1          {p[0]:.1f}\n\
                         t1           {1e3/p[1]:0.2f} ms\n\
                         A2          {p[3]:.1f}\n\
                         t2           {1e3/p[4]:0.2f} ms')               
    
    ax2.legend([ax2.lines[0], ax2.lines[1] ], ['integral', 'fit'])
    canvas.draw()
    
def _fitfunc(x, y) -> list: 
    '''
    Fit function to data using minimize() method from scipy.optimize library
     
    :return: best-fit parameters for func
    '''
    # initialize function parameters
    try: 
       # solve Az = b
        A = np.array((np.ones(x.size//2), x[0:x.size//2]))    
        b = np.log(np.abs(y[0:x.size//2] -np.mean(y[-2:])))     
        z = np.linalg.lstsq(A.T, b, rcond=None)
        ampl = np.exp(z[0][0])
        rate = -z[0][1]
    except:
        ampl = np.abs(y[0])
        rate = 1/(x[-1]-x[0])
    offset = y[-1]
    beta = 1.0
    if y[0] < y[-1]:
        ampl *= -1

    if model.get() == 'mono':
        p0 = [ampl, rate, offset] 
        bnds = ((None, None), (0, None), (None, None)) 
    elif model.get() == 'stretched':
        p0 = [ampl, rate, offset, beta] 
        bnds = ((None, None), (0, None), (None, None), (0, 1)) 
    elif model.get() == 'two':
        p0 = [.5*ampl, rate, offset, .5*ampl, rate] 
        bnds = ((None, None), (0, None), (None, None), (None, None), (0, None))     
             
    return minimize(objective, p0, args=(x, y), method='powell', bounds=bnds).x

def objective(p, x, y):
    residual = y - func(p, x)   
    return np.sum(residual * residual)        

def func(p, x):
    if len(p) == 3:		
        return p[0] * np.exp(-p[1]*x) + p[2]
    elif len(p) == 4:		
        return p[0] * np.exp(-(p[1]*x) ** p[3]) + p[2]
    elif len(p) == 5:		
        return p[0] * np.exp(-p[1]*x) + p[2] + p[3] * np.exp(-p[4]*x) 
        
def save():
    xdata = ax2.lines[0].get_xdata()
    ydata = ax2.lines[0].get_ydata()
    yfit = ax2.lines[1].get_ydata()

    if np.any(yfit):
        header = f'tau       integral       model       integration region: {ppm[intrng[1]]:.2f}-{ppm[intrng[0]]:.2f} ppm'
        data = np.column_stack((xdata, ydata, yfit))
    else:
        header = f'tau       integral       integration region: {ppm[intrng[1]]:.2f}-{ppm[intrng[0]]:.2f} ppm'
        data = np.column_stack((xdata, ydata))    

    path = hsqc.getIdentifier() + '/ct1t2.001'
    np.savetxt(path, data, fmt='%.4e', header=header, comments='', newline='\r\n')    

    root.destroy()

def get_baseline(y, mask): 
    ''' Baseline interpolation function'''
    w = np.array(mask, dtype='d')  
    lrange = array.array('d', np.linspace(2,13,111))
    z, _ = ws2doptv(y, w, lrange)
    return z
    
# create TopSpin interface object
top = Topspin()
dp = top.getDataProvider()

# get the dataset visible on the screen as a NMRDataSet object...
hsqc = dp.getCurrentDataset()
if not hsqc or hsqc.getDimension() != 2: 
    top.showError("Please open a 2D data set ")
    sys.exit(0)

attr = NMRDataSetAttributes(hsqc)
if not attr.is_processed():
    top.showError("Processed data not found.")
    sys.exit(0)
   
# and as a brukerIO object
dta = brukerIO.dataset([attr.name, attr.expno, attr.procno, attr.dir])
ppm = dta.getprocxppm() 

# read 2D data, remove zero-valued rows if any
rows = dta.readspect2d() + 1j*dta.readspect2d(name="2ii")
rows = rows[np.any(rows, axis=1)]

# read variable delay list (VDLIST):
try:
    xdata = np.loadtxt(str(attr._path.parents[1]) + '/vdlist')
    vdlist_flag = True
except:
    top.showMessage("Cannot retrieve variable delay list (VDLIST). Will be using a numerical range instead.")
    xdata = np.arange(len(rows))
    vdlist_flag = False
    
# make sure xdata is sorted
if not np.all(xdata[:-1] <= xdata[1:]):
    rows = rows[xdata.argsort()]
    xdata.sort()
    
# define a global variable for integration region
intrng = ()
    
# set up a matplotlib figure
fig = Figure(figsize=(1.25*6.4, 4.4*4.8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
fig.subplots_adjust(left=0.09, bottom=0.05, right=0.98, top=0.95, wspace=0.4, hspace=0.17) 

ax1.grid(alpha = 0.25)
ax2.grid(alpha = 0.25)         
ax1.invert_xaxis()

ax1.set_xlabel('ppm')
ax1.xaxis.set_label_coords(.98, -0.026) 
ax2.set_xlabel('sec')
ax2.xaxis.set_label_coords(.98, -0.026) 

annotation1 = AnchoredText('Drag mouse over the top graph to select an interval.', loc='center',  frameon=False)
ax2.add_artist(annotation1) 

annotation2 = AnchoredText('', loc=5,  frameon=False)
ax2.add_artist(annotation2)

for y in rows:
    ax1.plot(ppm, y.real)
ax2.plot([None], [None], '*')
ax2.plot([None], [None], '-')
    
# raise tkinter window, await user's commands             
root = tk.Tk()
root.wm_title(NMRDataSetAttributes(hsqc).get_ts_title())

controls = tk.Frame(root)
noise_filter = tk.IntVar(value=0)  
baseline = tk.IntVar(value=0)
semilogx = tk.IntVar(value=0)
semilogy = tk.IntVar(value=0)
model = tk.StringVar(value='mono')      

tk.Checkbutton(controls, text=" Baseline corrector", variable=baseline, command=integrate).grid(row=0, column=0, padx=(15,0), pady=(10,0), sticky='w')
tk.Checkbutton(controls, text=" Noise filter", variable=noise_filter, command=integrate).grid(row=0, column=1, padx=(5,0), pady=(10,0), sticky='w')

ttk.Separator(controls, orient=tk.VERTICAL).grid(row=0, column=2, rowspan=1, padx=(5,5), pady=(10,0), sticky='ns')

tk.Label(controls, text='Exponential model:').grid(row=0, column=3, pady=(10,0), sticky='w')
tk.Radiobutton(controls, text='Mono', variable=model, value='mono').grid(row=0, column=3, padx=(130,0), pady=(10,0), sticky='w')  
tk.Radiobutton(controls, text='Two', variable=model, value='two').grid(row=0, column=3, padx=(195,0), pady=(10,0), sticky='w') 
tk.Radiobutton(controls, text='Stretched', variable=model, value='stretched').grid(row=0, column=3, padx=(250,0), pady=(10,0), sticky='w')

# TODO
ilt_radio = tk.Radiobutton(controls, text='ILT', variable=model, value='ilt')
ilt_radio.grid(row=0, column=3, padx=(342,0), pady=(10,0), sticky='w')
ilt_radio.configure(state = tk.DISABLED)

fit_btn = ttk.Button(controls, text='Fit Me!', command=fitfunc)
fit_btn.grid(row=0, column=4, padx=(20,0), pady=(10,0), sticky='w')
fit_btn.configure(state = tk.DISABLED)

tk.Label(controls, text='Log-scale:').grid(row=1, column=0, padx=(25,0), pady=(5,0), sticky='w')
tk.Checkbutton(controls, text=" X Axis", variable=semilogx, command=set_scale).grid(row=1, column=0, padx=(95,0), pady=(5,0), sticky='w')
tk.Checkbutton(controls, text=" Y Axis", variable=semilogy, command=set_scale).grid(row=1, column=1, padx=(5,0), pady=(5,0), sticky='w')

save_btn = ttk.Button(controls, text='Save& Quit', command=save)
save_btn.grid(row=1, column=4, padx=(20,0), pady=(5,0), sticky='w')
save_btn.configure(state = tk.DISABLED)

canvas = FigureCanvasTkAgg(fig, root)
canvas_widget = canvas.get_tk_widget()
        
handler = ScrollEventHandler(ax1)
canvas.mpl_connect('scroll_event', handler.on_scroll)
      
NavigationToolbar2Tk(canvas, root)
controls.pack(side='bottom', anchor='w')
canvas_widget.pack(side='top', fill='both', expand=True) 
        
# use span selector
span_selector = SpanSelector(ax1, onselect, 'horizontal', useblit=True, props=dict(alpha=0.5, facecolor='tab:blue'))
#canvas.mpl_connect('key_press_event', span_selector)

root.mainloop()    
