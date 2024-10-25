"""
Interactive zero-order phase correction in time domain

Rotates a 1D time-domain signal visible under the 'FID' tab by an arbitrary
angle and saves it under a new EXPNO (next to maximum EXPNO for a given NAME).

The phasing in time domain is particularly useful in CPMG experiments where 
real parts of the echo signals are maximized upon quantitation. 

Usage requires Python 3.8+ environment

 - from shell (terminal):     python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of May 25, 2023 
 
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
from utils import NMRDataSetAttributes, ScrollEventHandler

def update_fid(val):
    global fid
    fid = fid0 * np.exp(-1j* np.pi/180 * angle.get()) 
    figure.axes[0].lines[0].set_ydata(fid.real)
    figure.axes[0].lines[1].set_ydata(fid.imag)
    canvas.draw_idle()
           
def reset():
    angle.set(0.) 
    update_fid(angle) 
        
def ok():  
    global ok_flag
    ok_flag = True
    root.destroy()

# create TopSpin interface object
top = Topspin()
dp = top.getDataProvider()

# get the dataset visible on the screen as an NMRDataSet object
proton = dp.getCurrentDataset()
if not proton or proton.getDimension() > 1: 
    top.showError("Open an active 1-D data set ")
    sys.exit(0)

# instantiate a brukerIO dataset object
a = NMRDataSetAttributes(proton)
dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])

'''
TODO: Check whether it is a complex signal
'''

# read fid with the digital filter part removed:  
fid = dta.readfidc(rmGRPDLY=False, applyNC=True)
fid0 = np.copy(fid)

# create a matplotlib figure
figure = Figure(figsize=(8.2, 4.6))
figure.subplots_adjust(left=0.08, bottom=0.1, right=0.98, top=0.95)
ax = figure.add_subplot(111)
                    
sl = np.log10(ax.get_ylim()[1]).__floor__()
ax.ticklabel_format(axis='y', scilimits=(sl,sl))

# plot fid
ax.plot(fid.real, label='Re')  
ax.plot(fid.imag, label='Im')

ax.grid(alpha = 0.25)
ax.legend(frameon=False)

handler = ScrollEventHandler(ax)
figure.canvas.mpl_connect('scroll_event', handler.on_scroll)

# raise a window, await user's commands
root = tk.Tk()
root.wm_title(a.get_ts_title())
controls = tk.Frame(root)

angle = tk.DoubleVar()
zeroim = tk.IntVar()

tk.Label(controls, text='    PH_ref [degree]').grid(row=0, column=0, pady=10, padx=5, sticky='e')
tk.Scale(controls, from_=-180., to=180., resolution=.1, variable=angle, orient='h', 
         length=500, showvalue=0, command=update_fid).grid(row=0, column=1, pady=5)
tk.Label(controls, textvariable=angle, width=5).grid(row=0, column=2, sticky='w')        
ttk.Button(controls, text="Reset", command=reset).grid(row=0, column=4, padx=10, pady=5, sticky='e')
ttk.Button(controls, text="Done", command=ok).grid(row=1, column=4, padx=10, sticky='e')
tk.Checkbutton(controls, text=' Zero Im', variable=zeroim).grid(row=0, column=3, padx=10, sticky='e')
                
canvas = FigureCanvasTkAgg(figure, root)
canvas_widget = canvas.get_tk_widget()

NavigationToolbar2Tk(canvas, root)
controls.pack(side='bottom', anchor='w')
canvas_widget.pack(side='top', fill='both', expand=True) 

ok_flag = False
root.mainloop()

# on closing the window:
if ok_flag:
    # save fid under new EXPNO
    new_expno = NMRDataSetAttributes(proton).get_max_expno() + 1
    proton.launch('wra '+ str(new_expno) + ' y')
    proton.launch('re '+ str(new_expno))
    proton = dp.getCurrentDataset()
    dta = brukerIO.dataset([a.name, str(new_expno), a.procno, a.dir])
    if zeroim.get():
        fid.imag = 0.
    dta.writefidc(fid)

