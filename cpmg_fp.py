"""
This script is intended for qcpmg type pulse programs generating continuously sampled 
echo decays (aka echo trains). The script will collapse the train into one average echo 
and process it into a regular 1D spectrum, following whole-echo processing steps.

Suitable for qcpmg type pulse programs with continuous sampling provided that
the echoes are separated by one or more zero data points. This is necessery because 
the script splits the signal into individual echoes by looking for zero data points. 
Enclosed is a pulse program op_qcpmg that ensures the signal the zero-point separators 
- consider using it instead of qcpmg. 

Usage requires Python 3.8+ environment

 - from shell (terminal):     python3 <script-name>
 - from TopSpin command line: xpy3 <script-name> 

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of July 8, 2023 
 
"""
import brukerIO
import numpy as np
import sys

import tkinter as tk

from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from pathlib import Path
from pybaselines.classification import _find_peak_segments
from scipy.optimize import leastsq, minimize
from statistics import mode
from utils import NMRDataSetAttributes, fft

# The function called on changing the number of echoes in the sum
def update_nech(val):
    global ydata
    ydata = np.mean(y2[:nech.get()], axis=0) 
    figure.axes[1].lines[0].set_ydata(ydata.real)
    figure.axes[1].lines[1].set_ydata(ydata.imag)        
    figure.axes[1].relim()
    figure.axes[1].autoscale_view(scalex=False)
    canvas.draw() 

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
      
# employ brukerIO dataset class for I/O functionality
a = NMRDataSetAttributes(proton)
dta = brukerIO.dataset([a.name, a.expno, a.procno, a.dir])

# read the echo train
y = dta.readfidc(rmGRPDLY=False, applyNC=False)
x = dta.getxtime()[::2]

# split the echo train into individual echoes
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

#y2 = np.array([y[i:j+1] for i,j in zip(starts, ends)])
#x2 = np.array([x[i:j+1] for i,j in zip(starts, ends)])

# average y2 down to one whole echo signal 
ydata = np.mean(y2, axis=0)

# create matplotlib figure 
figure = Figure(figsize=(1.25*6.4, 2.2*4.8))
ax1 = figure.add_subplot(211)
ax2 = figure.add_subplot(212)
figure.subplots_adjust(left=0.1, bottom=0.05, right=0.95, 
                       top=0.95, wspace=0.4, hspace=0.17)
ax1.set_xlabel('s')
ax1.xaxis.set_label_coords(1.02, -0.026) 
ax2.set_xlabel('points')
ax2.xaxis.set_label_coords(1., -0.026) 

sl = np.log10(ax1.get_ylim()[1]).__floor__()
ax1.ticklabel_format(axis='y', scilimits=(sl,sl))
ax2.ticklabel_format(axis='y', scilimits=(sl,sl))

# plot 
ax1.plot(x, y.real, lw=0.9, label='Re')  
ax1.plot(x, y.imag, lw=0.9, label='Im')
ax2.plot(ydata.real, label='Re') 
ax2.plot(ydata.imag, label='Im') 

# add a second x-axis for echo numbering 
ax3 = ax1.twiny()
ax3.set_xlim(ax1.get_xlim())

ax3.set_xlabel('echo num.')
ax3.xaxis.set_label_coords(1.00, 1.035) 

xtick_skip = len(x2) // 8
xtick_locations = np.array(np.mean(x2, axis=1))[::xtick_skip]

ax3.set_xticks(xtick_locations)
ax3.set_xticklabels(np.arange(len(xtick_locations))*xtick_skip)

ax2.grid(alpha = 0.25)
ax3.grid(alpha = 0.25) 

ax1.legend(frameon=False) 
ax2.legend(frameon=False, title='Echo average', loc='upper right')  

# raise the window, await user's commands
root = tk.Tk()
root.title(a.get_ts_title())

controls = tk.Frame(root)
nech = tk.IntVar(value=len(y2))
tk.Label(controls, text='# echoes').grid(row=0, column=0, pady=(20,20), padx=(10,5), sticky='e')
tk.Label(controls, textvariable=nech, width=5).grid(row=0, column=2, ipadx=1, sticky='w')
tk.Scale(controls, from_=1, to=len(y2), variable=nech, orient='h', length=500, showvalue=0, command=update_nech).grid(row=0, column=1, pady=5)
tk.Button(controls, text="OK", command=ok).grid(row=0, column=3, padx=(10,0), sticky='w')
tk.Button(controls, text='Cancel',command=root.destroy).grid(row=0, column=4, padx=5, sticky='w')
                
canvas = FigureCanvasTkAgg(figure, root)
NavigationToolbar2Tk(canvas, root)
controls.pack(side='bottom', anchor='w')
canvas.get_tk_widget().pack(side='top', fill='both', expand=True) 

ok_flag = False
root.mainloop()

# on closing the window 
if ok_flag:
    '''
    Whole-echo acquisition data processing.

    Firstly, swap half-spaces around the echo maximum moving one part at the start,
    and the other to the end of the data array. Next, zero-fill the array at the
    center according to the parameter SI. Then, Fourier-transform it with FCOR=1.
    Finally, apply automatic phase correction by minimizing the difference
    between real and imaginary parts of the spectrum. 

    When processed correctly, the imaginary part is close to zero because of the echo
    symmetry. Regardless of whether the automatic phasing works out, one can always
    phase the spectrum manually, by means of Topspin, first setting PHC0 and PHC1  
    such as to bring the real part of the spectrum (as displayed in the Topspin
    window) to zero and then adding 90 degrees to PHC0 [.ph90].
    
    '''  
    # swap half-spaces around maximum 
    imax = np.argmax(abs(ydata))
    ydata = np.roll(ydata, -imax)
    
    # zero-filling 
    si = int(proton.getPar('SI'))
    if si > len(ydata):
        ydata = np.insert(ydata, -imax, np.zeros(si - len(ydata)))
     
    # Fourier transform 
    yfft = fft(ydata, fcor=1)

    # autophasing 
    pivot = len(yfft) // 2
    xdata = np.arange(-pivot, -pivot + yfft.size) / yfft.size
    skip = max(1, len(yfft) // 1024)

    def fun(p, yfft):
        re = (yfft * np.exp(1.j * (p[0] + (p[1] * xdata[::skip])))).real
        im = (yfft * np.exp(1.j * (p[0] + (p[1] * xdata[::skip])))).imag
        return np.sqrt(np.mean((re - im)**2))

    res = minimize(fun, [0., 0.], args=(yfft[::skip]), method='powell')

    if np.sum(yfft.real) > 0:
        phc0 = res.x[0] - 45*np.pi/180
    else:
        phc0 = res.x[0] + 135*np.pi/180
    phc1 = res.x[1]

    yfft *= np.exp(1.j * (phc0 + (phc1 * xdata )))
    yfft *= 2.0

    # assign phc0 and phc1 values to corresponding Topspin parameters
    proton.setPar('PHC0', f'{phc0*180/np.pi:0.3f}')
    proton.setPar('PHC1', f'{phc1*180/np.pi:0.3f}')
    proton.setPar('status PH_mod', 1)
    
    # write the spectrum to disk:
    dta.writespect1dri(yfft.real[::-1], yfft.imag[::-1])
    top.getDisplay().show(proton, newWindow = False)
