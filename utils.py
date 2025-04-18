'''
Contains common utility functions and classes.

written by Oleg Petrov (oleg.petrov@matfyz.cuni.cz)
version of August 17, 2023

'''
import brukerIO
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import re
import os

import tkinter as tk
import tkinter.ttk as ttk
import threading, _thread

from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from contourpy import contour_generator, LineType
from contourpy.util.mpl_renderer import MplRenderer
from matplotlib import gridspec 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from pathlib import Path

class ContourPlot():
    ''' 2D spectrum display with basic options: zoom, scale, add, ... '''

    figsize = (8.0, 8.0)   
    palette = [['#0000FF', '#008080'], ['#ED1C24', '#B2009A'],
               ['#736900', '#80331A'], ['#008000', '#334D00']]          
    
    def __init__(self, ds:NMRDataSet, level_sign=0):
        self.root = tk.Tk() 
        self.attr = NMRDataSetAttributes(ds) 
        self.root.title(self.attr.get_ts_title())   
 
        self.data = [] # data store (cache), for faster redrawing
        self.add_data(ds)
        self.level_sign = level_sign # 0: positive & negative, 1: positive, -1: negative
              
        self.renderer = MplRenderer(figsize=self.figsize, show_frame=False)
        gs = gridspec.GridSpec(2, 2, figure=self.renderer._fig, width_ratios=(1,9), 
            height_ratios=(1,9), left=.03, right=.92, bottom=.08, top=.97, wspace=.01, hspace=.01)
                         
        self.ax_x = self.renderer._fig.add_subplot(gs[0, 1])
        self.ax_y = self.renderer._fig.add_subplot(gs[1, 0])
        self.ax_z = self.renderer._fig.add_subplot(gs[1, 1], sharex=self.ax_x, sharey=self.ax_y)
        
        xlabel = f'{ds.getPar("2s AXNUC")} (ppm)'
        ylabel = f'{ds.getPar("1s AXNUC")} (ppm)'
        self.ax_z.set(xlabel=xlabel, ylabel=ylabel)

        self.ax_z.yaxis.tick_right()
        self.ax_z.yaxis.set_label_position('right')
        plt.setp(self.ax_z.get_yticklabels(), rotation=90, va='center')
        
        self.ax_y.invert_xaxis()
        self.ax_y.axis('off')
        self.ax_x.axis('off')

        # retrieve xyz values from cache 
        x = self.data[0]['x']
        y = self.data[0]['y']
        z = self.data[0]['z']
        zmax = max(map(max, z))
            
        self.ax_z.set_xlim(x[0], x[-1])
        self.ax_z.set_ylim(y[0], y[-1])
        
        self.f2from = tk.DoubleVar()
        self.f1from = tk.DoubleVar()
        self.f2to = tk.DoubleVar()
        self.f1to = tk.DoubleVar()
        self.scale = tk.IntVar()

        self.f2from.set(f'{x[0]:0.2f}')
        self.f1from.set(f'{y[0]:0.2f}')

        self.f2to.set(f'{x[-1]:0.2f}')
        self.f1to.set(f'{y[-1]:0.2f}') 
        
        # create a canvas
        self.canvas = FigureCanvasTkAgg(self.renderer._fig, self.root) 
        self.canvas.draw()

        # save empty axes as backgrounds 
        self.bg_x = self.canvas.copy_from_bbox(self.ax_x.bbox) 
        self.bg_y = self.canvas.copy_from_bbox(self.ax_y.bbox) 
        self.bg_z = self.canvas.copy_from_bbox(self.ax_z.bbox)

        # plot contour lines and 1D projections
        self.contour(x, y, z, level_sign=self.level_sign)  
        self.ax_x.plot(x, np.max(z, axis=0)/zmax, 'b', lw=0.75)
        self.ax_y.plot(np.max(z, axis=1)/zmax, y, 'b', lw=0.75) 

        # create widgets
        canvas_widget = self.canvas.get_tk_widget()     
        self.ctrl1 = tk.Frame()
        self.ctrl2 = tk.Frame()        
        
        tk.Label(self.ctrl1, text='F2 [ppm]').grid(row=0, column=1, pady=(7,0), sticky='e')
        tk.Label(self.ctrl1, text='F2 [ppm]').grid(row=0, column=1, pady=(7,0), sticky='e')
        tk.Label(self.ctrl1, text='F1 [ppm]').grid(row=0, column=2, pady=(7,0), sticky='e')
        tk.Label(self.ctrl1, text='Intensity scale').grid(row=0, column=4, pady=(7,0), padx=(25,0), columnspan=3)
        tk.Label(self.ctrl1, text='External projection').grid(row=0, column=8, pady=(7,0), padx=(25,0), columnspan=2)
        tk.Label(self.ctrl1, text='Multiple display').grid(row=0, column=11, pady=(7,0), padx=(25,0))
        tk.Label(self.ctrl1, text='From').grid(row=1, column=0, padx=(25,0), sticky='w')
        tk.Label(self.ctrl1, text='To').grid(row=2, column=0, padx=(25,0), sticky='w')
        
        tk.Button(self.ctrl1, text='*2', command=self.scale_up).grid(row=1, column=4, padx=(25,0))
        tk.Button(self.ctrl1, text='/2', command=self.scale_down).grid(row=1, column=5)
        tk.Button(self.ctrl1, text='0', command=self.scale_not).grid(row=1, column=6, ipadx=2)    
        tk.Button(self.ctrl1, text=' F2 ', command=lambda: self.add_plot_dialog(projec='f2')).grid(row=1, column=8, padx=(25,0))
        tk.Button(self.ctrl1, text=' F1 ', command=lambda: self.add_plot_dialog(projec='f1')).grid(row=1, column=9)    
        tk.Button(self.ctrl1, text='Add spectrum', command=self.add_plot_dialog).grid(row=1, column=11, padx=(25,0))
        tk.Entry(self.ctrl1, textvariable=self.f2from, width=11).grid(row=1, column=1, sticky='w')
        tk.Entry(self.ctrl1, textvariable=self.f1from, width=11).grid(row=1, column=2, sticky='w')
        tk.Entry(self.ctrl1, textvariable=self.f2to, width=11).grid(row=2, column=1, sticky='w')
        tk.Entry(self.ctrl1, textvariable=self.f1to, width=11).grid(row=2, column=2, sticky='w')        
        tk.Label(self.ctrl2, text='Exact zoom: ').grid(row=0, column=0, padx=(25,0), pady=(4,7), sticky='w')
        tk.Button(self.ctrl2, text="Apply", command=self.zoom).grid(row=0, column=1, pady=(4,7), padx=(5,1))
        tk.Button(self.ctrl2, text="Reset", command=self.reset).grid(row=0, column=2, pady=(4,7), padx=1 )
              
        # (spacers)
        tk.Label(self.ctrl1, width=1).grid(row=0, column=3)
        tk.Label(self.ctrl1, width=1).grid(row=0, column=7)
        tk.Label(self.ctrl1, width=1).grid(row=0, column=10) 

        # pack the widgets
        NavigationToolbar2Tk(self.canvas, self.root)
        self.ctrl2.pack(side='bottom', anchor='w') 
        self.ctrl1.pack(side='bottom', anchor='w')
        canvas_widget.pack(fill='both', expand=True)
               
#        self.root.bind("<Configure>", self.on_resize)  
        self.root.resizable(0, 0)      
        self.root.mainloop()       
                      
    def add_plot_dialog(self, **kwargs):
        ''' Specify a new dataset for either 2D or 1D plotting'''
        self.top = tk.Toplevel(self.root)
        self.top.geometry("360x150")
        self.top.title('Add Spectrum')
        self.top.resizable(0, 0)
        
        fields = tk.Frame(self.top, pady=2)
        buttons = tk.Frame(self.top)

        tk.Label(fields, text="NAME = ").grid(column=0, row=0, sticky='w', padx=5, pady=1)
        tk.Label(fields, text="EXPNO = ").grid(column=0, row=1, sticky='w', padx=5, pady=1)
        tk.Label(fields, text="PROCNO = ").grid(column=0, row=2, sticky='w', padx=5, pady=1)
        tk.Label(fields, text="DIR = ").grid(column=0, row=3, sticky='w', padx=5, pady=1)
     
        self.name = tk.StringVar(value=self.attr.name)
        self.expno = tk.StringVar(value=self.attr.expno)
        self.procno = tk.StringVar(value=self.attr.procno)
        self.dir = tk.StringVar(value=self.attr.dir) 

        tk.Entry(fields, textvariable=self.name, width=29).grid(column=1, row=0, sticky='e', padx=5, pady=1)
        tk.Entry(fields, textvariable=self.expno, width=29).grid(column=1, row=1, sticky='e', padx=5, pady=1)
        tk.Entry(fields, textvariable=self.procno, width=29).grid(column=1, row=2, sticky='e', padx=5, pady=1)
        tk.Entry(fields, textvariable=self.dir, width=29).grid(column=1, row=3, sticky='e', padx=5, pady=1)
        
        tk.Button(buttons, text="Cancel", command=self.top.destroy).pack(side='right', padx=5, pady=5)
        tk.Button(buttons, text="OK", command=lambda: self.add_plot(projec=kwargs.get('projec'))).pack(side='right', pady=5)

        fields.pack(anchor='e') 
        buttons.pack(anchor='e')
        
        self.top.grab_set()
        self.top.wait_window()
        
    def add_data(self, ds:NMRDataSet):
        ''' Append xyz values and clevels from dataset to data store (cache) '''
        v = ds.getSpecDataMatrix()
        rng = v['physicalRanges']
        data = {}
        data['z'] = v['dataMatrix']
        data['x'] = np.linspace(rng[0]['start'], rng[0]['end'], num=data['z'].shape[1]) 
        data['y'] = np.linspace(rng[1]['start'], rng[1]['end'], num=data['z'].shape[0])  
        data['lvls'] = np.array(re.split("[|]", ds.getPar('LEVELSPAR LEVELS')), dtype='float')                      
        self.data.append(data)        
            
    def contour(self, x, y, z, levels=None, level_sign=0, level_color=0):
        ''' Generate contour lines and pass 'em to the renderer.'''
        print(level_sign)
        if levels is None: 
            levels = self.data[0]['lvls']
        if level_sign == 1:
            levels = levels[levels > 0]
        elif level_sign == -1:
            levels = levels[levels < 0]    
        cl = levels * 2**self.scale.get()

        colors = []
        for l in cl:
            if l > 0.: 
                colors.append(self.palette[level_color][0])
            else:
                colors.append(self.palette[level_color][1])

        # configure contour generator
        cont_gen = contour_generator(x=x, y=y, z=z, 
            name="serial", line_type=LineType.ChunkCombinedCode)
                
        for i, l in enumerate(cl):
            lines = cont_gen.lines(l)
            self.renderer.lines(lines, cont_gen.line_type, 
                ax=self.ax_z, color=colors[i], linewidth=0.7)
                                                           
    def scale_up(self):
        self.scale.set(self.scale.get() - 1)
        self.update_plot()

    def scale_down(self):
        self.scale.set(self.scale.get() + 1)
        self.update_plot()
        
    def scale_not(self):
        self.scale.set(0)
        self.update_plot()  
        
    def update_plot(self):         
        self.canvas.restore_region(self.bg_z)
        for coll in self.ax_z.collections: 
            coll.remove() 
        for i, data in enumerate(self.data):
            self.contour(data['x'], data['y'], data['z'], data['lvls'], level_sign=self.level_sign, level_color=i)            
        for coll in self.ax_z.collections: 
            self.ax_z.draw_artist(coll) 
        self.canvas.blit(self.ax_z.bbox)

    def zoom(self):
        self.ax_z.set_xlim(self.f2from.get(), self.f2to.get())
        self.ax_z.set_ylim(self.f1from.get(), self.f1to.get())  
        self.canvas.draw()

    def reset(self):
        x = self.data[0]['x']
        y = self.data[0]['y']
        z = self.data[0]['z']      
        self.f2from.set(f'{x[0]:0.4f}')
        self.f1from.set(f'{y[0]:0.4f}')        
        self.f2to.set(f'{x[-1]:0.4f}')
        self.f1to.set(f'{y[-1]:0.4f}')     
        self.zoom()
                
    def add_plot(self, projec=None):
        ''' Draw the data specified in add_data_dialog on the canvas'''
        ts = Topspin()
        path = self.dir.get()+'/'+self.name.get() + '/' + self.expno.get()+'/pdata/'+self.procno.get()    
        ds = ts.getDataProvider().getNMRData(path)       

        if ds == None: 
            ts.showError("No such file or directory ")
            self.top.destroy()

        if projec is None:
            if ds.getDimension() == 1: 
                top.showError("2D data set expected (1D given) ")
                self.top.destroy()                
                
            self.add_data(ds)
            data = self.data[-1]             
            self.contour(data['x'], data['y'], data['z'], data['lvls'], level_sign=self.level_sign, level_color=len(self.data)-1)      
            for i in range(len(data['lvls'])): 
                coll = self.ax_z.collections[-i]
                self.ax_z.draw_artist(coll) 
            self.canvas.blit(self.ax_z.bbox)
            
        else:
            if ds.getDimension() > 1: 
                ts.showError("1D data set expected (2D given) ")
                self.top.destroy() 
                                
            dv = ds.getSpecDataPoints()
            rng, = dv['physicalRanges']  
            z = dv['dataPoints']
            z = np.asarray(z)/max(z)
            x = np.linspace(rng['start'], rng['end'], len(z)) 
            
            if projec == 'f2': 
                self.canvas.restore_region(self.bg_x) 
                ln = self.ax_x.lines[0] 
                ln.set_data(x, z)
                self.ax_x.set_ylim(bottom=min(z)*2.)
                self.ax_x.draw_artist(ln)
                self.canvas.blit(self.ax_x.bbox)
            
            elif projec == 'f1':   
                self.canvas.restore_region(self.bg_y) 
                ln = self.ax_y.lines[0]
                ln.set_data(z, x)  
                self.ax_y.set_xlim(right=min(z))
                self.ax_y.draw_artist(ln)
                self.canvas.blit(self.ax_y.bbox)                
        self.top.destroy()    

class NMRDataSetAttributes: 

    def __init__(self, dataset:NMRDataSet):
        self._path = Path(dataset.getIdentifier())
        self.dir = str(self._path.parents[3])
        self.name = str(self._path.parts[-4])
        self.expno = str(self._path.parts[-3])
        self.procno = str(self._path.parts[-1])
            
    def is_processed(self):
        if list(self._path.glob('1r')) or list(self._path.glob('2rr')):
            return True
        else:
            return False
            
    def get_max_expno(self):
        p = self._path.parents[2]
        l = [x for x in p.iterdir() if x.is_dir()]
        l = [int(x.parts[-1]) for x in l if x.parts[-1].isdigit()]
        return max(l)   
        
    def get_ts_title(self):
        return f'{self.name}  {self.expno}  {self.procno}  {self.dir}'  
          
                
class PCA:
    ''' 
    Signal quantitation and approximation by means of principal component analysis
    
    Parameters: 
        X:  2D array of (complex) data   
        center:  whether the data to be centered before calling svd function
    '''
    def __init__( self, X, center=False):

        self.data = np.asarray(X)        
        if center:
            self.data = np.array(X)  
            self.mean = np.mean(X, axis=0) 
            self.data-= self.mean
            
        # SVD            
        self.U, self.d, self.Vt = np.linalg.svd( self.data, full_matrices=False )
        
        self.eigen = self.d**2		# PC eigenvalues
        self.pcs = self.Vt		# PC vectors (columns of V)
        self.scores = self.U * self.d	# PC scores 
        self.npc = len(self.d)		# the number of PCs
                      
    def evaluate( self, n=1 ):
        ''' evaluate signal integrals using first n principal components 
        '''
        assert 1 <= n <= self.npc
        data_approx = self.approximate(n=n)        
        return np.sum(data_approx, axis=1)
                           
    def evaluate_at_index( self, index=0, n=1 ):
        ''' evaluate signal at a given point using first n principal components 
        '''
        assert 1 <= n <= self.npc
        data_approx = self.approximate(n=n) 
        return data_approx[:, index]        
        
    def approximate( self, n=1):
        ''' reconstruct signal using first n principal components '''
        assert 1 <= n <= self.npc
        data_approx = self.scores[:,:n].dot(self.pcs[:n,:])
        if hasattr(self, 'mean'):
            data_approx += self.mean
        return data_approx

    def malinowski_indicator(self): 
        '''
        Malinowski's indicator of principal component significance.

        See G. Laurent et al. "Denoising applied to spectroscopies: parts I and II -
        programs and datas.", http://doi.org/10.5281/zenodo.1406172
     
        :return: number of significant principal components, for data approximation
        '''   
        r = max(self.data.shape)
        c = min(self.data.shape)

        eigenvalues = np.sqrt(self.eigen)
        real_error = np.sqrt(eigenvalues[::-1].cumsum()[::-1] /
            (r * (c - np.arange(c))))

        indicator = real_error / (c - np.arange(c)) **2
        n = np.argmin(indicator[:int(c / 2)])

        return max(1, int(n))
   
    def gavish_donoho_indicator(self, sigma=None): 
        '''
        Gavish and Donoho's indicator of principal component significance that guarantees 
        best asymptotic mean squared error between denoised and hypothetical noiseless data.

        See D. L. Donoho and M. Gavish, "The Optimal Hard Threshold for Singular
        Values is 4/sqrt(3)", http://arxiv.org/abs/1305.5870.
        
        :sigma: a scalar in X = X0 + sigma*Z, where X0 is noiseless data and Z the
                noise matrix; defaults to None, meaning 'an unknown noise level'
                
        :return: number of AMSE-optimal principal components, for data approximation
        '''  
        m, n = self.data.shape
        assert m <= n
        beta = m / n
        
        if sigma is not None:
            if isinstance(self.data[0,0], complex):
                sigma *= np.sqrt(2)
            w = (8* beta) / (beta + 1 + np.sqrt(beta**2 + 14 * beta + 1))
            lambda_star = np.sqrt(2* (beta + 1) + w)
            x = self.d > (lambda_star * np.sqrt(n) * sigma)
        else:
            omega = 0.56* beta**3 - 0.95* beta**2 + 1.82* beta + 1.43               
            x = self.d > (omega * np.median(self.d))
            
        return len(self.d[x]) if x.any() else 1


class ProgressbarThread(threading.Thread):
    def __init__(self, title='', maximum=None):
        threading.Thread.__init__(self, daemon=True)
        self.title = title 
        self.stop_event = threading.Event()
        
        if maximum:
            self.maximum = int(maximum)
            self.mode = 'determinate'
        else:    
            self.mode = 'indeterminate'   
             
        self.lock = threading.Lock()    # (1)
        self.lock.acquire()             # (2)
        self.start()

        # (1) Makes sure progressbar is fully loaded before executing anything
        with self.lock:
            return

    def close(self):
        self.stop_event.set()
        if self.mode == 'indeterminate':    
            os._exit(1)

    def run(self):
        self.root = tk.Tk()
        self.root.resizable(0,0)
        self.root.eval('tk::PlaceWindow . center')
        self.root.title(self.title)
        self.root.protocol("WM_DELETE_WINDOW", self.close)

        self.pb = ttk.Progressbar(self.root, length=400, mode=self.mode)
        if self.mode == 'determinate':  
            self.pb['value'] = 0
            self.pb['maximum'] = self.maximum
        else:
            self.pb.start(interval=20)

        self.pb.pack()
        self.lock.release()             # (2) Will release lock when finished
        self.root.mainloop()    

    def update(self, value: int):
        self.pb['value'] = value

class ProgressbarToplevel:
    def __init__(self, master, title='', maximum=None):
        self.master = master
        self.title = title
        
        if maximum:
            self.maximum = int(maximum)
            self.mode = 'determinate'
        else:    
            self.mode = 'indeterminate' 

        self.top = tk.Toplevel(self.master)
        self.top.title(self.title)

        x = self.master.winfo_x() + self.master.winfo_width()//2 - self.top.winfo_width()//2
        y = self.master.winfo_y() + self.master.winfo_height()//2 - self.top.winfo_height()//2
        self.top.geometry(f"400x20+{x}+{y}")
        
        self.top.resizable(0, 0)
        self.top.wait_visibility()
       
        self.pb = ttk.Progressbar(self.top, length=400, mode=self.mode)
        if self.mode == 'determinate':  
            self.pb['value'] = 0
            self.pb['maximum'] = self.maximum
        else:
            self.pb.start(interval=20)

        self.pb.pack()
        self.top.protocol("WM_DELETE_WINDOW", lambda: 0) # not able to close the pop

 
    def update(self, value: int):
        self.pb['value'] = value
        self.top.update()
        
    def close(self):
        self.top.destroy()
        

class ScrollEventHandler:

    def __init__(self, ax):
        self.ax = ax
        self.ymin = self.ax.get_ylim()[0]
        self.ymax = self.ax.get_ylim()[1]
        self.update()

    def on_scroll(self, event):
        scale = 1.25 if event.button == 'down' else 0.8
        self.ymin *= scale
        self.ymax *= scale
        self.update()

    def update(self):
        self.ax.set_ylim(self.ymin, self.ymax)
        self.ax.figure.canvas.draw_idle()    


def fft(data, fcor=0.5, axis=-1):
    ''' Fast Fourier transformation, an NMR style '''
    if axis < 0:
        axis += data.ndim
    if not 0 <= axis < data.ndim:
        raise IndexError("Not a valid axis for spectrum")
        
    slicing = (slice(None), ) * axis + (0, )
    data[slicing] *= fcor
    data_fft = np.fft.fftshift(np.fft.fft(data, axis=axis), axes=axis)
    data[slicing] /= fcor

    return data_fft 


