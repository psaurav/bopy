
import sys
import os
import re
from matplotlib.widgets import Cursor, Button, Slider
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import bopy as bp


class Viewer2d_xsection(object):
    def __init__(self,z,x=None, y=None):
        """
        Shows a given array in a 2d-viewer.
        Input: z, an 2d array.
        x,y coordinters are optional.
        """
        from matplotlib import cbook
        plt.rcParams['font.size']=10
        if x==None:
            #self.x=np.arange(z.shape[0])
            self.x=np.arange(z.shape[1])
        else:
            self.x=x
        if y==None:
            #self.y=np.arange(z.shape[1])
            self.y=np.arange(z.shape[0])
        else:
            self.y=y
        self.z=z
        zmin = np.min(z)
        zmax = np.max(z)
        rcmax = np.max(z.shape) - 1
        print >>sys.stderr, "len(x):", len(self.x)
        print >>sys.stderr, "len(y):", len(self.y)
        self.fig=plt.figure()
        self.fig.canvas.set_window_title('gim2ds')
        #Doing some layout with subplots:
        self.overview=plt.subplot2grid((8,4),(0,0),rowspan=6,colspan=2)
        self.cimg = self.overview.imshow(self.z)
        self.overview.set_xlabel('column')
        self.overview.set_ylabel('row')
        self.overview.set_title('Value')
        self.cbar = self.fig.colorbar(self.cimg, use_gridspec=True, orientation='horizontal')
        self.cbar.locator = MaxNLocator(nbins = 5)
        self.x_subplot=plt.subplot2grid((8,4),(0,2),rowspan=4,colspan=2)
        self.x_subplot.set_xlabel('column')
        self.x_subplot.set_ylabel('Value')
        self.x_subplot.set_xlim([0, rcmax])
        #self.x_subplot.(MaxNLocator(5))
        self.y_subplot=plt.subplot2grid((8,4),(4,2),rowspan=4,colspan=2)
        self.y_subplot.set_xlabel('row')
        self.y_subplot.set_ylabel('Value')
        self.y_subplot.set_xlim([0, rcmax])

        #Adding widgets, to not be gc'ed, they are put in a list:
        cursor=Cursor(self.overview, useblit=True, color='black', linewidth=2 )
        self.max_ax=plt.subplot2grid((8,4),(6,0),colspan=1)
        self.min_ax=plt.subplot2grid((8,4),(7,0),colspan=1)
        self.smin = Slider(self.min_ax, 'Min', zmin, zmax, valinit=zmin)
        self.smax = Slider(self.max_ax, 'Max', zmin, zmax, valinit=zmax)
        self.smin.on_changed(self.update)
        self.smax.on_changed(self.update)
        self._widgets=[cursor, self.cbar]
        self.tlo = plt.tight_layout()
        #connect events
        self.fig.canvas.mpl_connect('button_press_event',self.click)
        self.fig.canvas.mpl_connect('key_press_event', self.press)
        plt.show()

    def show_legend(self, event):
        """Shows legend for the plots"""
        for pl in [self.x_subplot,self.y_subplot]:
            if len(pl.lines)>0:
                pl.legend()
        plt.draw()

    def clear_xy_subplots(self, event):
        """Clears the subplots."""
        for j in [self.overview,self.x_subplot,self.y_subplot]:
            j.lines=[]
            j.legend_ = None
        plt.draw()

    def press(self, event):
        print 'pressed:', event.key
        if event.key=='r':
            self.clear_xy_subplots(event)
        if event.key=='v':
            s = raw_input('Value text:')
            self.overview.set_title(s)
            self.x_subplot.set_ylabel(s)
            self.y_subplot.set_ylabel(s)
            plt.draw()
        #    s = raw_input("Title: ")
        #    ax.set_title(s)
        #    fig.canvas.draw()


    def click(self,event):
        """
        What to do, if a click on the figure happens:
            1. Check which axis
            2. Get data coord's.
            3. Plot resulting data.
            4. Update Figure
        """
        if event.inaxes==self.overview:
            #Get nearest data
            #xpos=np.argmin(np.abs(event.xdata-self.x))
            #ypos=np.argmin(np.abs(event.ydata-self.y))
            xpos=int(event.xdata)
            ypos=int(event.ydata)
            if event.button==3:
                #Plot it                
                c,=self.y_subplot.plot(self.y, self.z[:,xpos],label=str(self.x[xpos]))
                self.overview.axvline(self.x[xpos],color=c.get_color(),lw=2)
                c,=self.x_subplot.plot(self.x, self.z[ypos,:],label=str(self.y[ypos]))
                self.overview.axhline(self.y[ypos],color=c.get_color(),lw=2)

        if event.inaxes==self.y_subplot:
            ypos=np.argmin(np.abs(event.xdata-self.y))
            c=self.x_subplot.plot(self.x, self.z[ypos,:],label=str(self.y[ypos]))
            self.overview.axhline(self.y[ypos],color=c.get_color(),lw=2)

        if event.inaxes==self.x_subplot:
            xpos=np.argmin(np.abs(event.xdata-self.x))
            c,=self.y_subplot.plot(self.y, self.z[:,xpos],label=str(self.x[xpos]))
            self.overview.axvline(self.x[xpos],color=c.get_color(),lw=2)
        #Show it

        self.show_legend(event)

        plt.draw()

    def update(self, val):
        self.cimg.set_clim([self.smin.val,self.smax.val])
        self.fig.canvas.draw()

if __name__=='__main__':
    #Build some strange looking data:
    try:
        filename = sys.argv[1]
    except:
        print >>sys.stderr, "Usage:", sys.argv[0], "<filename>"
        sys.exit(1)

    d = bp.io.read_frame(filename)
    fig_v=Viewer2d_xsection(d)
    #Show it
    #plt.show()
