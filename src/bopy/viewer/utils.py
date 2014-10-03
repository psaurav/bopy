import os
import numpy as np
import matplotlib.pyplot as plt

cdict3 = {'red':  ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.8, 1.0),
                   (0.75,1.0, 1.0),
                   (1.0, 0.4, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.9, 0.9),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.4),
                   (0.25,1.0, 1.0),
                   (0.5, 1.0, 0.8),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }

def showbluered(data):
    plt.register_cmap(name='BlueRed3', data=cdict3)
    plt.clf()
    fig = plt.figure()
    ax = plt.subplot(111)
    mx = np.abs(np.max(data))
    mn = np.abs(np.min(data)) 
    vmax = np.max((mx, mn))
    ax.imshow(data, interpolation='nearest', vmax=vmax, vmin=-vmax)
    ax.set_cmap('BlueRed3')
    plt.show()

def myimshow(data, interpolation='nearest', filename='',
             vmin='', vmax='', title='', xlabel='', ylabel='', cbarunits=''):
    #import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    plt.clf()
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    divider = make_axes_locatable(ax)
    
    if vmin is not '' and vmax is not '':
        im = ax.imshow(data, interpolation=interpolation, vmin=vmin, vmax=vmax)
    else:
        im = ax.imshow(data, interpolation=interpolation)

    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)

    if cbarunits is not '':
        cbar.set_label(cbarunits)

    if filename is not '':
        plt.savefig(filename, bbox_inches='tight')
    else:

        def onkey(event):
            if event.key == 't':
                tstr = raw_input('Input title: ')
                print tstr
                ax.set_title(tstr)
                plt.draw()
            if event.key == 'x':
                xstr = raw_input('Input xlabel: ')
                ax.set_xlabel(xstr)
                plt.draw()
            if event.key == 'y':
                ystr = raw_input('Input xlabel: ')
                ax.set_xlabel(ystr)
                plt.draw()
        kid = fig.canvas.mpl_connect('key_press_event', onkey)

        numrows, numcols = data.shape
        def format_coord(x, y):
            col = int(x+0.5)
            row = int(y+0.5)
            if col>=0 and col<numcols and row>=0 and row<numrows:
                z = data[row,col]
                return '(%4d, %4d) z=%1.3e'%(row, col, z)
            else:
                return 'x=%1.3f, y=%1.3f'%(x, y)
        ax.format_coord = format_coord
        plt.show()


def show_ppdata(ddir, rootext, wl, kind='phi'):
    import os
    import numpy as np
    ddir = ddir.rstrip('/')


    #from mpl_toolkits.axes_grid1 import AxesGrid
    import matplotlib.gridspec as gridspec
    fig = plt.figure()

    nrow = 11
    ncol = 19
    gs = gridspec.GridSpec(nrow, ncol+1)
    for r in range(nrow):
        for c in range(ncol):
            i = c + ncol*r + 1
            fpath = '{0}_wl{1}_s{2}_{3}.npy'.format(rootext, wl, i, kind)
            fpath = os.path.join(ddir, fpath)
            d = np.load(fpath)
            ax = plt.subplot(gs[r, c])
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_title('%d'%(i), size='x-small')
            im = ax.imshow(d, interpolation='nearest')
    
    plt.tight_layout()
    plt.show()

def transill2(ddir, rootext, ddirr, rootextr, wl, kind='A', collapse='mean', mask=None, srcs=None, filename=''):
    from bopy.gen3pp.utils import get_transill2
    d = get_transill2(ddir, rootext, ddirr, rootextr, wl, kind=kind, collapse=collapse, mask=mask, srcs=srcs)
    if filename is not '':
        npyfname = os.path.join('{0}.npy'.format(os.path.splitext(filename)[0]))
        np.save(npyfname, d)
    myimshow(d, title='{0} ({1}nm) ({2} - {3})'.format(kind, wl, rootext, rootextr), filename=filename)

