import numpy as np
import matplotlib.pyplot as plt

def im_show2d(d, interpolation='nearest'):
    """Customize imshow

    Args:
        d: numpy.array
        interpolation: imshow option (default=nearest)
    
    """
    fig = plt.figure()
    ax = plt.gca()
    im = plt.imshow(d, interpolation=interpolation)

    #
    # format output string
    #
    def report_pixel(x, y):
        # get the pixel value 
        x = int(x+0.5)
        y = int(y+0.5)
        v = d[y,x]
        return "%f    r=%d, c=%d" % (v, y, x)
    ax.format_coord = report_pixel

    #
    # align colorbar
    #
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%")
    cb = plt.colorbar(im, cax=cax)

    #
    # show 
    #
    plt.tight_layout()
    plt.show()
