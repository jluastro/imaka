import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from astropy.io import fits

def first_light_images():
    """
    Plot up the first light images for the press release.
    """
    work_dir = '/Users/jlu/work/imaka/pleiades/press_release/'

    img_op = fits.getdata(work_dir + 'west_stack_open.fits')
    img_tf = fits.getdata(work_dir + 'west_stack_ttf.fits')
    img_cl = fits.getdata(work_dir + 'west_stack_closed.fits')


    # norm = plt.Normalize(img_cl.min(), img.max())
    norm = colors.LogNorm(vmin=30, vmax=2000)

    img_height = float(img_cl.shape[0])
    img_width = float(img_cl.shape[1])

    # Make the large scale closed-loop image.
    fig = plt.figure(1)
    fig.set_size_inches(3*img_width/img_height, 3, forward=False)
    plt.clf()
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    ax = plt.axes([0, 0, 1, 1])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(img_cl+50, cmap=cm.Greys_r, norm=norm)
    plt.savefig(work_dir + 'full_closed.png', pad_inches=0.05, dpi=300)


    # Make the large scale open-loop image.
    fig = plt.figure(1)
    fig.set_size_inches(3*img_width/img_height, 3, forward=False)
    plt.clf()
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    ax = plt.axes([0, 0, 1, 1])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(img_op+50, cmap=cm.Greys_r, norm=norm)
    plt.savefig(work_dir + 'full_open.png', pad_inches=0.05, dpi=300)

    
    
    

    
