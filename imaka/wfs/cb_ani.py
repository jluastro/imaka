#Animated WFS

import matplotlib
matplotlib.use('Agg')
from glob import glob
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.animation as animation

def animate_cb(cb_dir, out_dir):
    date = cb_dir.split('/')[-3]
    files = glob(cb_dir+"aocb*")
    frames = []
    times = []
    names = []
    for file in files:

        name = file.split('cb')[-1].split('.')[0]

        hdul = fits.open(file)
        data_raw = hdul[1].data
        hdr = hdul[1].header
        hdul.close()
        stack = np.median(data_raw, axis=0)

        WFS0 = stack[0]
        WFS1 = stack[1]
        WFS2 = stack[2]
        WFS3 = stack[3]

        t0 = hdr['TSTAMPA0']
        t1 = hdr['TSTAMPA1']
        t2 = hdr['TSTAMPA2']
        t3 = hdr['TSTAMPA3']
        t = [t0, t1, t2, t3]

        top = np.hstack([WFS0, WFS1])
        bottom = np.hstack([WFS2, WFS3])
        img = np.vstack([top, bottom])
        frames.append(img)
        times.append(t)
        names.append(name)

    fig = plt.figure(figsize=(8,8))

    vmin = np.min(frames[0])
    vmax = np.max(frames[0])

    im = plt.imshow(frames[0], vmin=vmin, vmax=vmax, animated=True)
    file_text = plt.text(85, 96, '', color='r', fontsize=16)

    WFS0_text = plt.text(3, 15, '', color='r', fontsize=12)
    WFS1_text = plt.text(98, 15, '', color='r', fontsize=12)
    WFS2_text = plt.text(3, 111, '', color='r', fontsize=12)
    WFS3_text = plt.text(98, 111, '', color='r', fontsize=12)

    plt.title(date)
    plt.axvline(96, color='r')
    plt.axhline(96, color='r')
    plt.text(3, 5, 'WFS0', color='r')
    plt.text(98, 5, 'WFS1', color='r')
    plt.text(3, 101, 'WFS2', color='r')
    plt.text(98, 101, 'WFS3', color='r')
    plt.colorbar()

    def updatefig(frame):
        i=int(frame)
        im.set_array(frames[i])
        file_text.set_text(names[i])
        WFS0_text.set_text(times[i][0])
        WFS1_text.set_text(times[i][1])
        WFS2_text.set_text(times[i][2])
        WFS3_text.set_text(times[i][3])

        return im, file_text, WFS0_text, WFS1_text, WFS2_text, WFS3_text,

    ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True)
    title = out_dir+"cb_ani_"+date+".mp4"
    ani.save(title)
    
    return
