#!/usr/bin/env python

"""
Visualize crowd simulation results

NB: Requires a modern Matplotlib version; also needs
 either FFMPeg (for MP4) or ImageMagick (for GIF)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as manimation
import sys




# Required Metadata: Number of Particles, Bounding Box

def main(infile="particles.csv", outfile="out.mp4"):


    # Parse CSV File
    df = pd.read_csv(infile)
    PDirX = df['PDirX']
    PDirY = df['PDirY']
    PTag = df['PTag']
    PId = df['PId']
    PLocX = df['PLocX']
    PLocY = df['PLocY']

    num_particles = len(set(PId))
    nframe = len(PLocX)//num_particles
    print("nframe = {0}".format(nframe))
    
    # Define map from Particle ID to Colors (assume only two unique particle tags) (Optional)
    particle_tags = set(PId)
    particle_tags = list(particle_tags)
    dict_colors = {}
    for i in range(num_particles):
        pid = PId[i]
        if(PTag[i] == particle_tags[0]):
            dict_colors[pid] = 'r';
        if(PTag[i] == particle_tags[1]):
            dict_colors[pid] = 'b'

    def plot_frame(i):
        ax = fig.add_subplot(111)
        ax.set_xlim([0,1])
        ax.set_ylim([0,1])
        xp = PLocX[i*num_particles:(i+1)*num_particles]
        yp = PLocY[i*num_particles:(i+1)*num_particles]
        xd = PDirX[i*num_particles:(i+1)*num_particles]
        yd = PDirY[i*num_particles:(i+1)*num_particles]
        pids = PId[i*num_particles:(i+1)*num_particles]
        ax.quiver(xp,yp,xd,yd)
        ax.scatter(xp, yp, s = 200, c = PTag[i*num_particles:(i+1)*num_particles] + 1)
        # , c=map(dict_colors.get, pids)
        return ax
    
    fig = plt.figure(figsize=(20,20))



    metadata = dict(title='Crowd Simulation', artist='Matplotlib')
    if outfile[-4:] == ".mp4":
        Writer = manimation.writers['ffmpeg']
        writer = Writer(fps=15, metadata=metadata,
                        extra_args=["-r", "30",
                                    "-c:v", "libx264",
                                    "-pix_fmt", "yuv420p"])
    elif outfile[-4:] == ".gif":
        Writer = manimation.writers['imagemagick']
        writer = Writer(fps=15, metadata=metadata)

    with writer.saving(fig, outfile, nframe):
        for i in range(nframe):
            ax = plot_frame(i)
            writer.grab_frame()
            plt.delaxes(ax)


if __name__ == "__main__":
    main(*sys.argv[1:])
