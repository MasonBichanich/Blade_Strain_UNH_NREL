# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 14:41:34 2023

@author: Mason Bichanich
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as spio

# Close all previous plots
plt.close('all')

# Load data from .mat files
mat = spio.loadmat('pa_shifted.mat', squeeze_me=True, simplify_cells=True)
pa_shifted = mat['pa_shifted']

mat2 = spio.loadmat('phase_averaged.mat', squeeze_me=True, simplify_cells=True)
pa = mat2['phase_averaged']

# Define a color mapping function to rescale 'r' values
def rescale_color(r):
    return (r - np.min(r)) / (np.max(r) - np.min(r))

# Create polar bar plots for each dataset
for i in np.arange(0, 8):
    theta = np.deg2rad(pa[i][:, 0])  # Convert degrees to radians
    r = pa_shifted[i]

    # Create a new figure with a polar subplot
    plt.figure()
    ax = plt.subplot(projection='polar')
    ax.set_ylim([-50, 50])
    ax.set_yticks([-40,-30,-20,-10,0,10,20,30,40])
    ax.set_rlabel_position(138)
    ax.set_theta_zero_location('W')
    ax.set_theta_direction(-1)
    ax.text(2.2, 2, 'Microstrain',
        rotation=45, ha='center', va='center')
    ax.set_xlabel('Degree of rotation')

    # Apply the color mapping to 'r' values and create colored bars
    cbar = plt.cm.viridis(rescale_color(r))
    ax.bar(theta, r, width=5 * np.pi / 180, color=cbar)

    # You can customize labels, titles, and other plot settings here if needed

# Show the plots
plt.show()
mins = np.arange(0,8,dtype = 'float')
maxs = np.arange(0,8,dtype = 'float')
means = np.arange(0,8,dtype = 'float')
stds = np.arange(0,8,dtype = 'float')
for i in range(0,8):
    stds[i] = np.std(pa_shifted[i])
    mins[i] = min(pa_shifted[i])
    maxs[i] = max(pa_shifted[i])
    means[i] = np.mean(pa_shifted[i])
