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
    r = pa[i][:, 1]

    # Create a new figure with a polar subplot
    plt.figure()
    ax = plt.subplot(projection='polar')
    ax.set_ylim([-50, 40])

    # Apply the color mapping to 'r' values and create colored bars
    cbar = plt.cm.viridis(rescale_color(r))
    ax.bar(theta, r, width=5 * np.pi / 180, color=cbar)

    # You can customize labels, titles, and other plot settings here if needed

# Show the plots
plt.show()
