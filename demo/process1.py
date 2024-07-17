# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 17:33:25 2024

@author: dell
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from parameter_list import *
from function_library import *
import multiprocessing
import time


phiP, thetaP = np.meshgrid(phiP_list, thetaP_list)
x = R2 * np.cos(phiP) * np.cos(thetaP)
y = R2 * np.cos(phiP) * np.sin(thetaP)
z = R2 * np.sin(phiP)

# Plotting the sphere
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot the surface with intensity as color
mappable = ax.plot_surface(x, y, z, facecolors=plt.cm.gray(Intensity.T /np.max(Intensity)), rstride=1, cstride=1, antialiased=False)

# Plot the incident and reflected vectors
ax.quiver(-(2000 + R2) * np.cos(Theta), -(2000 + R2) * np.sin(Theta), 0, np.cos(Theta), np.sin(Theta), 0, color='r', length=2000.0, normalize=True)
ax.quiver(R2 * camera[0], R2 * camera[1], R2 * camera[2], camera[0], camera[1], camera[2], color='g', length=2000.0, normalize=True, linestyle='dashed')

# Set axis labels
ax.set_xlabel('X (km)') 
ax.set_ylabel('Y (km)') 
ax.set_zlabel('Z (km)') 

# Set the view angle 
elev = np.arcsin(camera[2])
if camera[2] == 1:
    azim = 0
else:
    azim = np.arccos(camera[0]/np.cos(elev))
ax.view_init(elev=np.rad2deg(elev) , azim=np.rad2deg(azim))

# Show the plot
plt.show()
