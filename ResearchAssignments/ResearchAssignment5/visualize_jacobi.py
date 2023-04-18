# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile
from GalaxyMass import ComponentMass


def scatterplot_galaxy(filename):
    time, total, data = Read(filename)
    index = np.where(data['type'] == 2)
    xs = data['x'][index]
    ys = data['y'][index]
    zs = data['z'][index]
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(xs,ys,zs)
    ax.set_xlabel('x (kpc)')
    ax.set_ylabel('y (kpc)')
    ax.set_zlabel('z (kpc)')
scatterplot_galaxy('M31_270.txt')
plt.show()

