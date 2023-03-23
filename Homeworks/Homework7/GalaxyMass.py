import numpy as np
from ReadFile import Read
def ComponentMass(filename, particle_type):
    '''
    This function sorts through a data array of galactic components,
    and adds up the total mass of a specified type, either Halo,
    Disk, or Bulge.

    Inputs:
    - filename (string), the name of the data file to be parsed
    - particle_type (int), either 1, 2, or 3, depending on the
    particle type of interest.

    Outputs:
    - total mass (float) in units of 10^12 Msun, rounded to 3 decimal places
    '''
    data_array = Read(filename)[2] #constructing data array from Read function in ReadFile.py
    index = np.where(data_array['type'] == particle_type)  #array of indices of all specified particle types
    mass_list = data_array['m'][index] #array of the mass values  at each type index, each val in 10^10 Msun
    sum_mass = np.around(np.sum(mass_list)/100, 3) #total mass in 10^12 Msun (hence /100) and rounded to 3 decimal places
    return sum_mass

'''
below are 3 print statements I use to extract different masses and baryonic ratios

print(ComponentMass('M33_000.txt', 1),  ComponentMass('M33_000.txt', 2),  ComponentMass('M33_000.txt', 3))
print(ComponentMass('M33_000.txt', 1) + ComponentMass('M33_000.txt', 2) + ComponentMass('M33_000.txt', 3))
print(np.around(((ComponentMass('M33_000.txt', 2) + ComponentMass('M33_000.txt', 3))/(ComponentMass('M33_000.txt', 1) + ComponentMass('M33_000.txt', 2) + ComponentMass('M33_000.txt', 3))), 3))
'''
