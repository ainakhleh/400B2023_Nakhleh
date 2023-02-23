import numpy as np
import astropy.units as u
def Read(filename):
    '''
    opens a datafile and finds information about the time, number of particles,
    then constructs the data entries as a NumPy data array

    Inputs: filename as a string

    Returns: timestamp associated with file, number of particles in the file,
    and the data array
    '''
    info_file = open(filename, "r") #opening the file to extract data
    line1 = info_file.readline() #first line of file
    label, value = line1.split() #splitting "Time" string and numerical time into two variables
    time = float(value) * u.Myr #assigns time Myr units via AstroPy
    line2 = info_file.readline() #2nd line of file
    label, value = line2.split() #splitting "Total" string and numerical particle number into two variables
    num_particles = float(value) #no units
    info_file.close()
    data = np.genfromtxt(filename, dtype = None, names = True, skip_header = 3)
    #turns space seperated data entries into condensed NumPy array

    return time, num_particles, data

