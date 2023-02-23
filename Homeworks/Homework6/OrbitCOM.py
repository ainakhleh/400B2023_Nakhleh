

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from modded_COM import CenterOfMass




def OrbitCOM(galaxy, start, end, n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
    galaxy (string): the name of the galaxy of interest, i.e. MW, M31, M33
    start (int): number of the first snapshot in the simulation to be read in
    end (int): number of the last snapshot in the simulation to be read in
    n (int): the interval size over which COM will be returned
          
    outputs: a text file with 7 columns of information -> time, x position, y position, 
    z position, x velocity, y velocity, z velocity
    """
    
    # compose the filename for output
    fileout = "Orbit" + galaxy + "_" + str(start) + "_" + str(end) + ".txt"
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1 #the tolerance level for our COM calculation method
    volDec = 2 #stands for volume devreased. used for MW and M31 for the shrinking spheres method
    volDec_M33 = 4 #for M33 because it becomes tidally stripped, so the COM needs to be found more quickly
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start, end + 1, n)
    if len(snap_ids) == 0:
        print('provide a non-empty time list') #adding a checker to make sure the array list isn't zero
        exit()
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids), 7]) #7 values, so our array needs to have 7 columns for each value

    
    # a for loop 
    for  i, snap_id in enumerate(snap_ids):# loop over all files
        
        # compose the data filename (be careful about the folder)
        ilbl = '000' + str(snap_id)
        ilbl = ilbl[-3:]
        filename = "timestamps/" + galaxy + '_' + str(ilbl) + '.txt' #my folder holding all snapshots is called timestamps
        # Initialize an instance of CenterOfMass class, using disk particles
        # NOTE, to save space on my computer I had to delete the snapshot folders, so I will need to download timestamps to reinvestigate

        # Store the COM pos and vel. Remember that now COM_P required VolDec
        COM = CenterOfMass(filename, 2) #looking at disk particles right now
        if galaxy == 'M33':
            COM_pos = COM.COM_P(delta, volDec_M33) #exception for shrinking spheres method of M33 due to tidal stripping
            COM_vel = COM.COM_V(COM_pos[0], COM_pos[1], COM_pos[2])
        else:
            COM_pos = COM.COM_P(delta, volDec)
            COM_vel = COM.COM_V(COM_pos[0], COM_pos[1], COM_pos[2])

        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        time = COM.time / 1000
        orbit[i][0] = time.value
        orbit[i][1] = COM_pos[0].value
        orbit[i][2] = COM_pos[1].value
        orbit[i][3] = COM_pos[2].value
        orbit[i][4] = COM_vel[0].value
        orbit[i][5] = COM_vel[1].value
        orbit[i][6] = COM_vel[2].value
        # print snap_id to see the progress
        #print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))





# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 

#OrbitCOM("M33", 0, 800, 5) 
#OrbitCOM("M31", 0, 800, 5)
#OrbitCOM("MW", 0, 800, 5)

#above three files have already been made, so no need to call it anymore






# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt

MW_orbit_data = np.genfromtxt("OrbitMW_0_800.txt",dtype=None,names=True,skip_header=0) #reading in COM data evolution for MW
M31_orbit_data = np.genfromtxt("OrbitM31_0_800.txt",dtype=None,names=True,skip_header=0) #reading in COM data evolution for M31
M33_orbit_data = np.genfromtxt("OrbitM33_0_800.txt",dtype=None,names=True,skip_header=0) #reading in COM data evolution for M33



# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  

def separation(x1,x2,y1,y2,z1,z2):
    return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) #typical distance formula for 3d space
    #could make this n-dimensional but not needed for this problem




# Determine the magnitude of the relative position and velocities 

# of MW and M31

# of M33 and M31

time_list = MW_orbit_data['t'] #same time for all three galaxies, so this can be a shared x-axis for all plots
MW_M31_separations = np.array([]) #will be filled in with separation distances between MW and M31
M33_M31_separations = np.array([]) #will be filled in with separation distances between M33 and M31
MW_M31_vel_separations = np.array([]) #will be filled in with relative velocity between MW and M31
M33_M31_vel_separations = np.array([]) #will be filled in with relative velocity between M33 and M31

for i in range(0, len(MW_orbit_data['x'])):
    MW_M31_separations = np.append(MW_M31_separations, separation(MW_orbit_data['x'][i], M31_orbit_data['x'][i], 
    MW_orbit_data['y'][i], M31_orbit_data['y'][i], MW_orbit_data['z'][i], M31_orbit_data['z'][i]))

    # MW_orbit_data['x'][i] represents the x value at every time incriment i for MW
    # same pattern for the rest of the variables to be fed into separation function
    # repeat below for all other data arrays of interest
    # all arrays have the same length so they can all be filled in within the same for loop

    M33_M31_separations = np.append(M33_M31_separations, separation(M33_orbit_data['x'][i], M31_orbit_data['x'][i], 
    M33_orbit_data['y'][i], M31_orbit_data['y'][i], M33_orbit_data['z'][i], M31_orbit_data['z'][i]))
    

    MW_M31_vel_separations = np.append(MW_M31_vel_separations, separation(MW_orbit_data['vx'][i], M31_orbit_data['vx'][i], 
    MW_orbit_data['vy'][i], M31_orbit_data['vy'][i], MW_orbit_data['vz'][i], M31_orbit_data['vz'][i]))

    M33_M31_vel_separations = np.append(M33_M31_vel_separations, separation(M33_orbit_data['vx'][i], M31_orbit_data['vx'][i], 
    M33_orbit_data['vy'][i], M31_orbit_data['vy'][i], M33_orbit_data['vz'][i], M31_orbit_data['vz'][i]))




# Plot the Orbit of the galaxies 
#################################

fig, (ax1, ax2) = plt.subplots(2) #creating a subplot to include separations between MW, M31, and M31, M33
fig.suptitle('Galactic Separations throughout 12 Gyr')
ax1.plot(time_list, MW_M31_separations, label = 'Separation between MW and M31')
ax1.set_ylabel('Separation (kpc)')
ax1.legend()
ax2.plot(time_list, M33_M31_separations, label = 'Separation between M33 and M31')
ax2.set_ylabel('Separation (kpc)')
ax1.set_yscale('log') #logged to get better visualization at far times when galaxies have merged
ax2.legend()
ax2.set_xlabel('time (Gyr)') #both plots can share the same x-axis (time)
plt.show()


# Plot the orbital velocities of the galaxies 
#################################

fig, (ax1, ax2) = plt.subplots(2) #creating another subplot to include veolcities between MW, M31, and M31, M33
fig.suptitle('Galactic Relative Velocities throughout 12 Gyr')
ax1.plot(time_list, MW_M31_vel_separations, label = 'Relative Velocity between MW and M31')
ax1.set_ylabel('Relative Velocity (km/s)')
ax1.legend()
ax2.plot(time_list, M33_M31_vel_separations, label = 'Relative Velocity between M33 and M31')
ax2.set_ylabel('Relative Velocity (km/s)')
ax2.legend()
ax2.set_xlabel('time (Gyr)') #both plots can share the same x-axis (time)
plt.show()