import numpy as np
import astropy.units as u
import astropy.table as tbl
from ReadFile import Read
from CenterofMass import CenterOfMass
from astropy.constants import G
import matplotlib.pyplot as plt
from GalaxyMass import ComponentMass

class MassProfile:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot

    def __init__(self, galaxy, snap):
        ''' Class to calculate the Mass and velocity profiles of galaxies MW, M31, and M33
        at any snapshot
            
            PARAMETERS
            ----------
            galaxy : `str`
                galaxy name (MW, M31, M33)
            ptype : `snap`
                timestamp of galaxies in the calculation
        '''
        ilbl = '000' + str(snap)
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy) + ilbl + '.txt'
        #this strategy is used to extract the file name from the galaxy name and snapshot number (time 0, 1, 100, etc)
        self.time, self.total, self.data = Read(self.filename)
        self.x = self.data['x'] * u.kpc #adding proper units
        self.y = self.data['y'] * u.kpc #adding proper units
        self.z = self.data['z'] * u.kpc #adding proper units
        self.m = self.data['m']
        self.gname = galaxy 
    
    def MassEnclosed(self, ptype, radii):
        '''
        This function returns the enclosed mass of a specified particle type within a 
        certain distance from the COM of the galaxy of interest, with proper units of Msun

        inputs:
            ptype (int): the type of particle we want to find the enclosed mass of, either
            1, 2, or 3

            radii (astropy quantity): an array of radii with units of kpc, provided to find the total enclosed mass
            within several distances of the COM in series for the Mass and velocity profile.

        outputs:
            mass_array: an array of masses, representing the total enclosed mass of the specified 
            particle type within the radius at the matching index of the radii input array.
        '''
        COM = CenterOfMass(self.filename, ptype) #calling center of mass function
        COM_P_array = COM.COM_P(0.1) #finding center of mass within the CenterOfMass class

        x_COM = COM_P_array[0] #splitting array into individual x coord
        y_COM = COM_P_array[1] #splitting array into individual y coord
        z_COM = COM_P_array[2] #splitting array into individual z coord

        type_index = np.where(self.data['type'] == ptype)[0] #we only want to look at coords of one type of particle at a time
        x_new = self.x[type_index] - x_COM  #finding relative x distance from COM x coordinate
        y_new = self.y[type_index] - y_COM  #finding relative y distance from COM y coordinate
        z_new = self.z[type_index] - z_COM #finding relative z distance from COM z coordinate
        m_new = self.m[type_index]  #only want masses inside the radius of the specified particle type
        mass_array = np.array([]) #empty array that will be filled will mass sums

        distance = np.sqrt((x_new)**2 + (y_new)**2 + (z_new)**2) 
        #applying distance formula to find if particle is within a given radius of the COM
        
        for i in range(0, len(radii)): #looping through radius values of interest
            index_within = np.where(distance < radii[i]) #only selecting indices of particles within specified radius of COM
            masses = m_new[index_within] #finding the masses of specified particle type associated with these indices
            tot_mass = np.sum(masses) #summing all masses within the specified radius
            mass_array = np.append(mass_array, tot_mass * 1e10) #adding summed mass to the array for the specified radius
        return mass_array * u.Msun #adding appropriate units of Msun

    def MassEnclosedTotal(self, radii):
        '''
        This function returns the total mass of a all particle types within a 
        certain distance from the COM of the galaxy of interest, with proper units of Msun.

        inputs:
            radii (astropy quantity): an array of radii with units of kpc, provided to find the total enclosed mass
            within several distances of the COM in series for the Mass and velocity profile.

        outputs:
            tot_mass_array: an array of masses, representing the total enclosed mass of all specified 
            particle types within the radius at the matching index of the radii input array.
        '''
        halo_array = self.MassEnclosed(1, radii) #type 1 matter is halo matter
        disk_array = self.MassEnclosed(2, radii) #type 2 matter is disk matter
        if self.gname == 'M33':
            tot_mass_array = halo_array + disk_array 
            # M33 has no bulge matter, so we only add halo matter and disk matter together (at every radius incriment), 
            # in order to get total enclosed mass within a radius away from the COM
            
            return tot_mass_array
        else:
            bulge_array = self.MassEnclosed(3, radii) #type 3 matter is bulge matter
            tot_mass_array = halo_array + disk_array + bulge_array #MW and M31 have bulge matter
            return tot_mass_array
    def HernquistMass(self, radius, a, Mhalo):
        '''
        This function returns the Hernquist mass at a specified radius, depending on the
        galaxy's scale factor, and total dark matter halo mass, which are inputs in the 
        Hernquist mass function.

        inputs:
            radius (astropy quantity): the radius at which we wish to find the Hernquist mass, with units of kpc
            a (astropy quantity): scale factor of the galaxy of interest, with units of kpc
            Mhalo (astropy quantity): total dark matter halo mass of the galaxy of interest, with units of Msun

        outputs:
            Hernquist_mass (astropy quantity): a float mass of the Hernquist mass calculated using the Hernquist formula
            at the specified radius, with the scale factor a and dark matter mass Mhalo, in units of Msun
        '''
        Hernquist_mass = Mhalo * radius**2 / (a+radius)**2 #using Hernquist mass function
        return Hernquist_mass
    def CircularVelocity(self, ptype, radii):
        '''
        This function returns the circular velocity using the Newtonian force equation and assuming 
        spherical symmetry. (GMm / r^2 = v^2/r) -> (v = sqrt(GM/r)), where M is the enclosed mass any particle outside
        would be susceptible to the gravity of. The enclosed mass is only enclosed mass of a specified particle type

        inputs:
            ptype (int): the type of particle we want to find the enclosed mass of, either
            1, 2, or 3

            radii (astropy quantity): an array of radii with units of kpc, provided to find the total enclosed mass
            within several distances of the COM in series for the Mass and velocity profile.

        outputs:
            vel_array (astropy array): an array of velocities, based on the total enclosed mass of the specified 
            particle type within the radius at the matching index of the radii input array.
        '''
        Grav = G.to(u.kpc*u.km**2/(u.s**2*u.Msun)) #converting G to our units of distance, mass, and time
        vel_array = np.sqrt(Grav * self.MassEnclosed(ptype, radii)/radii) #using Newtonian and circular velocity force balacing
        return vel_array
    def CircularVelocityTotal(self, radii):
        '''
         This function returns the circular velocity using the Newtonian force equation and assuming 
        spherical symmetry. (GMm / r^2 = v^2/r) -> (v = sqrt(GM/r)), where M is the enclosed mass any particle outside
        would be susceptible to the gravity of. The enclosed mass M is enclosed mass of all particle types in the galaxy.

        inputs:
            radii (astropy quantity): an array of radii with units of kpc, provided to find the total enclosed mass
            within several distances of the COM in series for the Mass and velocity profile.

        outputs:
            vel_array (astropy array): an array of velocities, based on the total enclosed mass of all 
            particle types within the radius at the matching index of the radii input array.
        '''
        Grav = G.to(u.kpc*u.km**2/(u.s**2*u.Msun)) #converting G to our units of distance, mass, and time
        vel_array = np.sqrt(Grav * self.MassEnclosedTotal(radii) / radii) 
        #same force balancing equation, but calling total enclosed mass to do it
        return vel_array

    def HernquistVCirc(self, radius, a, Mhalo):
        '''
        This function returns the circular velocity using the Newtonian force equation and assuming 
        spherical symmetry. (GMm / r^2 = v^2/r) -> (v = sqrt(GM/r)), where M is the Hernquist mass we 
        calculated earlier.
        '''
        Grav = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        velocity = np.sqrt(Grav * self.HernquistMass(radius, a, Mhalo)/ radius)
        return velocity.value #must be unitless in the return array when plotting



if __name__ == '__main__' : 
    MW = MassProfile("MW", 0) #getting MW galaxy data
    M31 = MassProfile("M31", 0) #getting M31 galaxy data
    M33 = MassProfile("M33", 0) #getting M33 galaxy data
    r = np.arange(0.1, 30.1, 0.1) * u.kpc #structuring my radius list with lots of points for smooth curves
    
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, sharey=True)

    MW_Halo_M = ComponentMass("MW_000.txt", 1) * 10**12 #using ComponentMass function to find total mass of dark matter in galaxy
    MW_scale_factor_guess = 63 * u.kpc #guess + assigned units
    MW_Hernquist_Mass = np.array([]) #will be filled
    for i in range(0, len(r)): #need to break the array into one component at a time to find the individual hernquist mass
        M_element = MW.HernquistMass(r[i], MW_scale_factor_guess, MW_Halo_M) #calculating the Hernquist mass element for one radius at a time
        MW_Hernquist_Mass = np.append(MW_Hernquist_Mass, M_element) #filling out MW_Hernquist_Mass
    ax1.semilogy(r, MW.MassEnclosed(1,r), linestyle="dotted", label = 'Halo Mass Profile')
    ax1.semilogy(r, MW.MassEnclosed(2,r), linestyle="--", label = 'Disk Mass Profile')
    ax1.semilogy(r, MW.MassEnclosed(3,r), linestyle=":", label = 'Bulge Mass Profile')
    ax1.semilogy(r, MW.MassEnclosedTotal(r), linestyle="-.", label = 'Total Mass Profile')
    ax1.semilogy(r, MW_Hernquist_Mass, label = 'Hernquist Mass Profile, a = ' + str(MW_scale_factor_guess))
    ax1.set_title("Milky Way Mass Profiles")
    ax1.set_xlabel("Radius (kpc)")
    ax1.set_ylabel("Mass (Msun)")
    ax1.set_ylim(10**4, 10**12) #setting a limit to see the interesting parts of the curve (small radius values have very small mass values)
    ax1.legend()

    ##### basically copy pasting the method to plot all of MW data, but now with M31 values. Most quantities will look the same. 
    ##### Maybe write this method into its own function?
    M31_Halo_M = ComponentMass("M31_000.txt", 1) * 10**12 
    M31_scale_factor_guess = 60 * u.kpc #different guess for M31 as opposed to MW
    M31_Hernquist_Mass = np.array([])
    for i in range(0, len(r)):
        M_element = M31.HernquistMass(r[i], M31_scale_factor_guess, M31_Halo_M)
        M31_Hernquist_Mass = np.append(M31_Hernquist_Mass, M_element)
    ax2.semilogy(r, M31.MassEnclosed(1,r), linestyle="dotted", label = 'Halo Mass Profile')
    ax2.semilogy(r, M31.MassEnclosed(2,r), linestyle="--", label = 'Disk Mass Profile')
    ax2.semilogy(r, M31.MassEnclosed(3,r), linestyle=":", label = 'Bulge Mass Profile')
    ax2.semilogy(r, M31.MassEnclosedTotal(r), linestyle="-.", label = 'Total Mass Profile')
    ax2.semilogy(r, M31_Hernquist_Mass, label = 'Hernquist Mass Profile, a = ' + str(M31_scale_factor_guess))
    ax2.set_title("M31 Mass Profiles")
    ax2.set_xlabel("Radius (kpc)")
    ax2.set_ylabel("Mass (Msun)")
    ax2.set_ylim(10**4, 10**12)
    ax2.legend()

    M33_Halo_M = ComponentMass("M33_000.txt", 1) * 10**12 
    M33_scale_factor_guess = 25 * u.kpc #different guess for M33 scale factor
    M33_Hernquist_Mass = np.array([])
    for i in range(0, len(r)):
        M_element = M33.HernquistMass(r[i], M33_scale_factor_guess, M33_Halo_M)
        M33_Hernquist_Mass = np.append(M33_Hernquist_Mass, M_element)
    ax3.semilogy(r, M33.MassEnclosed(1,r), linestyle="dotted", label = 'Halo Mass Profile')
    ax3.semilogy(r, M33.MassEnclosed(2,r), linestyle="--", label = 'Disk Mass Profile')
    ax3.semilogy(r, M33.MassEnclosedTotal(r), linestyle="-.", label = 'Total Mass Profile')
    ax3.semilogy(r, M33_Hernquist_Mass, color = 'tab:purple', label = 'Hernquist Mass Profile, a = ' + str(M33_scale_factor_guess))
    # note, we skip callnig M33.MassEnclosed(3,r) since M33 has no bulge particles
    ax3.set_title("M33 Mass Profiles")
    ax3.set_xlabel("Radius (kpc)")
    ax3.set_ylabel("Mass (Msun)")
    ax3.set_ylim(10**4, 10**12)
    ax3.legend()
    plt.show()
    
    
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, sharey=True) # a new set of 3 subplots
    MW_Hernquist_V = np.array([])
    for i in range(0, len(r)):
        V_element = MW.HernquistVCirc(r[i], MW_scale_factor_guess, MW_Halo_M) #instead of calling HernquistMass, we now call HernquistVCirc
        MW_Hernquist_V = np.append(MW_Hernquist_V, V_element)
    ax1.semilogy(r, MW.CircularVelocity(1,r), linestyle="dotted", label = 'Halo Rotation Curve') #replacing mass profiles with rotation curves
    ax1.semilogy(r, MW.CircularVelocity(2,r), linestyle="--", label = 'Disk Rotation Curve')
    ax1.semilogy(r, MW.CircularVelocity(3,r), linestyle=":", label = 'Bulge Rotation Curve')
    ax1.semilogy(r, MW.CircularVelocityTotal(r), linestyle="-.", label = 'Total Rotation Curve')
    ax1.semilogy(r, MW_Hernquist_V, label = 'Hernquist Rotation Curve, a = ' + str(MW_scale_factor_guess))
    ax1.set_title("Milky Way Rotation Curves")
    ax1.set_xlabel("Radius (kpc)")
    ax1.set_ylabel("Velocity (km/s)")
    ax1.set_ylim(1, 10**3) #different y limit to optimize viewing
    ax1.legend()

    # repeat above but for M31 instead of MW
    M31_Hernquist_V = np.array([])
    for i in range(0, len(r)):
        V_element = M31.HernquistVCirc(r[i], M31_scale_factor_guess, M31_Halo_M)
        M31_Hernquist_V = np.append(M31_Hernquist_V, V_element)
    ax2.semilogy(r, M31.CircularVelocity(1,r), linestyle="dotted", label = 'Halo Rotation Curve')
    ax2.semilogy(r, M31.CircularVelocity(2,r), linestyle="--", label = 'Disk Rotation Curve')
    ax2.semilogy(r, M31.CircularVelocity(3,r), linestyle=":", label = 'Bulge Rotation Curve')
    ax2.semilogy(r, M31.CircularVelocityTotal(r), linestyle="-.", label = 'Total Rotation Curve')
    ax2.semilogy(r, M31_Hernquist_V, label = 'Hernquist Rotation Curve, a = ' + str(M31_scale_factor_guess))
    ax2.set_title("M31 Rotation Curves")
    ax2.set_xlabel("Radius (kpc)")
    ax2.set_ylabel("Velocity (km/s)")
    ax2.set_ylim(1, 10**3)
    ax2.legend()
    
    #repeat above but for M33 instead of M31 (don't forget about the lack of bulge particles!)
    M33_Hernquist_V = np.array([])
    for i in range(0, len(r)):
        V_element = M33.HernquistVCirc(r[i], M33_scale_factor_guess, M33_Halo_M)
        M33_Hernquist_V = np.append(M33_Hernquist_V, V_element)
    ax3.semilogy(r, M33.CircularVelocity(1,r), linestyle="dotted", label = 'Halo Rotation Curve')
    ax3.semilogy(r, M33.CircularVelocity(2,r), linestyle="--", label = 'Disk Rotation Curve')
    ax3.semilogy(r, M33.CircularVelocityTotal(r), linestyle="-.", label = 'Total Rotation Curve')
    ax3.semilogy(r, M33_Hernquist_V, color = 'tab:purple', label = 'Hernquist Rotation Curve, a = ' + str(M33_scale_factor_guess))
    #bulge particle plot not included for M33 rotation curves
    ax3.set_title("M33 Rotation Curves")
    ax3.set_xlabel("Radius (kpc)")
    ax3.set_ylabel("Velocity (km/s)")
    ax3.set_ylim(1, 10**3)
    ax3.legend()
    plt.show()





        
