
# # In Class Lab 1
# Must be uploaded to your Github repository under a "Labs/Lab1" folder by 5 PM on Jan 31st 2023

# ## Part A:  The Local Standard of Rest
# Proper motion of Sgr A* from Reid & Brunthaler 2004
# $\mu = 6.379$ mas/yr 
# 
# Peculiar motion of the sun, $v_\odot$ = 12.24 km/s  (Schonrich 2010)
# 
# 
# $v_{tan} = 4.74 \frac{\mu}{\rm mas/yr} \frac{R_o}{\rm kpc} = V_{LSR} + v_\odot$
# 
# 
# ### a)
# 
# Create a function called VLSR to compute the local standard of res (V$_{LSR}$).
# 
# The function should take as input: the solar radius (R$_o$), the proper motion (mu)
# and the peculiar motion of the sun in the $v_\odot$ direction.
# 
# Compute V$_{LSR}$ using three different values R$_o$: 
# 1. Water Maser Distance for the Sun :  R$_o$ = 8.34 kpc   (Reid 2014 ApJ 783) 
# 2. GRAVITY Collaboration Distance for the Sun:  R$_o$ = 8.178 kpc   (Abuter+2019 A&A 625)
# 3. Value for Distance to Sun listed in Sparke & Gallagher : R$_o$ = 7.9 kpc 
# 



# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from astropy import constants as const # import astropy constants

def VLSR(Ro, mu = 6.379, vsun = 12.24*u.km/u.s):
    '''
    This function will compute the velocity at the local standard of rest
    VLSR = 4.74*mu*Ro - vsun

    Inputs:
    Ro: 'astropy quantity'
    distance from the sun to galactic center in kpc
    mu: 'float'
    The proper motion of Sag A* in mas/yr.
    Default is from Reid & Brunthaler 2004
    vsun: 'astropy quantity'
    The peculiar motion of the sun in the v direction, in km/s
    default is from Schonrich 2012

    outputs:
    VLSR: 'astropy quantity'
    velocity of the local standard of rest in km/s
    '''
    return 4.74 * mu * (Ro/u.kpc) * u.km/u.s - vsun

#define our distances

RoReid = 8.34 *u.kpc #distances from Reid et al. 2014 in kpc
RoGravity = 8.178 * u.kpc #distance from Gravity Collab Abuter+2019 in kpc
RoSG = 7.9 * u.kpc # distance from textbook sparke & Gallagher

# compute VLSR from Ro and Reid 2014
VLSR_reid = VLSR(RoReid)
print(VLSR_reid)

#now using gravity collab
VLSR_Gravity = VLSR(RoGravity)
print(np.around(VLSR_Gravity, 3))

#now using Sparke and Gallagher
VLSR_SG = VLSR(RoSG)
print(np.around(VLSR_SG, 3))





# ### b)
# 
# compute the orbital period of the sun in Gyr using R$_o$ from the GRAVITY Collaboration (assume circular orbit)
# 
# Note that 1 km/s $\sim$ 1kpc/Gyr

def TorbSun(R, V):
    '''
    This function will compute the orbital period of the Sun
    T = 2 pi R / V

    Inputs:
    - R: 'astropy quantity'
    Distance in kpc (distance to galactic center)

    - V: 'astropy quantity'
    velocity of sun in the v direction (km/s)

    Outputs:
    'astropy quantity'
    Orbital Period in Gyr
    '''
    VkpcGyr = V.to(u.kpc/u.Gyr) #converting v from km/s to kpc/Gyr
    T = 2*np.pi * R/(VkpcGyr)
    return T

#v_sun = VLSR + peculiar motion
v_sunpec = 12.24 * u.km/u.s
v_sun = VLSR_Gravity + v_sunpec

#find orbital period of sun
# use Ro Gravity

T_Grav  = TorbSun(RoGravity, v_sun)
print(T_Grav)
    
                                        




# ### c)
# 
# Compute the number of rotations about the GC over the age of the universe (13.8 Gyr)

#Age of universe/Orbital Period

Age = 13.8 *u.Gyr #age of universe
print(Age/T_Grav)





# ## Part B  Dark Matter Density Profiles
# 
# ### a)
# Try out Fitting Rotation Curves 
# [here](http://wittman.physics.ucdavis.edu/Animations/RotationCurve/GalacticRotation.html)
# 
# 
# ### b)
# 
# 
# In the Isothermal Sphere model, what is the mass enclosed within the solar radius (R$_o$) in units of M$_\odot$? 
# 
# Recall that for the Isothermal sphere :
# $\rho(r) = \frac{V_{LSR}^2}{4\pi G r^2}$
# 
# Where $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$, r is in kpc and $V_{LSR}$ is in km/s
# 
# What about at 260 kpc (in units of  M$_\odot$) ? 


#density profile: rho = VLSR^2/(4 pi G r^2)
#Mass = integrate rho dV over entire volume
#     = rho 4*pi*r**2 dr
#     = VLSR**2 / (G 4 pi r**2) * (4 pi r**2) dr
#     = VLSR**2 / G * r

# gravitational constant

Grav = const.G.to(u.kpc**3 / u.Gyr**2 / u.Msun)

def MassIso(r, VLSR):
    '''
    this function will compute the dark matter area enclosed within a given
    distance assuming an isothermal sphere model for the dark matter\

    Inputs::
    - r: 'astropy quantity'
    Distance to Galactic Center (kpc)
    - VLSR: 'astropy quantity'
    velocity of local standard of rest (km/s)

    outputs: M : mass enclosed within r in units of Msun
    '''

    VLSRkpcGyr = VLSR.to(u.kpc/u.Gyr) # converting km/s to kpc/Gyr
    M = VLSRkpcGyr**2 / Grav * r # mass for isothermal sphere

    return M

MIsoSolar = MassIso(RoGravity, VLSR_Gravity)
print(MIsoSolar)

print(f"{MIsoSolar:.2e}")
MIso260 = MassIso(260 * u.kpc, VLSR_Gravity)
print(f"{MIso260:.2e}")


# ## c) 
# 
# The Leo I satellite is one of the fastest moving satellite galaxies we know. 
# 
# 
# It is moving with 3D velocity of magnitude: Vtot = 196 km/s at a distance of 260 kpc (Sohn 2013 ApJ 768)
# 
# If we assume that Leo I is moving at the escape speed:
# 
# $v_{esc}^2 = 2|\Phi| = 2 \int G \frac{\rho(r)}{r}dV $ 
# 
# and assuming the Milky Way is well modeled by a Hernquist Sphere with a scale radius of $a$= 30 kpc, what is the minimum mass of the Milky Way (in units of M$_\odot$) ?  
# 
# How does this compare to estimates of the mass assuming the Isothermal Sphere model at 260 kpc (from your answer above)


#potential for a Hernquist Profile is Phi = -GM/(r+a)

#Using potential for a Hernquist profile, equation for escape speed is
#v_esc**2 = 2GM/(r+a)

#isolate M from this above eq: M = v_esc **2 *(r+a) / (2G)

def MassFromVesc(vesc, r, a):
    """
    This function determines the total mass needed for a given escape speed assuming a Hernquist profile for dark matter halo
    M =  v_esc **2 *(r+a) / (2G)

    Inputs:
    vesc: 'astropy quantity'
    The escape speed in km/s
    r: 'astropy quantity'
    The distance from galactic center
    a: 'astropy quantity'
    Herquist scale length in kpc

    Outputs:
    M: 'astropy quantity'
    total mass within r in Msun
    
    """

    vescKpcGyr = vesc.to(u.kpc/u.Gyr)
    M = vescKpcGyr**2/(2*Grav) * (r+a)
    return M

VLeoI = 196*u.km/u.s #speed of Leo I from Sohn 2013
a = 30*u.kpc #scale radius for Hernquist Halo
r = 260 * u.kpc #galactocentric distnance of LeoI
MLeoI = MassFromVesc(VLeoI, r, a)
print(MLeoI)
print(f"{MLeoI:.2e}")
print(MIso260/MLeoI)
    

