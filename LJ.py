"""
Lennard Jones (LJ) by Lucas Schön (s1718492). February 2020.
Lennard Jones Project Computer Modelling.

Velocity Verlet time integration of two bonded particles,
with the aim of determining
vibrational frequency of the bond.

Produces plots of the position of the particles
and the energy, both as function of time. Also
saves both to file.

Uses the LJ potential (and by extension,
force) with user-provided parameters of XXX
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3D import P3D
import MDUtilities

def parameters():
    # Open file with parameters
    parameterfile=open("LJparameters.txt", "r")
    list=P3D.parameterfilereader(parameterfile)
    e_k_b=float(list[0])
    sigma=float(list[1])
    m=float(list[2])
    N=int(list[3])
    temp=float(list[4])
    rho=float(list[5])
    epsilon=e_k_b*(1.38064852*10**(-23))
    return e_k_b, sigma, m, N, temp, rho, epsilon

def initialiser():
    # Sets particle objects with appropriate positions and randomised velocities
    _, _, _, N, temp, rho, _ = parameters()
    particles=[P3D(0, 0, 0, 0, 0, 0, 1, "particle" + str(i)) for i in range(N)]
    box, particles= MDUtilities.set_initial_positions(N, rho, particles)
    MDUtilities.set_initial_velocities(temp, particles)

    return particles, box

 # No need for pbc, only mic
def pbc():
    particles, box, N = initialiser()
    for i in range(N):
        particles[i].position=np.mod(particles[i].position,box)
    return particles # May not even need this return particles according to Miguel

def mic():
    particles, box, N = initialiser()
    for i in range(N):
        particles[i].position=np.mod(particles[i].position+box/2,box)-box/2
    print(particles.position)
    return particles

def LJforce(rsep, rmag):
    """
    Returns force acting on particle as Numpy array
    """
    if rmag < 2.5 :
        force= 48*((1/rmag**14)-(1/(2*rmag**8)))*rsep
    else:
        force = 0


    return force

def LJpot(r_red):
    """
    Returns potential energy of particle as float
    """
    if rmag_red < 2.5:
        pot=4*(1/rmag_red**(12)-1/rmag_red**6)
    else:
        pot = 0

        return pot

"""
Some rough work

"""

e_k_b, sigma, m, N, temp, rho, epsilon = parameters()
particles, box =initialiser()

#Quite possibly works -- marvel at my genius
#note that what prints is the force acting on each particle ie should get N forces


#verlet integrator
#we have a list of particle seperations sep_list
#we have a the forces acting on each particle total_force







"""
"""

"""

def MSD():
    MSD=[]
    for i in range(N):
        displacement= particles[i].position-

    MSD=
    1/N=
    return x

def RDF():

    return x

def pot_tot():

    return x

def KE_tot():

    return


"""
#Some rough work


e_k_b, sigma, m, N, temp, rho, epsilon = parameters()
particles, box =initialiser()

#Quite possibly works -- marvel at my genius
#note that what prints is the force acting on each particle ie should get N forces
def totalforce():
    for j in range(N):
        force_list = []
        sep_list = []
        for i in range(N-1):
            if (i != j):
                rsep = (P3D.separation(particles[j],particles[i]))
                rmag = np.linalg.norm(rsep)
                force = LJforce(rsep,rmag)
                force_list.append(force)
        total_force = sum(force_list)
        print(total_force)
        return total_force
totalforce()
#verlet integrator
#we have a list of particle seperations sep_list
#we have the forces acting on each particle total_force
def timeintegrator():
    dt = 0.01
    time = 0.0
    while time<1:
        for k in range (N):
            # Update particle position
            particles[k].leap_pos_2nd(dt, force)

            # Update force
            force_new = totalforce()

            print(force_new)

            # Update particle velocity by averaging current and new forces
            particles[k].leap_velocity(dt, 0.5*(force+force_new))

            # Re-define separation after change in position
            r_red = P3D.separation(p1,p2)

            # Re-define force value
            force = force_new

            # Increase time
            time += dt

    return particles


# Begin main code
def simulation():

    # Read name of output file from command line
    if len(sys.argv)!=2:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <output file>")
        quit()
    else:
        outfile_name = sys.argv[1]


    _, _, _, N, temp, rho, _ = parameters()

    #Intialise Particle Seperations

    for i in range (N):
        rsep=P3D.separation(particle1,particle2)
        rmag = np.linalg.norm(rsep)
    print(rsep)
    print(rmag)



    # Define reduced units
    rsep_red=rsep/sigma
    rmag_red=rmag/sigma
    E_red=E/epsilon
    T_red=e_k_b
    t_red=sigma*sqrt(m/epsilon)

    # Open output file
    outfile = open(outfile_name, "w")

    #Define initial force
    energy = p1.KE() + p2.KE() + LJpot(r)
    outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time,r,energy))

    # Define initial force
    force = LJforce(rsep_red, r_red)

    # Plot particle separation to screen
    pyplot.title('Velocity Verlet: Separation vs time')
    pyplot.xlabel('Time (10.18 fs)')
    pyplot.ylabel('Separation (Å)')
    pyplot.plot(time_list, pos_list, label="particle separation")
    pyplot.legend(loc="upper left")
    pyplot.show()

    # Plot particle energy to screen
    pyplot.title('Velocity Verlet: total energy vs time')
    pyplot.xlabel('Time (10.18 fs)')
    pyplot.ylabel('Energy (eV)')
    pyplot.plot(time_list, energy_list, label="particle energy")
    pyplot.legend(loc="upper left")
    pyplot.show()
