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

N=4
particles=[P3D(0, 0, 0, 0, 0, 0, 1, "particle" + str(i)) for i in range(N)]
box=MDUtilities.set_initial_positions()

def pbc(box):
    # Finding the projected position within box of sides l
    for i in range(N):
        particles_pos=particles[i].position
        print(particles_pos)
        pbc_pos=np.mod(particles_pos[i],box)
    print("The image of x projected onto the periodic box with range 0 to l is ", pbc_pos)
    return pbc_pos
pbc(box)

"""def mic(box):
    # INCOMPLETE, TRYING TO MAKE PBC WORK FIRST BEFORE THIS ONE.
    # Finding the projected x within box centered on origin
    particles_pos=[particles.position for i in range (N)]
    mic_pos=[]
    for i in range(len(particles)):
        append.mic_pos[i]=np.mod(particles_pos[i]+box/2,box)-box/2
    print("The image of x projected onto the periodic box centered at the origin with range -l/2 to l/2 is  ", pxo)
    return mic_pos"""

def LJforce(rsep_red, r_red):
    """
    Returns force acting on particle as Numpy array
    """
    force=48*(1/r_red**(14)-1/(2*r_red**8))*rsep_red

    return force

def LJpot(r_red):
    """
    Returns potential energy of particle as float
    """

    pot=4*(1/r_red**(12)-1/r_red**6)

    return pot

# Begin main code
def main():
    # Read name of output file from command line
    if len(sys.argv)!=2:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <output file>")
        quit()
    else:
        outfile_name = sys.argv[1]

    # Open file with parameters
    parameterfile=open("LJparameters.txt", "r")
    list=P3D.parameterfilereader(parameterfile)
    e_k_b=float(list[0])
    sigma=float(list[1])
    m=float(list[2])
    epsilon=e_k_b*(1.38064852*10**(-23))

    # Write out initial conditions
    # Need to find separation between N particles including boundary conditions
    # From this point onwards, we need to change a 2 particle system to an N particle system
    rsep=P3D.separation(p1,p2)
    r = np.linalg.norm(rsep)

    # Define reduced units
    r_red=r/sigma
    rsep_red=rsep/sigma
    E_red=E/epsilon
    T_red=e_k_b
    t_red=sigma*sqrt(m/epsilon)

    # Open output file
    outfile = open(outfile_name, "w")

    # Set up simulation parameters
    dt = 0.078
    time = 0.0


    #Define initial force
    energy = p1.KE() + p2.KE() + LJpot(r)
    outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time,r,energy))

    # Define initial force
    force = LJforce(rsep_red, r_red)

    # Initialise data lists for plotting later
    time_list = [time]
    pos_list = [r_red]
    energy_list = [energy]


    # Start the time integration loop
    while time<=20:

        # Update particle position
        p1.leap_pos_2nd(dt, force)
        p2.leap_pos_2nd(dt, -force)

        # Update force
        force_new = LJforce(rsep_red, r_red)

        # Update particle velocity by averaging current and new forces
        p1.leap_velocity(dt, 0.5*(force+force_new))
        p2.leap_velocity(dt, -0.5*(force+force_new))

        # Re-define separation ater change in position
        r_red = P3D.separation(p1,p2)

        # Re-define force value
        force = force_new

        # Increase time
        time += dt

        # Output particle information
        energy = p1.KE() +p2.KE() + LJforce(rsep_red, r_red)
        outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time, r_red, energy))

        # Append information to data lists
        time_list.append(time)
        pos_list.append(r_red)
        energy_list.append(energy)


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

# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
