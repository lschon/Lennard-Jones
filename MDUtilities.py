"""
CMod Project B: auxiliary MD methods
"""

import random
import numpy as np
from Particle3D import P3D


def set_initial_positions(rho,particles):

    # Determine number of particles
    natoms = len(particles)

    # Set box dimensions
    box_size = (natoms/rho)**(1./3.)

    # Number or particles in each direction
    ndim = int(float((natoms-1)/4.0)**(1./3.))+1

    # Give warning if fcc lattice will not be fully occupied
    if 4*ndim**3 != natoms:
        print("Atoms will not fill a fcc lattice completely.\n")

    # Separation between particles
    delta = box_size / ndim

    # Set particle positions
    i_atom = 0
    for ix in range(ndim):
        for iy in range(ndim):
            for iz in range(ndim):
                if i_atom<natoms:
                    pos = np.array([ix*delta, iy*delta, iz*delta])
                    particles[i_atom].position = pos
                    i_atom += 1
                if i_atom<natoms:
                    pos = np.array([(ix+0.5)*delta, (iy+0.5)*delta, iz*delta])
                    particles[i_atom].position = pos
                    i_atom += 1
                if i_atom<natoms:
                    pos = np.array([(ix+0.5)*delta, iy*delta, (iz+0.5)*delta])
                    particles[i_atom].position = pos
                    i_atom += 1
                if i_atom<natoms:
                    pos = np.array([ix*delta, (iy+0.5)*delta, (iz+0.5)*delta])
                    particles[i_atom].position = pos
                    i_atom += 1

    # Some output
    print(N, "atoms placed on a face-centered cubic lattice.")
    print("Box dimensions: ", box_size, box_size, box_size)
    box=np.array([box_size, box_size, box_size])
    # Return the box size as Numpy array and the particle positions

    return box, particles


def set_initial_velocities(temp, particles):

    # Determine number of particles
    natoms = len(particles)

    # Zero the accumulators
    xv0 = 0.0
    yv0 = 0.0
    zv0 = 0.0
    vsq = 0.0

    # Loop over particles, set velocities
    for i in range(natoms):
        # Random inital velocities
        xvt = random.random() - 0.5
        yvt = random.random() - 0.5
        zvt = random.random() - 0.5

        particles[i].velocity = np.array([xvt, yvt, zvt])

        # Add to total velocity
        xv0 += xvt
        yv0 += yvt
        zv0 += zvt
        vsq += xvt**2 + yvt**2 + zvt**2

        print(xvt,yvt,zvt)

    # Centre-of-mass motion
    xv0 /= natoms
    yv0 /= natoms
    zv0 /= natoms

    # Boltzmann factor
    kB = (3*natoms*temp/vsq)**(1./2.)
    print("kB is: ", kB)
    # Zero the probe accumulators
    xv0_tot = 0.0
    yv0_tot = 0.0
    zv0_tot = 0.0
    v0sq = 0.0

    # Rescale all velocities
    for i in range(natoms):
        vtemp = particles[i].velocity
        print("vtemp is: ", vtemp)
        xvt = kB*(vtemp[0] - xv0)
        yvt = kB*(vtemp[1] - yv0)
        zvt = kB*(vtemp[2] - zv0)
    # Note to self: After this step, the component velocities of the particles add up to zero every time
    # this does not make sense
        print(xvt,yvt,zvt)

        particles[i].velocity = np.array([xvt, yvt, zvt])

        xv0_tot += xvt
        yv0_tot += yvt
        zv0_tot += zvt
        v0sq += xvt**2 + yvt**2 + zvt**2

    # Output
    com_vel=[xv0_tot/natoms, yv0_tot/natoms, zv0_tot/natoms]
    print("Temperature = ", temp, "Kelvin")
    print("Centre-of-mass velocity is: ", com_vel)

    return particles
