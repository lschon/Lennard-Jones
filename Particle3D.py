"""PARTICLE 3D by Lucas Schön (s1718492). 28 November 2019. Exercise 2 Computer Modelling.
Aim: to create a class which provides the position, velocity, separation, and energy details
of a 3D particle given text files that hold the initial conditions. This code is to be referenced
in another time integration program and will make use of the functions provided here at each time step"""

import numpy as np
import csv

class P3D(object):

    def __init__(self, x, y, z, vx, vy, vz, mass, label):
        "Initialise a Particle3D instance"
        "param x,y,z: floats of positions in x, y, and z directions combined into numpy array"
        "param vx,vy,vz: floats of velocity in x, y, and z directions combined into numpy array"
        "param mass: mass as float"
        "param label: label as string"
        self.position=np.array([x, y, z])
        self.velocity=np.array([vx, vy, vz])
        self.mass=float(mass)
        self.label=str(label)

    def __str__(self):
        "Define output format."
        return str(self.label) + " " + str(self.position[0]) + " " +str(self.position[1]) + " " +str(self.position[2])

    def KE(self):
        "Return kinetic energy as 1/2*mass*vel^2"
        return 0.5*self.mass*(np.sum(np.square(self.velocity)))

    def leap_velocity(self, dt, force):
        "First-order velocity update, v(t+dt) = v(t) + dt*F(t)"
        "param dt: timestep as float"
        "param force: force on particle as float"
        self.velocity += dt*force/self.mass

    def leap_pos_1st(self, dt):
        "First-order position update, x(t+dt) = x(t) + dt*v(t)"
        ":param dt: timestep as float"
        self.position += dt*self.velocity

    def leap_pos_2nd(self, dt, force):
        "Second-order position update, x(t+dt) = x(t) + dt*v(t) + (1/2m)*dt^2*F(t)"
        "param dt: timestep as float"
        "param force: current force as float"
        self.position += dt*self.velocity + 0.5*dt**2*force/self.mass


    @staticmethod
    def parameterfilereader(parameterfile):
        line=parameterfile.readline()
        tokens=line.split(",")
        e_k_b=float(tokens[0])
        sigma=float(tokens[1])
        m=float(tokens[2])
        N=float(tokens[3])
        temp=float(tokens[4])
        rho=float(tokens[5])
        return e_k_b, sigma, m, N, temp, rho

    @staticmethod
    def separation(p1,p2):
        """Static method for determining the separation between two particles"""
        return (p1.position-p2.position)
