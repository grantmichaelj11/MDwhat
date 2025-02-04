
import numpy as np
from itertools import combinations

class VelocityVerlet():

    def __init__(self, dt, atoms, pairwise_list):

        self.dt = dt
        self.atoms = atoms
        self.pairwise_list = pairwise_list

    def iterate(self, timesteps, abinitio=False):

        per_atom_forces = np.zeros((len(self.atoms), 3))

        for step in range(timesteps):
            #If using Classical Molecular Dynamics
            if not abinitio:
                if step == 0:
                    #Calculate Pairwise Forces
                    for i in range(len(self.pairwise_list)):

                        atom_id_1 = self.pairwise_list[i][0]['id']
                        atom_position_1 = np.array([self.pairwise_list[i][0]['x'], self.pairwise_list[i][0]['y'], self.pairwise_list[i][0]['z']])

                        atom_id_2 = self.pairwise_list[i][1]['id']
                        atom_position_2 = np.array([self.pairwise_list[i][1]['x'], self.pairwise_list[i][1]['y'], self.pairwise_list[i][1]['z']])

                        # Uses forcefield object to calculate force on atom 1
                        initial_force_update = self.pairwise_list[i][2].calculate_force(atom_position_1, atom_position_2)

                        # Update Forces (Use 3rd Newtons Law)

                        per_atom_forces[atom_id_1-1] += initial_force_update
                        per_atom_forces[atom_id_2-1] -= initial_force_update

                #Update Positions
                for i in range(len(self.atoms)):

                    self.atoms[i]['x'] += self.dt*self.atoms[i]['vx'] + (self.dt**2)/(2*self.atoms[i]['mass'])*per_atom_forces[i][0]
                    self.atoms[i]['y'] += self.dt*self.atoms[i]['vy'] + (self.dt**2)/(2*self.atoms[i]['mass'])*per_atom_forces[i][1]
                    self.atoms[i]['z'] += self.dt*self.atoms[i]['vz'] + (self.dt**2)/(2*self.atoms[i]['mass'])*per_atom_forces[i][2]

                # Update Pairwise List
                new_pairwise_list = list(combinations(self.atoms,2))

                for i in range(len(self.pairwise_list)-1):
                    self.pairwise_list[i][0] = new_pairwise_list[i][0]
                    self.pairwise_list[i+1][1] = new_pairwise_list[i+1][1]

                #Calculate Force (At next timestep)
                per_atom_forces_dt = np.zeros((len(self.atoms), 3))

                for i in range(len(self.pairwise_list)):

                    atom_id_1 = self.pairwise_list[i][0]['id']
                    atom_position_1 = np.array([self.pairwise_list[i][0]['x'], self.pairwise_list[i][0]['y'], self.pairwise_list[i][0]['z']])

                    atom_id_2 = self.pairwise_list[i][1]['id']
                    atom_position_2 = np.array([self.pairwise_list[i][1]['x'], self.pairwise_list[i][1]['y'], self.pairwise_list[i][1]['z']])

                    # Uses forcefield object to calculate force on atom 1
                    initial_force_update = self.pairwise_list[i][2].calculate_force(atom_position_1, atom_position_2)

                    # Update Forces (Use 3rd Newtons Law)

                    per_atom_forces_dt[atom_id_1-1] += initial_force_update
                    per_atom_forces_dt[atom_id_2-1] -= initial_force_update

                #Update Velocities
                for i in range(len(self.atoms)):

                    self.atoms[i]['vx'] += self.dt/(2*self.atoms[i]['mass'])*(per_atom_forces[i][0] + per_atom_forces_dt[i][0])
                    self.atoms[i]['vy'] += self.dt/(2*self.atoms[i]['mass'])*(per_atom_forces[i][1] + per_atom_forces_dt[i][1])
                    self.atoms[i]['vz'] += self.dt/(2*self.atoms[i]['mass'])*(per_atom_forces[i][2] + per_atom_forces_dt[i][2])

                # Update Pairwise List
                new_pairwise_list = list(combinations(self.atoms,2))

                for i in range(len(self.pairwise_list)-1):
                    self.pairwise_list[i][0] = new_pairwise_list[i][0]
                    self.pairwise_list[i+1][1] = new_pairwise_list[i+1][1]

                #The next timesteps forces are updated
                per_atom_forces = per_atom_forces_dt

            return self.atoms, self.pairwise_list

        #If using DFT based Abinitio
        if abinitio:

            pass


            