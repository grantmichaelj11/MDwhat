import numpy as np
from itertools import combinations

class CreateBox():

    def __init__(self, Lx, Ly, Lz, periodicity='periodic'):

        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz

        self.periodicity = 'periodic'

        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedral = []
        self.imporpers = []
        self.coulomb = []
        self.masses = {'H': 1.0078, 'He': 4.0026, 'Li': 6.9410, 'Be': 9.0122, 'B': 10.811, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180}


    def insert_atom(self, x, y, z, atom_type):

        """Places an atom in a specified location with no velocity"""

        atom_mass = self.masses[atom_type] 

        # atom_type, mass, x, y, z, velocity
        self.atoms.append({'id': len(self.atoms)+1, 'atom': atom_type, 'mass': atom_mass, 'x': x, 'y': y, 'z': z, 'vx': 0, 'vy': 0, 'vz': 0, 'fx': 0, 'fy': 0, 'fz': 0})

    def insert_molecule_from_xyz(self, xyz_file):

        pass

    def insert_molecule_from_smiles(self, smiles_file):

        pass

class InitializeSimulation():

    def __init__(self, box, forcefields, dt=2):

        assert isinstance(box, CreateBox), "You need a valid 'CreateBox()' object to initialize the simulation"

        #Build Pairwise Interaction List
        self.details = box
        self.pairwise_list = list(combinations(box.atoms, 2))
        self.dt = dt #femtoseconds

        for forcefield in forcefields:
            for i in range(len(self.pairwise_list)):
                if self.pairwise_list[i][0]['atom'] == forcefield.atom_type_1 and self.pairwise_list[i][1]['atom'] == forcefield.atom_type_2:
                    self.pairwise_list[i] = np.append(self.pairwise_list[i], forcefield)

        self.pairwise_list = np.array(self.pairwise_list)

    def set_velocities(self, box, T_desired):

        #Randomly initialize velocities
        for i in range(len(self.details.atoms)):
            self.details.atoms[i]['vx'] = np.random.uniform(-1, 1)
            self.details.atoms[i]['vy'] = np.random.uniform(-1, 1)
            self.details.atoms[i]['vz'] = np.random.uniform(-1, 1)

        #Calculate the current Temperature
        T = 0
        Nf = 3 * len(box.atoms) - 3
        kb = 0.001987 #kcal/mol * K
        for i in range(len(box.atoms)):
            
            T += self.details.atoms[i]['mass']/Nf/kb*(self.details.atoms[i]['vx']**2 + self.details.atoms[i]['vy']**2 + self.details.atoms[i]['vz']**2)

        #Rescale Temperatures
        for i in range(len(self.details.atoms)):
            self.details.atoms[i]['vx'] *= np.sqrt(T_desired/T)
            self.details.atoms[i]['vy'] *= np.sqrt(T_desired/T)
            self.details.atoms[i]['vz'] *= np.sqrt(T_desired/T)
    

    def run(self, integrator, timesteps):

        self.details.atoms, self.pairwise_list = integrator.iterate(timesteps)
        
        
    