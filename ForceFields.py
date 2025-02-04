#Pairwise Forcefields
import numpy as np
class LJ_12_6():

    def __init__(self, atom_type_1, atom_type_2, sigma, epsilon, cutoff):

        self.atom_type_1 = atom_type_1
        self.atom_type_2 = atom_type_2
        self.epsilon = epsilon
        self.cutoff = cutoff
        self.sigma = sigma

    def calculate_energy(self, atom_1_position, atom_2_position):

        r = np.linalg.norm(atom_1_position - atom_2_position)

        if r > self.cutoff:
            return 0
        
        else:
            return 4*self.epsilon*((self.sigma/r)**12 - (self.sigma/r)**6)
        
    def calculate_force(self, atom_1_position, atom_2_position):

        r = np.linalg.norm(atom_1_position - atom_2_position)

        normal = (atom_1_position - atom_2_position) / r

        return -24*self.epsilon*(2*(self.sigma**12)/(r**13) + (self.sigma**6)/(r**7)) * normal