from ovito.data import CutoffNeighborFinder, DataCollection
from ovito.pipeline import ModifierInterface
from traits.api import Int, Float
import numpy as np 

class CalculateHydrogenBonds(ModifierInterface):
    donor = Int()
    accep = Int()
    H_atom = Int()
    da_rcut = Float()
    dha_rangle = Float()

    def __init__(self, donor, accep, H_atom, da_rcut, dha_rangle):
        self.donor = donor
        self.accep = accep
        self.H_atom = H_atom
        self.da_rcut = da_rcut
        self.dha_rangle = dha_rangle

    # Function to calculate periodic boundary conditions (PBC)-corrected distances
    def pbc_distance(self, pos1, pos2, box_lengths):
        half_box_length = box_lengths / 2  # Half the box length for wrapping
        delta = pos1 - pos2  # Calculate raw distance
        # Correct distances that exceed half the box length
        delta = np.where(np.abs(delta) > half_box_length, delta - np.sign(delta) * box_lengths, delta)
        return delta
    
    # Function to calculate angle (donor-hydrogen-acceptor) with PBC handling
    def calculate_angle(self, d, h, a, box_lengths):
        vec_hd = self.pbc_distance(h, d, box_lengths)  # Vector from hydrogen to donor
        vec_ha = self.pbc_distance(h, a, box_lengths)  # Vector from hydrogen to acceptor
        # Compute the cosine of the angle using dot product
        cos_ang = np.dot(vec_hd, vec_ha) / (np.linalg.norm(vec_hd) * np.linalg.norm(vec_ha))
        angle_d = np.arccos(cos_ang)  # Calculate the angle in radians
        angle_d_deg = np.round(np.degrees(angle_d), 3)  # Convert to degrees and round to 3 decimals
        return angle_d_deg
    
    def modify(self, data, **kwargs):
        finder = CutoffNeighborFinder(self.da_rcut, data)  # Initialize a neighbor finder with a cutoff distance
        particle_positions = data.particles.position  # Positions of all particles
        particle_types = data.particles.particle_types  # Types of all particles
        box_lengths = np.array([data.cell[0][0], data.cell[1][1], data.cell[2][2]])
        
        # Identify donor and acceptor particles
        donor_indices = np.where(particle_types == self.donor)[0]
        acceptor_indices = np.where(particle_types == self.accep)[0]
        
        oh_rcut = 1.05  # Cutoff for donor-hydrogen bonds in H2O
        nhbond = 0  # Initialize hydrogen bond count
        h_donor_array = np.zeros(data.particles.count, dtype=int)
        
        for donor_index in donor_indices:
            donor_O = particle_positions[donor_index]  # Get donor oxygen position
    
            # Find neighbors of the donor oxygen
            neighbor_indices = np.array([neigh.index for neigh in finder.find(donor_index)])
            neighbor_distances = np.array([neigh.distance for neigh in finder.find(donor_index)])
            neighbor_positions = particle_positions[neighbor_indices]
    
            # Identify donor hydrogens (distance ≤ 1.0 Å) and acceptor atoms
            donor_hyd_mask = neighbor_distances <= oh_rcut
            acceptor_mask = (particle_types[neighbor_indices] == self.accep) & (neighbor_distances > oh_rcut)
            donor_hyd_positions = neighbor_positions[donor_hyd_mask]  # Positions of donor hydrogens
            acceptor_positions = neighbor_positions[acceptor_mask]  # Positions of acceptor atoms
            
            # Count hydrogen bonds based on angle and distance criteria
            angles = len([angle for donor_hyd in donor_hyd_positions
                       for acceptor in acceptor_positions
                       if (angle := self.calculate_angle(donor_O, donor_hyd, acceptor, box_lengths)) >= self.dha_rangle])
            nhbond += angles  
            h_donor_array[donor_index] = 1 if angles > 0 else 0
        
        # Store the hydrogen bond count as an attribute of the current frame
        data.attributes["Hbond_count"] = nhbond
        data.attributes["donor_count"] = len(donor_indices)
         # Assign the `HDonor` property to all particles
        data.particles_.create_property("HDonor", data=h_donor_array)
