import os
from pymol import cmd
import math

class PISAPlugin:
    def __init__(self):
        self.output_directory = r"C:\path\to\output"
        self.output_filename = "interactions_list.txt"
        self.output_path = os.path.join(self.output_directory, self.output_filename)
        self.distance_cutoff = 4.0
        self.interactions = {
            "hydrophobic": [], "hbond": [], "salt_bridge": [], "covalent": [],
            "vdw": [], "aromatic": [], "cation_pi": [], "other": []
        }

    def analyze_interactions(self, rec_file, lig_file):
        print(f"Running interaction analysis between {rec_file} and {lig_file}...")
        cmd.select("close_contacts", f"({rec_file} within {self.distance_cutoff} of {lig_file}) or ({lig_file} within {self.distance_cutoff} of {rec_file})")
        cmd.select("interacting_residues", "byres close_contacts")

        rec_atoms = cmd.get_model(f"{rec_file} and close_contacts").atom
        lig_atoms = cmd.get_model(f"{lig_file} and close_contacts").atom

        for rec_atom in rec_atoms:
            for lig_atom in lig_atoms:
                dist = cmd.get_distance(f"{rec_file}//{rec_atom.chain}/{rec_atom.resi}/{rec_atom.name}",
                                        f"{lig_file}//{lig_atom.chain}/{lig_atom.resi}/{lig_atom.name}")
                if dist <= self.distance_cutoff:
                    self.classify_interaction(rec_atom, lig_atom, dist)

        self.write_results()
        self.visualize_interactions(rec_file, lig_file)

    def classify_interaction(self, rec_atom, lig_atom, dist):
        if self.is_hbond(rec_atom, lig_atom, dist):
            self.interactions["hbond"].append((rec_atom, lig_atom, dist))
        elif self.is_salt_bridge(rec_atom, lig_atom, dist):
            self.interactions["salt_bridge"].append((rec_atom, lig_atom, dist))
        elif self.is_hydrophobic(rec_atom, lig_atom, dist):
            self.interactions["hydrophobic"].append((rec_atom, lig_atom, dist))
        elif self.is_aromatic(rec_atom, lig_atom, dist):
            self.interactions["aromatic"].append((rec_atom, lig_atom, dist))
        elif self.is_cation_pi(rec_atom, lig_atom, dist):
            self.interactions["cation_pi"].append((rec_atom, lig_atom, dist))
        elif self.is_covalent(dist):
            self.interactions["covalent"].append((rec_atom, lig_atom, dist))
        elif self.is_vdw(rec_atom, lig_atom, dist):
            self.interactions["vdw"].append((rec_atom, lig_atom, dist))
        else:
            self.interactions["other"].append((rec_atom, lig_atom, dist))

    def write_results(self):
        with open(self.output_path, "w") as f:
            for interaction_type, pairs in self.interactions.items():
                f.write(f"{interaction_type.upper()} INTERACTIONS:\n")
                for atom1, atom2, dist in pairs:
                    f.write(f"{atom1.chain}/{atom1.resn}`{atom1.resi}/{atom1.name} -- {atom2.chain}/{atom2.resn}`{atom2.resi}/{atom2.name} : {dist:.2f} Å\n")
                f.write("\n")
        print(f"Results saved to file: {self.output_path}")

    def visualize_interactions(self, rec_file, lig_file):
        cmd.show("sticks", "interacting_residues")
        for interaction_type, color in zip(["hbond", "salt_bridge", "covalent", "hydrophobic"],
                                           ["cyan", "magenta", "yellow", "orange"]):
            if self.interactions[interaction_type]:
                for atom1, atom2, dist in self.interactions[interaction_type]:
                    cmd.distance(f"{interaction_type}_bonds",
                                 f"{rec_file}//{atom1.chain}/{atom1.resi}/{atom1.name}",
                                 f"{lig_file}//{atom2.chain}/{atom2.resi}/{atom2.name}")
                cmd.show("dashes", f"{interaction_type}_bonds")
                cmd.color(color, f"{interaction_type}_bonds")

    def is_hydrophobic(self, atom1, atom2, dist):
        hydrophobic_atoms = {"C", "S"}
        return dist <= 4.0 and atom1.name[0] in hydrophobic_atoms and atom2.name[0] in hydrophobic_atoms

    def is_hbond(self, atom1, atom2, dist):
        if dist > 3.5:
            return False
        hbond_donors_acceptors = {"N", "O", "F"}
        if atom1.name[0] in hbond_donors_acceptors and atom2.name[0] in hbond_donors_acceptors:
            donor_hydrogen = self.find_hydrogen_near(atom1)
            if donor_hydrogen:
                angle = self.calculate_angle(atom1, donor_hydrogen, atom2)
                return angle > 120
        return False

    def is_salt_bridge(self, atom1, atom2, dist):
        positive_atoms = {"NZ", "NH1", "NH2", "NE", "ND1", "ND2"}
        negative_atoms = {"OD1", "OD2", "OE1", "OE2"}
        return dist <= 4.0 and (
            (atom1.name in positive_atoms and atom2.name in negative_atoms) or
            (atom1.name in negative_atoms and atom2.name in positive_atoms)
        )

    def is_vdw(self, atom1, atom2, dist):
        return 1.6 < dist <= 4.0 and not (self.is_hbond(atom1, atom2, dist) or self.is_hydrophobic(atom1, atom2, dist))

    def is_aromatic(self, atom1, atom2, dist):
        aromatic_atoms = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "CH"}
        return dist <= 5.0 and atom1.name in aromatic_atoms and atom2.name in aromatic_atoms

    def is_cation_pi(self, atom1, atom2, dist):
        cation_atoms = {"NZ", "NH1", "NH2", "NE"}
        pi_atoms = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "CH"}
        return dist <= 6.0 and (
            (atom1.name in cation_atoms and atom2.name in pi_atoms) or
            (atom2.name in cation_atoms and atom1.name in pi_atoms)
        )

    def is_covalent(self, dist):
        return dist < 1.6

    def calculate_angle(self, atom1, atom2, atom3):
        vector1 = [
            atom1.coord[0] - atom2.coord[0],
            atom1.coord[1] - atom2.coord[1],
            atom1.coord[2] - atom2.coord[2]
        ]
        vector2 = [
            atom3.coord[0] - atom2.coord[0],
            atom3.coord[1] - atom2.coord[1],
            atom3.coord[2] - atom2.coord[2]
        ]
        dot_product = sum(v1 * v2 for v1, v2 in zip(vector1, vector2))
        magnitude1 = math.sqrt(sum(v**2 for v in vector1))
        magnitude2 = math.sqrt(sum(v**2 for v in vector2))
        if magnitude1 == 0 or magnitude2 == 0:
            return 0
        cosine_angle = dot_product / (magnitude1 * magnitude2)
        return math.degrees(math.acos(cosine_angle))

    def find_hydrogen_near(self, atom, max_dist=1.2):
        for neighbor in cmd.get_model(f"({atom.segi}//{atom.chain}/{atom.resi}) and name H*").atom:
            dist = cmd.get_distance(f"{atom.segi}//{atom.chain}/{atom.resi}/{atom.name}",
                                    f"{neighbor.segi}//{neighbor.chain}/{neighbor.resi}/{neighbor.name}")
            if dist < max_dist:
                return neighbor
        return None