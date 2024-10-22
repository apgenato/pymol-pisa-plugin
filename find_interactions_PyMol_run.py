import os
from pymol import cmd
import math

# Path to the directory where the file will be saved
output_directory = r"C:\path\to\output"
output_filename = "interactions_list.txt"
output_path = os.path.join(output_directory, output_filename)

# Names of receptor and ligand objects (don't change if ClusPro model)
rec_file = "rec"

lig_file = ""
for obj in cmd.get_object_list():
    if "lig" in obj.lower():
        lig_file = obj
        break

# Maximum distance cutoff for interaction analysis
distance_cutoff = 4.0

# Functions to calculate angles
def calculate_angle(atom1, atom2, atom3):
    """Calculate the angle (in degrees) formed by three atoms (atom1-atom2-atom3)."""
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

def find_hydrogen_near(atom, max_dist=1.2):
    """Find the hydrogen atom bonded to the given donor atom."""
    for neighbor in cmd.get_model(f"({atom.segi}//{atom.chain}/{atom.resi}) and name H*").atom:
        dist = cmd.get_distance(f"{atom.segi}//{atom.chain}/{atom.resi}/{atom.name}",
                                f"{neighbor.segi}//{neighbor.chain}/{neighbor.resi}/{neighbor.name}")
        if dist < max_dist:
            return neighbor
    return None

# Functions to determine various types of interactions
def is_hydrophobic(atom1, atom2, dist):
    # Hydrophobic interactions between carbons and sulfur at distances up to 4 Å
    hydrophobic_atoms = {"C", "S"}
    return dist <= 4.0 and atom1.name[0] in hydrophobic_atoms and atom2.name[0] in hydrophobic_atoms

def is_hbond(atom1, atom2, dist):
    # Hydrogen bonds consider distances up to 3.5 Å and an appropriate angle (>120°)
    if dist > 3.5:
        return False
    hbond_donors_acceptors = {"N", "O", "F"}
    if atom1.name[0] in hbond_donors_acceptors and atom2.name[0] in hbond_donors_acceptors:
        # Find a hydrogen atom attached to the donor (atom1)
        donor_hydrogen = find_hydrogen_near(atom1)
        if donor_hydrogen:
            angle = calculate_angle(atom1, donor_hydrogen, atom2)
            return angle > 120  # Only consider it a hydrogen bond if the angle is greater than 120°
    return False

def is_salt_bridge(atom1, atom2, dist):
    # Salt bridges between positively and negatively charged groups
    positive_atoms = {"NZ", "NH1", "NH2", "NE", "ND1", "ND2"}
    negative_atoms = {"OD1", "OD2", "OE1", "OE2"}
    return dist <= 4.0 and (
        (atom1.name in positive_atoms and atom2.name in negative_atoms) or
        (atom1.name in negative_atoms and atom2.name in positive_atoms)
    )

def is_vdw(atom1, atom2, dist):
    # Van der Waals interactions, if they are neither hydrophobic nor hydrogen bonds
    return 1.6 < dist <= 4.0 and not (is_hbond(atom1, atom2, dist) or is_hydrophobic(atom1, atom2, dist))

def is_aromatic(atom1, atom2, dist):
    # Aromatic interactions
    aromatic_atoms = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "CH"}
    return dist <= 5.0 and atom1.name in aromatic_atoms and atom2.name in aromatic_atoms

def is_cation_pi(atom1, atom2, dist):
    # Cation-π interactions
    cation_atoms = {"NZ", "NH1", "NH2", "NE"}
    pi_atoms = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "CH"}
    return dist <= 6.0 and (
        (atom1.name in cation_atoms and atom2.name in pi_atoms) or
        (atom2.name in cation_atoms and atom1.name in pi_atoms)
    )

def is_covalent(dist):
    # Covalent bonds, if the distance is less than 1.6 Å
    return dist < 1.6

# Collect all interactions
interactions = {
    "hydrophobic": [], "hbond": [], "salt_bridge": [], "covalent": [],
    "vdw": [], "aromatic": [], "cation_pi": [], "other": []
}

cmd.feedback("push")
cmd.feedback("disable", "all", "actions")
cmd.feedback("enable", "cmd", "results")

print(f"Running interaction analysis between {rec_file} and {lig_file}...")

# Select interacting atoms from the receptor and ligand
cmd.select("close_contacts", f"({rec_file} within {distance_cutoff} of {lig_file}) or ({lig_file} within {distance_cutoff} of {rec_file})")
cmd.select("interacting_residues", "byres close_contacts")

rec_atoms = cmd.get_model(f"{rec_file} and close_contacts").atom
lig_atoms = cmd.get_model(f"{lig_file} and close_contacts").atom

# Determine interactions
for rec_atom in rec_atoms:
    for lig_atom in lig_atoms:
        dist = cmd.get_distance(f"{rec_file}//{rec_atom.chain}/{rec_atom.resi}/{rec_atom.name}",
                                f"{lig_file}//{lig_atom.chain}/{lig_atom.resi}/{lig_atom.name}")
        if dist <= distance_cutoff:  # Ensure only interactions up to 4 Å are considered
            # Pass the actual atom objects instead of names
            if is_hbond(rec_atom, lig_atom, dist):
                interactions["hbond"].append((rec_atom, lig_atom, dist))
            elif is_salt_bridge(rec_atom, lig_atom, dist):
                interactions["salt_bridge"].append((rec_atom, lig_atom, dist))
            elif is_hydrophobic(rec_atom, lig_atom, dist):
                interactions["hydrophobic"].append((rec_atom, lig_atom, dist))
            elif is_aromatic(rec_atom, lig_atom, dist):
                interactions["aromatic"].append((rec_atom, lig_atom, dist))
            elif is_cation_pi(rec_atom, lig_atom, dist):
                interactions["cation_pi"].append((rec_atom, lig_atom, dist))
            elif is_covalent(dist):
                interactions["covalent"].append((rec_atom, lig_atom, dist))
            elif is_vdw(rec_atom, lig_atom, dist):
                interactions["vdw"].append((rec_atom, lig_atom, dist))
            else:
                interactions["other"].append((rec_atom, lig_atom, dist))

# Write results to a file
with open(output_path, "w") as f:
    for interaction_type, pairs in interactions.items():
        f.write(f"{interaction_type.upper()} INTERACTIONS:\n")
        for atom1, atom2, dist in pairs:
            f.write(f"{atom1.chain}/{atom1.resn}`{atom1.resi}/{atom1.name} -- {atom2.chain}/{atom2.resn}`{atom2.resi}/{atom2.name} : {dist:.2f} Å\n")
        f.write("\n")

print(f"Results saved to file: {output_path}")

# Visualization
cmd.show("sticks", "interacting_residues")

for interaction_type, color in zip(["hbond", "salt_bridge", "covalent", "hydrophobic"],
                                   ["cyan", "magenta", "yellow", "orange"]):
    if interactions[interaction_type]:
        for atom1, atom2, dist in interactions[interaction_type]:
            cmd.distance(f"{interaction_type}_bonds",
                         f"{rec_file}//{atom1.chain}/{atom1.resi}/{atom1.name}",
                         f"{lig_file}//{atom2.chain}/{atom2.resi}/{atom2.name}")
        cmd.show("dashes", f"{interaction_type}_bonds")
        cmd.color(color, f"{interaction_type}_bonds")

cmd.feedback("pop")
