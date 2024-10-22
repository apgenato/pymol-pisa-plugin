import os
from pymol import cmd

# Path to the directory where the file will be saved
output_directory = r"C:\path\to\output"
output_filename = "interactions_list.txt"
output_path = os.path.join(output_directory, output_filename)

rec_file = "rec"

lig_file = ""
for obj in cmd.get_object_list():
    if "lig" in obj.lower():
        lig_file = obj
        break


distance_cutoff = 4.0

# Functions to determine all possible types of interactions
def is_hydrophobic(atom1, atom2, dist):
    # Hydrophobic interactions between carbons and sulfur at distances up to 4 Å
    hydrophobic_atoms = {"C", "S"}
    return dist <= 4.0 and atom1[0] in hydrophobic_atoms and atom2[0] in hydrophobic_atoms

def is_hbond(atom1, atom2, dist):
    # Hydrogen bonds consider distances up to 3.5 Å and the presence of suitable donors/acceptors.
    if dist > 3.5:
        return False
    hbond_donors = {"N", "O", "F"} 
    return (atom1 in hbond_donors and atom2 in hbond_donors)

def is_salt_bridge(atom1, atom2, dist):
    positive_atoms = {"NZ", "NH1", "NH2", "NE", "ND1", "ND2"}
    negative_atoms = {"OD1", "OD2", "OE1", "OE2"}
    return dist <= 4.0 and (
        (atom1 in positive_atoms and atom2 in negative_atoms) or
        (atom1 in negative_atoms and atom2 in positive_atoms)
    )

def is_vdw(atom1, atom2, dist):
    # Van der Waals interactions are defined if they are neither hydrophobic nor hydrogen bonds
    return 1.6 < dist <= 4.0 and not (is_hbond(atom1, atom2, dist) or is_hydrophobic(atom1, atom2, dist))

def is_aromatic(atom1, atom2, dist):
    aromatic_atoms = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "CH"}
    return dist <= 5.0 and atom1 in aromatic_atoms and atom2 in aromatic_atoms

def is_cation_pi(atom1, atom2, dist):
    cation_atoms = {"NZ", "NH1", "NH2", "NE"}
    pi_atoms = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "CH"}
    return dist <= 6.0 and (
        (atom1 in cation_atoms and atom2 in pi_atoms) or
        (atom2 in cation_atoms and atom1 in pi_atoms)
    )

def is_covalent(dist):
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

# Iterate over all atoms in lig and rec to find interactions
cmd.select("close_contacts", f"({rec_file} within {distance_cutoff} of {lig_file}) or ({lig_file} within {distance_cutoff} of {rec_file})")
cmd.select("interacting_residues", "byres close_contacts")

rec_atoms = cmd.get_model(f"{rec_file} and close_contacts").atom
lig_atoms = cmd.get_model(f"{lig_file} and close_contacts").atom

# Collect all interactions in one pass
for rec_atom in rec_atoms:
    rec_name = rec_atom.name
    for lig_atom in lig_atoms:
        lig_name = lig_atom.name
        dist = cmd.get_distance(f"{rec_file}//{rec_atom.chain}/{rec_atom.resi}/{rec_atom.name}",
                                f"{lig_file}//{lig_atom.chain}/{lig_atom.resi}/{lig_atom.name}")
        if dist <= distance_cutoff:  # Ensure only interactions up to 4 Å are considered
            if is_hbond(rec_name, lig_name, dist):
                interactions["hbond"].append((rec_atom, lig_atom, dist))
            elif is_salt_bridge(rec_name, lig_name, dist):
                interactions["salt_bridge"].append((rec_atom, lig_atom, dist))
            elif is_hydrophobic(rec_name, lig_name, dist):
                interactions["hydrophobic"].append((rec_atom, lig_atom, dist))
            elif is_aromatic(rec_name, lig_name, dist):
                interactions["aromatic"].append((rec_atom, lig_atom, dist))
            elif is_cation_pi(rec_name, lig_name, dist):
                interactions["cation_pi"].append((rec_atom, lig_atom, dist))
            elif is_covalent(dist):
                interactions["covalent"].append((rec_atom, lig_atom, dist))
            elif is_vdw(rec_name, lig_name, dist):
                interactions["vdw"].append((rec_atom, lig_atom, dist))
            else:
                interactions["other"].append((rec_atom, lig_atom, dist))

# Output results
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
