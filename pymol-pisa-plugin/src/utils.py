def calculate_distance(atom1, atom2):
    """Calculate the distance between two atoms."""
    return math.sqrt(
        (atom1.coord[0] - atom2.coord[0]) ** 2 +
        (atom1.coord[1] - atom2.coord[1]) ** 2 +
        (atom1.coord[2] - atom2.coord[2]) ** 2
    )

def write_interaction_results(output_path, interactions):
    """Write interaction results to a specified output file."""
    with open(output_path, "w") as f:
        for interaction_type, pairs in interactions.items():
            f.write(f"{interaction_type.upper()} INTERACTIONS:\n")
            for atom1, atom2, dist in pairs:
                f.write(f"{atom1.chain}/{atom1.resn}`{atom1.resi}/{atom1.name} -- {atom2.chain}/{atom2.resn}`{atom2.resi}/{atom2.name} : {dist:.2f} Å\n")
            f.write("\n")

def get_ligand_object(cmd):
    """Retrieve the ligand object from the PyMOL session."""
    for obj in cmd.get_object_list():
        if "lig" in obj.lower():
            return obj
    return None

def select_interacting_residues(rec_file, lig_file, distance_cutoff):
    """Select residues that are within a specified distance cutoff."""
    cmd.select("close_contacts", f"({rec_file} within {distance_cutoff} of {lig_file}) or ({lig_file} within {distance_cutoff} of {rec_file})")
    cmd.select("interacting_residues", "byres close_contacts")