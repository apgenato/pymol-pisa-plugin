import os
import math
import csv
from pymol import cmd, Qt

MAX_VISUALIZED_INTERACTIONS = 2000

# --- Qt6/PySide6 compatibility patch ---
try:
    from pymol.Qt import QtWidgets
    # Try to access Qt6 enums to check version
    _ = QtWidgets.QMessageBox.StandardButton.Ok
    QT6 = True
except AttributeError:
    from pymol.Qt import QtWidgets
    QT6 = False
except ImportError:
    from PyQt6 import QtWidgets
    QT6 = True

METAL_ELEMENTS = {
    "LI", "NA", "K", "RB", "CS",
    "MG", "CA", "CAL", "SR", "BA",
    "AL", "GA", "IN", "TL",
    "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG",
    "MO", "W", "RU", "RH", "PD", "AG", "OS", "IR", "PT", "AU"
}

METAL_PARTNER_ELEMENTS = {"N", "O", "S", "F", "CL", "BR", "I"}

def calculate_angle(atom1, atom2, atom3):
    """
    Calculates the angle formed by three atoms in 3D space.

    Parameters:
        atom1 (object): The first atom with a 'coord' attribute containing its coordinates.
        atom2 (object): The second atom (vertex of the angle) with a 'coord' attribute.
        atom3 (object): The third atom with a 'coord' attribute containing its coordinates.

    Returns:
        float: The angle in degrees between the vectors formed by atom1-atom2 and atom3-atom2.
               Returns 0 if either vector has zero magnitude.
    """
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
    # Clamp value to avoid math domain error
    cosine_angle = max(-1.0, min(1.0, cosine_angle))
    return math.degrees(math.acos(cosine_angle))

def find_nearest_hydrogen(atom):
    """
    Finds a hydrogen atom near the specified atom within a given distance.

    Parameters:
        atom (object): The atom to search around.
        max_dist (float): Maximum distance to search for hydrogen atoms (default is 1.2 Å).

    Returns:
        object: The nearest hydrogen atom within the specified distance, or None if no hydrogen is found.
    """
    # Build selection string for the parent object and residue
    obj = atom.model
    sele = f'{obj} and chain {atom.chain} and resi {atom.resi} and elem H'
    hydrogens = cmd.get_model(sele).atom
    closest_neighbor = None
    closest_distance = float('inf')
    for neighbor in hydrogens:
        # Skip if it's the same atom
        if neighbor.name == atom.name:
            continue
        # Calculate distance using coordinates directly
        dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(atom.coord, neighbor.coord)))
        if dist < 1.2 and dist < closest_distance:
            closest_neighbor = neighbor
            closest_distance = dist
    return closest_neighbor

def is_hydrophobic(atom1, atom2, dist):
    hydrophobic_atoms = {"C", "S"}
    return dist <= 4.0 and atom1.name[0] in hydrophobic_atoms and atom2.name[0] in hydrophobic_atoms

def atom_element(atom):
    """Return an uppercase element symbol, with a conservative fallback."""
    symbol = getattr(atom, "symbol", "")
    if symbol:
        return symbol.upper()
    return atom.name[:2].strip().upper()

def is_metal_contact(atom1, atom2, dist):
    """Classify metal coordination-like contacts (metal-ligand proximity)."""
    if dist > 3.0:
        return False
    elem1 = atom_element(atom1)
    elem2 = atom_element(atom2)
    return (
        (elem1 in METAL_ELEMENTS and elem2 in METAL_PARTNER_ELEMENTS) or
        (elem2 in METAL_ELEMENTS and elem1 in METAL_PARTNER_ELEMENTS)
    )

def is_hbond(atom1, atom2, dist):
    hbond_donors_acceptors = {"N", "O", "F"}
    if dist < 3.5 and atom1.name[0] in hbond_donors_acceptors and atom2.name[0] in hbond_donors_acceptors:
        donor_hydrogen1 = find_nearest_hydrogen(atom1)
        donor_hydrogen2 = find_nearest_hydrogen(atom2)
        if donor_hydrogen1:
            angle1 = calculate_angle(atom1, donor_hydrogen1, atom2)
            if angle1 > 120:
                return True
            return False
        if donor_hydrogen2:
            angle2 = calculate_angle(atom2, donor_hydrogen2, atom1)
            if angle2 > 120:
                return True
            return False
        return False

def is_salt_bridge(atom1, atom2, dist):
    positive_atoms = {"NZ", "NH1", "NH2", "NE", "ND1", "ND2"}
    negative_atoms = {"OD1", "OD2", "OE1", "OE2", "OE11", "OE12", "OE21", "OE22"}
    return dist <= 4.0 and (
        (atom1.name in positive_atoms and atom2.name in negative_atoms) or
        (atom1.name in negative_atoms and atom2.name in positive_atoms)
    )

def is_vdw(atom1, atom2, dist):
    return 1.6 < dist <= 4.0 and not (is_hbond(atom1, atom2, dist) or is_hydrophobic(atom1, atom2, dist))

def is_aromatic(atom1, atom2, dist):
    aromatic_atoms = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "CH"}
    return dist <= 5.0 and atom1.name in aromatic_atoms and atom2.name in aromatic_atoms

def is_cation_pi(atom1, atom2, dist):
    cation_atoms = {"NZ", "NH1", "NH2", "NE"}
    pi_atoms = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "CH"}
    return dist <= 6.0 and (
        (atom1.name in cation_atoms and atom2.name in pi_atoms) or
        (atom2.name in cation_atoms and atom1.name in pi_atoms)
    )

def is_covalent(dist):
    return dist < 1.6

def atom_distance(atom1, atom2):
    """Compute Euclidean distance directly from atom coordinates."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(atom1.coord, atom2.coord)))

def atom_selection(obj_name, atom):
    """Build a robust atom selection using PyMOL's internal atom index."""
    return f"{obj_name} and index {atom.index}"

def run_analysis(rec_file, lig_file, output_path, distance_cutoff=4.0, ensure_hydrogens=True):
    # Optionally check for hydrogens in receptor and ligand, add if missing
    if ensure_hydrogens:
        rec_has_h = cmd.count_atoms(f"{rec_file} and elem H") > 0
        lig_has_h = cmd.count_atoms(f"{lig_file} and elem H") > 0
        if not rec_has_h:
            cmd.h_add(rec_file)
        if not lig_has_h:
            cmd.h_add(lig_file)

    interactions = {
        "hydrophobic": [], "hbond": [], "metal_contact": [], "salt_bridge": [], "covalent": [],
        "vdw": [], "aromatic": [], "cation_pi": [], "other": []
    }

    cmd.feedback("push")
    cmd.feedback("disable", "all", "actions")
    cmd.feedback("enable", "cmd", "results")

    cmd.select("close_contacts", f"({rec_file} within {distance_cutoff} of {lig_file}) or ({lig_file} within {distance_cutoff} of {rec_file})")
    cmd.select("interacting_residues", "byres close_contacts")

    rec_atoms = cmd.get_model(f"{rec_file} and close_contacts").atom
    lig_atoms = cmd.get_model(f"{lig_file} and close_contacts").atom

    if rec_file == lig_file:
        for i, rec_atom in enumerate(rec_atoms):
            for j, lig_atom in enumerate(lig_atoms):
                if j <= i:
                    continue  # skip self and duplicate pairs
                dist = atom_distance(rec_atom, lig_atom)
                if dist <= distance_cutoff:
                    if rec_file == lig_file and rec_atom.index == lig_atom.index:
                        continue  # skip self-self
                    if is_hbond(rec_atom, lig_atom, dist):
                        interactions["hbond"].append((rec_atom, lig_atom, dist))
                    elif is_metal_contact(rec_atom, lig_atom, dist):
                        interactions["metal_contact"].append((rec_atom, lig_atom, dist))
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
    else:
        for rec_atom in rec_atoms:
            for lig_atom in lig_atoms:
                dist = atom_distance(rec_atom, lig_atom)
                if dist <= distance_cutoff:
                    if rec_file == lig_file and rec_atom.index == lig_atom.index:
                        continue  # skip self-self
                    if is_hbond(rec_atom, lig_atom, dist):
                        interactions["hbond"].append((rec_atom, lig_atom, dist))
                    elif is_metal_contact(rec_atom, lig_atom, dist):
                        interactions["metal_contact"].append((rec_atom, lig_atom, dist))
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

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "interaction_type",
            "atom1_chain", "atom1_resn", "atom1_resi", "atom1_name", "atom1_element","--",
            "atom2_chain", "atom2_resn", "atom2_resi", "atom2_name", "atom2_element",
            "distance (Å)"
        ])
        for interaction_type, pairs in interactions.items():
            for atom1, atom2, dist in pairs:
                writer.writerow([
                    interaction_type,
                    atom1.chain, atom1.resn, atom1.resi, atom1.name, atom1.symbol,"--",
                    atom2.chain, atom2.resn, atom2.resi, atom2.name, atom2.symbol,
                    f"{dist:.2f}"
                ])

    cmd.show("sticks", "interacting_residues")
    for interaction_type, color in zip(["hbond", "metal_contact", "salt_bridge", "covalent", "hydrophobic"],
                                       ["cyan", "green", "magenta", "yellow", "orange"]):
        # Clear previous measurements from earlier runs to avoid accumulation.
        cmd.delete(f"{interaction_type}_bonds")
        if interactions[interaction_type]:
            for atom1, atom2, dist in interactions[interaction_type][:MAX_VISUALIZED_INTERACTIONS]:
                cmd.distance(f"{interaction_type}_bonds",
                             atom_selection(rec_file, atom1),
                             atom_selection(lig_file, atom2))
            cmd.show("dashes", f"{interaction_type}_bonds")
            cmd.color(color, f"{interaction_type}_bonds")
    cmd.feedback("pop")
    return output_path

class PisaPluginGUI(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle("PISA Interaction Analysis")
        layout = QtWidgets.QGridLayout()

        # Get object list from PyMOL
        self.objects = cmd.get_object_list()

        # Receptor
        layout.addWidget(QtWidgets.QLabel("Receptor object:"), 0, 0)
        self.rec_combo = QtWidgets.QComboBox()
        self.rec_combo.addItems(self.objects)
        layout.addWidget(self.rec_combo, 0, 1)

        # Ligand
        layout.addWidget(QtWidgets.QLabel("Ligand object:"), 1, 0)
        self.lig_combo = QtWidgets.QComboBox()
        self.lig_combo.addItems(self.objects)
        layout.addWidget(self.lig_combo, 1, 1)

        # Output directory
        layout.addWidget(QtWidgets.QLabel("Output directory:"), 2, 0)
        self.outdir_edit = QtWidgets.QLineEdit(os.path.expanduser("~"))
        layout.addWidget(self.outdir_edit, 2, 1)
        self.browse_btn = QtWidgets.QPushButton("Browse")
        self.browse_btn.clicked.connect(self.browse_dir)
        layout.addWidget(self.browse_btn, 2, 2)

        # Output filename
        layout.addWidget(QtWidgets.QLabel("Output filename:"), 3, 0)
        default_filename = f"{self.rec_combo.currentText()}_{self.lig_combo.currentText()}_interactions.csv"
        self.outname_edit = QtWidgets.QLineEdit(default_filename)
        self.rec_combo.currentIndexChanged.connect(self.update_filename)
        self.lig_combo.currentIndexChanged.connect(self.update_filename)
        layout.addWidget(self.outname_edit, 3, 1)

        # Add hydrogens checkbox
        self.hydrogen_checkbox = QtWidgets.QCheckBox("Ensure hydrogens are present (add if missing)")
        self.hydrogen_checkbox.setChecked(True)
        layout.addWidget(self.hydrogen_checkbox, 4, 0, 1, 2)

        # Allow intra-object analysis checkbox
        self.intra_checkbox = QtWidgets.QCheckBox("Allow intra-object (self) analysis")
        self.intra_checkbox.setChecked(False)
        layout.addWidget(self.intra_checkbox, 5, 0, 1, 2)

        # Run button
        self.run_btn = QtWidgets.QPushButton("Run Analysis")
        self.run_btn.clicked.connect(self.run_analysis_gui)
        layout.addWidget(self.run_btn, 6, 0, 1, 3)

        self.setLayout(layout)

    def update_filename(self):
        rec = self.rec_combo.currentText()
        lig = self.lig_combo.currentText()
        self.outname_edit.setText(f"{rec}_{lig}_interactions.csv")

    def browse_dir(self):
        d = QtWidgets.QFileDialog.getExistingDirectory(self, "Select Output Directory", self.outdir_edit.text())
        # Qt6 returns a tuple (dir, selectedFilter), Qt5 returns str
        if isinstance(d, tuple):
            d = d[0]
        if d:
            self.outdir_edit.setText(d)

    def run_analysis_gui(self):
        rec = self.rec_combo.currentText()
        lig = self.lig_combo.currentText()
        outdir = self.outdir_edit.text()
        outname = self.outname_edit.text()
        ensure_hydrogens = self.hydrogen_checkbox.isChecked()
        allow_intra = self.intra_checkbox.isChecked()
        if not rec or not lig:
            if QT6:
                QtWidgets.QMessageBox.warning(self, "Error", "Please select both receptor and ligand objects.", QtWidgets.QMessageBox.StandardButton.Ok)
            else:
                QtWidgets.QMessageBox.warning(self, "Error", "Please select both receptor and ligand objects.")
            return
        if rec == lig and not allow_intra:
            if QT6:
                QtWidgets.QMessageBox.warning(self, "Error", "Receptor and ligand must be different objects (or enable intra-object analysis).", QtWidgets.QMessageBox.StandardButton.Ok)
            else:
                QtWidgets.QMessageBox.warning(self, "Error", "Receptor and ligand must be different objects (or enable intra-object analysis).")
            return
        if not os.path.isdir(outdir):
            if QT6:
                QtWidgets.QMessageBox.warning(self, "Error", "Output directory does not exist.", QtWidgets.QMessageBox.StandardButton.Ok)
            else:
                QtWidgets.QMessageBox.warning(self, "Error", "Output directory does not exist.")
            return
        output_path = os.path.join(outdir, outname)
        try:
            run_analysis(rec, lig, output_path, ensure_hydrogens=ensure_hydrogens)
            if QT6:
                QtWidgets.QMessageBox.information(self, "Done", f"Results saved to {output_path}", QtWidgets.QMessageBox.StandardButton.Ok)
            else:
                QtWidgets.QMessageBox.information(self, "Done", f"Results saved to {output_path}")
        except Exception as e:
            if QT6:
                QtWidgets.QMessageBox.critical(self, "Error", str(e), QtWidgets.QMessageBox.StandardButton.Ok)
            else:
                QtWidgets.QMessageBox.critical(self, "Error", str(e))

def show_gui():
    global pisa_gui
    pisa_gui = PisaPluginGUI()
    pisa_gui.show()

def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt("PISA Interaction Analysis", show_gui)
