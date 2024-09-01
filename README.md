# find_interactions_PyMol_run.py

## Description

The `find_interactions_PyMol_run.py` script is designed to analyze molecular interactions between a receptor and a ligand within a structure loaded in PyMOL. It automatically identifies and classifies various types of interactions, such as hydrogen bonds, hydrophobic interactions, salt bridges, van der Waals interactions, and others.

This script is useful for post-docking analysis in molecular modeling, enabling quick identification and visualization of key interactions between molecules.

## Features

- **Hydrophobic Interactions:** Identified by the presence of carbon and sulfur atoms within a 4 Å distance.
- **Hydrogen Bonds:** Consider donors and acceptors of hydrogen bonds (N, O, F) within a 3.5 Å distance.
- **Salt Bridges:** Interactions between positively and negatively charged atoms within a 4 Å distance.
- **Van der Waals Interactions:** Determined based on distance (1.6-4.0 Å) if not classified as hydrophobic or hydrogen bonds.
- **Aromatic and Cation-π Interactions:** Analyzed based on specific atoms and distances.
- **Visualization:** Automatically displays interacting residues and draws bonds between them in PyMOL.

## Usage

1. Save the `find_interactions_PyMol_run.py` script in the desired directory.
2. Open PyMOL and load your receptor and ligand structure.
3. Run the script in the PyMOL command line:
   ```python
   run C:/path/to/find_interactions_PyMol_run.py
