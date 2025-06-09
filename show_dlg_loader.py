from pymol import cmd, Qt
from pymol.Qt import QtWidgets
import os


def __init_plugin__(app=None):
    """Initialize the PyMOL plugin."""
    from pymol.plugins import addmenuitemqt
    addmenuitemqt("Load DLG File", show_gui)


class DLGLoaderGUI(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle("DLG File Loader")
        layout = QtWidgets.QVBoxLayout()

        self.file_label = QtWidgets.QLabel("Select .dlg file:")
        layout.addWidget(self.file_label)

        self.file_input = QtWidgets.QLineEdit()
        layout.addWidget(self.file_input)

        self.browse_button = QtWidgets.QPushButton("Browse")
        self.browse_button.clicked.connect(self.browseFile)
        layout.addWidget(self.browse_button)

        self.load_button = QtWidgets.QPushButton("Load into PyMOL")
        self.load_button.clicked.connect(self.loadFile)
        layout.addWidget(self.load_button)

        self.setLayout(layout)

    def browseFile(self):
        options = QtWidgets.QFileDialog.Options()
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open DLG File", "",
                                                            "DLG Files (*.dlg);;All Files (*)", options=options)
        if filename:
            self.file_input.setText(filename)

    def loadFile(self):
        filename = self.file_input.text()
        if filename:
            load_dlg(filename)
        else:
            QtWidgets.QMessageBox.warning(self, "Error", "Please select a file first.")


def load_dlg(filename):
    """
    Load a .dlg file, extract docked models, and add them to PyMOL.

    Args:
        filename (str): Path to the .dlg file.
    """
    name = os.path.splitext(os.path.basename(filename))[0]
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Extract docked models
    docked_lines = [line[8:] for line in lines if line.startswith("DOCKED:")]
    models = []
    current_model = []

    for line in docked_lines:
        if line.startswith("MODEL"):
            if current_model:
                models.append(current_model)
            current_model = [line]
        elif line.startswith("ENDMDL"):
            current_model.append(line)
            models.append(current_model)
            current_model = []
        else:
            current_model.append(line)

    # Load extracted models into PyMOL
    for idx, model in enumerate(models):
        model_str = "".join(model)
        temp_file = f"{name}_pose{idx + 1}.pdb"
        with open(temp_file, 'w') as pdb_file:
            pdb_file.write(model_str)
        cmd.load(temp_file, f"{name}_pose{idx + 1}")
        os.remove(temp_file)  # Clean up temporary file

        # Hide everything and set representation to lines
        cmd.hide("everything", f"{name}_pose{idx + 1}")
        cmd.show("lines", f"{name}_pose{idx + 1}")

    print(f"Loaded {len(models)} docked poses from {filename}")

    # Group all loaded poses into a group named after the DLG file
    cmd.group(name, f"{name}_pose*")


# PyMOL command registration
cmd.extend("load_dlg", load_dlg)

# Drag-and-drop support is not available in this PyMOL version
# cmd.file_ext.append(("dlg", load_dlg))


def show_gui():
    global dlg_loader_gui
    dlg_loader_gui = DLGLoaderGUI()
    dlg_loader_gui.show()


cmd.extend("show_dlg_loader", show_gui)