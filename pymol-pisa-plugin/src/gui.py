from PyQt5 import QtWidgets, QtGui
import sys
from pymol import cmd
from .core import analyze_interactions

class PISAPluginGUI(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PyMOL PISA Plugin")
        self.setGeometry(100, 100, 400, 300)
        
        self.layout = QtWidgets.QVBoxLayout()
        
        self.rec_label = QtWidgets.QLabel("Receptor Object Name:")
        self.layout.addWidget(self.rec_label)
        
        self.rec_input = QtWidgets.QLineEdit()
        self.layout.addWidget(self.rec_input)
        
        self.lig_label = QtWidgets.QLabel("Ligand Object Name:")
        self.layout.addWidget(self.lig_label)
        
        self.lig_input = QtWidgets.QLineEdit()
        self.layout.addWidget(self.lig_input)
        
        self.analyze_button = QtWidgets.QPushButton("Analyze Interactions")
        self.analyze_button.clicked.connect(self.analyze_interactions)
        self.layout.addWidget(self.analyze_button)
        
        self.output_label = QtWidgets.QLabel("")
        self.layout.addWidget(self.output_label)
        
        self.setLayout(self.layout)

    def analyze_interactions(self):
        rec_name = self.rec_input.text()
        lig_name = self.lig_input.text()
        
        if rec_name and lig_name:
            output_path = analyze_interactions(rec_name, lig_name)
            self.output_label.setText(f"Results saved to: {output_path}")
        else:
            self.output_label.setText("Please enter both receptor and ligand names.")

def main():
    app = QtWidgets.QApplication(sys.argv)
    window = PISAPluginGUI()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()