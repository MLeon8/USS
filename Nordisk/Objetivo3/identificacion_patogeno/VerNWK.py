from Bio import Phylo
import matplotlib.pyplot as plt

tree = Phylo.read("arbol_filogenetico.nwk", "newick")
Phylo.draw(tree)  # Dibuja en ventana gráfica
