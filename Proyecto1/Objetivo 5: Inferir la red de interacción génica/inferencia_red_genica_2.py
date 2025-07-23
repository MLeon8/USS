import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

# Crear red de interacción
genes = ['TP53', 'MDM2', 'ATM', 'CHEK2', 'BRCA1', 'CDKN1A', 'RB1', 'E2F1', 'CDK2', 'EGFR', 'MYC', 'CDK4']
edges = [
    ('TP53', 'MDM2'), ('TP53', 'ATM'), ('ATM', 'CHEK2'), ('CHEK2', 'BRCA1'),
    ('BRCA1', 'CDKN1A'), ('CDKN1A', 'RB1'), ('RB1', 'E2F1'), ('E2F1', 'CDK2'),
    ('EGFR', 'MYC'), ('MYC', 'CDK4')
]
G = nx.Graph()
G.add_edges_from(edges)

# Métricas
print("Coeficiente de clustering promedio:", nx.average_clustering(G))
print("Grado promedio:", sum(dict(G.degree()).values()) / G.number_of_nodes())
print("Camino promedio:", nx.average_shortest_path_length(G))
print("Diámetro:", nx.diameter(G))
print("Conectada:", nx.is_connected(G))

# Geodésica
path = nx.shortest_path(G, source="TP53", target="CDK2")
print("Geodésica TP53 - CDK2:", path)
