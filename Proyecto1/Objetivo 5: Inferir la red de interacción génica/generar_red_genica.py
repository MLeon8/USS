import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

# Genes y conexiones hipotéticas
genes = ['TP53', 'MDM2', 'ATM', 'CHEK2', 'BRCA1', 'CDKN1A', 'RB1', 'E2F1', 'CDK2']
edges = [
    ('TP53', 'MDM2'), ('TP53', 'ATM'), ('ATM', 'CHEK2'), ('CHEK2', 'BRCA1'),
    ('BRCA1', 'CDKN1A'), ('CDKN1A', 'RB1'), ('RB1', 'E2F1'), ('E2F1', 'CDK2')
]

# Crear el grafo
G = nx.Graph()
G.add_edges_from(edges)

# Layouts y visualización
spring_pos = nx.spring_layout(G, seed=42)
kamada_pos = nx.kamada_kawai_layout(G)

plt.figure(figsize=(8,6))
nx.draw(G, pos=spring_pos, with_labels=True, node_color='skyblue', edge_color='gray')
plt.title("Spring Layout")
plt.savefig("red_spring_layout.png", dpi=300)
plt.clf()

nx.draw(G, pos=kamada_pos, with_labels=True, node_color='lightgreen', edge_color='gray')
plt.title("Kamada-Kawai Layout")
plt.savefig("red_kamada_layout.png", dpi=300)

# Métricas
print("Clustering:", nx.average_clustering(G))
print("Grado promedio:", sum(dict(G.degree()).values()) / G.number_of_nodes())
print("Camino largo promedio:", nx.average_shortest_path_length(G))
print("Diámetro:", nx.diameter(G))
print("Geodésica (TP53-CDK2):", nx.shortest_path(G, 'TP53', 'CDK2'))
