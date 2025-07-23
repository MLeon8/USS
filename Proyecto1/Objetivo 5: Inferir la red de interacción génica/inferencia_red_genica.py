#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import requests
from collections import defaultdict

# Cargar genes sobreexpresados
df = pd.read_csv("GO_annotated_overexpressed_genes.csv")
genes = df["query"].unique().tolist()[:50]  # Limitar a 50 genes para eficiencia

# Descargar interacciones desde STRING
def get_string_interactions(gene_list, species=9606, score_cutoff=700):
    interactions = []
    for gene in gene_list:
        url = f"https://string-db.org/api/tsv/network?identifiers={gene}&species={species}&required_score={score_cutoff}"
        try:
            r = requests.get(url)
            if r.ok:
                lines = r.text.strip().split("\n")[1:]
                for line in lines:
                    data = line.split("\t")
                    if len(data) > 5:
                        interactions.append((data[2], data[3], float(data[5])))
        except:
            continue
    return interactions

edges = get_string_interactions(genes)

# Construir red
G = nx.Graph()
for a, b, score in edges:
    G.add_edge(a, b, weight=score)

print(f"Nodos: {G.number_of_nodes()}  Aristas: {G.number_of_edges()}")

# Graficar red con spring_layout
plt.figure(figsize=(10, 8))
pos = nx.spring_layout(G, seed=42)
nx.draw_networkx_nodes(G, pos, node_size=200, node_color="skyblue")
nx.draw_networkx_edges(G, pos, alpha=0.5)
nx.draw_networkx_labels(G, pos, font_size=8)
plt.title("Red de interacción génica (Spring Layout)")
plt.axis("off")
plt.tight_layout()
plt.savefig("red_spring_layout.png")
plt.show()

# Graficar red con kamada_kawai_layout
plt.figure(figsize=(10, 8))
pos2 = nx.kamada_kawai_layout(G)
nx.draw(G, pos2, with_labels=True, node_size=200, node_color="salmon", edge_color="gray")
plt.title("Red de interacción génica (Kamada-Kawai Layout)")
plt.axis("off")
plt.tight_layout()
plt.savefig("red_kamada_layout.png")
plt.show()

# --- Análisis topológico
print("\n--- Análisis Topológico ---")
print("Coef. Clustering:", nx.average_clustering(G))
print("Grado promedio:", sum(dict(G.degree()).values()) / G.number_of_nodes())
print("Camino largo promedio:", nx.average_shortest_path_length(G))
print("Diámetro:", nx.diameter(G))

# Geodésica (camino más largo mínimo)
u, v = nx.periphery(G)[0:2]
geo = nx.shortest_path(G, u, v)
print(f"Geodésica ({u} → {v}): {geo}")
print("Nodos geodésicos anotados:")

# Anotar nodos de la geodésica
geo_df = pd.DataFrame(geo, columns=["gene"])
geo_annotated = geo_df.merge(df[["query", "name"]], left_on="gene", right_on="query", how="left")
print(geo_annotated[["gene", "name"]])

# Probabilidad de conectividad
connectivity = nx.is_connected(G)
print("¿Red conectada completamente?:", connectivity)
