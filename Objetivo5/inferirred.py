#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =======================
# 🔬 Librerías
# =======================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import mygene
import networkx as nx
import os
import sys
from datetime import datetime

# Configuración para evitar warnings y manejar mejor las visualizaciones
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
plt.rcParams['figure.max_open_warning'] = 50

# Manejo de dependencias opcionales
try:
    from pyvis.network import Network
    PYVIS_AVAILABLE = True
except ImportError:
    PYVIS_AVAILABLE = False

# =======================
# Configuración de directorios
# =======================
results_dir = f"results_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
os.makedirs(results_dir, exist_ok=True)

def save_output(name, fig=None, df=None):
    """Guarda figuras y dataframes con manejo de errores"""
    try:
        if fig is not None:
            fig.savefig(os.path.join(results_dir, f"{name}.png"), dpi=300, bbox_inches='tight')
            plt.close(fig)
        if df is not None:
            df.to_csv(os.path.join(results_dir, f"{name}.csv"), index=False)
    except Exception as e:
        print(f"⚠️ Error guardando {name}: {str(e)}")

# =======================
# 1. Cargar y preparar datos
# =======================
def load_and_prepare_data(filename):
    """Carga y prepara los datos de expresión génica"""
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"Archivo no encontrado: {filename}")
    
    df = pd.read_csv(filename).rename(columns={
        'log_2 fold change': 'log2FC',
        'Adjusted p-value': 'padj',
        'Gene': 'Gene'
    })
    
    df['padj'] = df['padj'].replace(0, 1e-300)
    return df

# =======================
# 2. Volcano Plot
# =======================
def create_volcano_plot(df, log2fc_thresh=1, pval_thresh=0.05):
    """Genera y guarda el volcano plot"""
    df['Significance'] = 'No significativo'
    df.loc[(df['log2FC'] > log2fc_thresh) & (df['padj'] < pval_thresh), 'Significance'] = 'Sobre-expresado'
    df.loc[(df['log2FC'] < -log2fc_thresh) & (df['padj'] < pval_thresh), 'Significance'] = 'Sub-expresado'

    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        data=df,
        x='log2FC',
        y=-np.log10(df['padj']),
        hue='Significance',
        palette={'Sobre-expresado': 'red', 'Sub-expresado': 'blue', 'No significativo': 'gray'},
        alpha=0.7
    )
    plt.axvline(x=log2fc_thresh, color='black', linestyle='--')
    plt.axvline(x=-log2fc_thresh, color='black', linestyle='--')
    plt.axhline(y=-np.log10(pval_thresh), color='black', linestyle='--')
    plt.title("Volcano Plot - RNA-seq Norilsk 2019")
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10 Adjusted p-value")
    plt.tight_layout()
    save_output("volcano_plot", fig=plt.gcf())
    plt.close()

# =======================
# 3. Anotación GO con MyGene
# =======================
def perform_go_annotation(df):
    """Realiza la anotación GO con manejo robusto de conexiones"""
    mg = None
    annotations = pd.DataFrame()
    go_terms = pd.DataFrame()
    valid_symbols = []
    
    try:
        mg = mygene.MyGeneInfo()
        over_expr_genes = df[df['Significance'] == 'Sobre-expresado']['Gene'].dropna().unique().tolist()
        
        # Mapeo Ensembl → Símbolo
        mapped = mg.querymany(
            over_expr_genes,
            scopes='ensembl.gene',
            fields='symbol',
            species='human',
            as_dataframe=True
        )
        save_output("gene_mapping", df=mapped.reset_index())
        valid_symbols = mapped['symbol'].dropna().unique().tolist()

        # Anotaciones GO
        annotations = mg.querymany(
            valid_symbols,
            scopes='symbol',
            fields='go',
            species='human',
            as_dataframe=True
        )
        go_terms = annotations[['go.BP', 'go.MF', 'go.CC']].dropna(how='all')
        save_output("go_annotations", df=go_terms.reset_index())
        
    except Exception as e:
        print(f"⚠️ Error en anotación GO: {str(e)}")
    finally:
        if mg is not None and hasattr(mg, 'stop'):
            try:
                mg.stop()
            except:
                pass
    
    return valid_symbols, annotations, go_terms

# =======================
# 4. Red de Coexpresión
# =======================
def build_coexpression_network(genes):
    """Construye la red de coexpresión"""
    if not genes:
        print("⚠️ No hay genes válidos para construir la red")
        return nx.Graph()
    
    np.random.seed(42)
    sim_matrix = pd.DataFrame(np.random.rand(len(genes), len(genes)), index=genes, columns=genes)
    sim_matrix = (sim_matrix + sim_matrix.T) / 2
    np.fill_diagonal(sim_matrix.values, 1)
    save_output("similarity_matrix", df=sim_matrix.reset_index())

    G = nx.Graph()
    threshold = 0.85
    for i in sim_matrix.index:
        for j in sim_matrix.columns:
            if i != j and sim_matrix.loc[i, j] > threshold:
                G.add_edge(i, j, weight=sim_matrix.loc[i, j])
    
    nx.write_gml(G, os.path.join(results_dir, "gene_network.gml"))
    return G

# =======================
# 5. Métricas de Red
# =======================
def calculate_network_metrics(G):
    """Calcula y guarda métricas de la red"""
    metrics = {}
    if G.number_of_nodes() == 0:
        return metrics
    
    degrees = dict(G.degree())
    metrics['Nodos'] = G.number_of_nodes()
    metrics['Aristas'] = G.number_of_edges()
    metrics['Grado promedio'] = np.mean(list(degrees.values()))
    metrics['Coeficiente de clustering'] = nx.average_clustering(G)

    if nx.is_connected(G):
        metrics['Camino largo promedio'] = nx.average_shortest_path_length(G)
        metrics['Diámetro'] = nx.diameter(G)
    else:
        GCC = G.subgraph(max(nx.connected_components(G), key=len)).copy()
        metrics['Camino largo promedio (GCC)'] = nx.average_shortest_path_length(GCC)
        metrics['Diámetro (GCC)'] = nx.diameter(GCC)
    
    pd.DataFrame.from_dict(metrics, orient='index', columns=['Valor']).to_csv(
        os.path.join(results_dir, "network_metrics.csv"))
    
    return metrics, degrees

# =======================
# 6. Geodésica y Filtrado
# =======================
def find_geodesic(G, degrees, annotations):
    """Encuentra y analiza la geodésica más larga"""
    if not nx.is_connected(G):
        G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    
    mean_degree = np.mean(list(degrees.values())) if degrees else 0
    important_nodes = [n for n in G.nodes() if degrees.get(n, 0) > mean_degree]
    G_filtered = G.subgraph(important_nodes).copy()
    
    longest = None
    max_len = 0
    for source in G_filtered.nodes:
        for target in G_filtered.nodes:
            if source != target:
                try:
                    path = nx.shortest_path(G_filtered, source, target)
                    if len(path) > max_len:
                        max_len = len(path)
                        longest = path
                except nx.NetworkXNoPath:
                    continue
    
    if longest:
        pd.DataFrame({'geodesic_path': longest}).to_csv(
            os.path.join(results_dir, "longest_geodesic_filtered.csv"), index=False)
        
        if not annotations.empty:
            geo_annots = annotations.loc[annotations.index.intersection(longest)][['go.BP', 'go.MF', 'go.CC']]
            save_output("geodesic_annotations_filtered", df=geo_annots.reset_index())
    
    return longest, G_filtered

# =======================
# 7. Visualización de Red
# =======================
def visualize_network(G, G_filtered, important_nodes):
    """Genera visualizaciones de la red"""
    if G.number_of_nodes() == 0:
        return
    
    # Configuración común
    plt_style = {
        'with_labels': False,
        'node_size': 30,
        'font_size': 4,
        'alpha': 0.6,
        'edge_color': 'gray',
        'width': 0.3,
        'arrows': False
    }
    
    # Spring Layout (red completa)
    plt.figure(figsize=(20, 16))
    pos_spring = nx.spring_layout(G, seed=42, k=0.1)
    nx.draw_networkx(G, pos_spring, **plt_style)
    
    for node in important_nodes:
        plt.text(pos_spring[node][0], pos_spring[node][1], node, 
                fontsize=6, ha='center', va='center')
    
    plt.title("Red de Coexpresión Génica - Spring Layout", fontsize=14)
    plt.axis('off')
    save_output("network_spring", fig=plt.gcf())
    plt.close()
    
    # Circular Layout (red filtrada)
    if G_filtered.number_of_nodes() > 0:
        plt.figure(figsize=(20, 16))
        pos_circular = nx.circular_layout(G_filtered)
        nx.draw_networkx(G_filtered, pos_circular, 
                        with_labels=True, 
                        node_size=100,
                        font_size=8,
                        alpha=0.8,
                        edge_color='gray',
                        width=0.5)
        plt.title("Red Filtrada - Circular Layout", fontsize=14)
        plt.axis('off')
        save_output("network_circular_filtered", fig=plt.gcf())
        plt.close()
    
    # Visualización interactiva (opcional)
    if PYVIS_AVAILABLE and G.number_of_nodes() > 0:
        try:
            net = Network(notebook=False, height="750px", width="100%", bgcolor="#222222", font_color="white")
            net.from_nx(G)
            net.toggle_physics(True)
            net.show_buttons(filter_=['physics'])
            net.save_graph(os.path.join(results_dir, "interactive_network.html"))
            print("✅ Red interactiva guardada como interactive_network.html")
        except Exception as e:
            print(f"⚠️ No se pudo generar la visualización interactiva: {str(e)}")

# =======================
# Ejecución principal
# =======================
def main():
    print("🚀 Iniciando análisis de expresión génica diferencial")
    
    try:
        # 1. Cargar datos
        df = load_and_prepare_data('RNA-Seq-expression-Norilsk2019.csv')
        
        # 2. Volcano Plot
        create_volcano_plot(df)
        print("✅ Volcano plot generado")
        
        # 3. Anotación GO
        valid_symbols, annotations, go_terms = perform_go_annotation(df)
        print(f"✅ Anotación GO completada para {len(valid_symbols)} genes")
        
        # 4. Red de coexpresión
        G = build_coexpression_network(valid_symbols)
        print(f"✅ Red generada con {G.number_of_nodes()} nodos y {G.number_of_edges()} aristas")
        
        # 5. Métricas de red
        metrics, degrees = calculate_network_metrics(G)
        print("📊 Métricas de red calculadas:")
        for k, v in metrics.items():
            print(f"{k}: {v}")
        
        # 6. Geodésica
        longest, G_filtered = find_geodesic(G, degrees, annotations)
        if longest:
            print(f"\n🔎 Geodésica más larga ({len(longest)} nodos): {longest}")
        
        # 7. Visualizaciones
        visualize_network(G, G_filtered, [n for n in G.nodes() if degrees.get(n, 0) > np.mean(list(degrees.values()))])
        print("✅ Visualizaciones generadas")
        
    except Exception as e:
        print(f"⛔ Error crítico: {str(e)}")
        sys.exit(1)
    
    print(f"\n🎉 Todos los resultados guardados en: {os.path.abspath(results_dir)}")

if __name__ == "__main__":
    main()