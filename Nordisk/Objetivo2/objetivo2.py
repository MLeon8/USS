import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
import gseapy as gp  # Para análisis de enriquecimiento GO

def cargar_datos(input_path):
    """Carga y prepara los datos con manejo de valores extremos"""
    df = pd.read_csv(input_path, sep=',')
    df = df[['Gene', 'log_2 fold change', 'Adjusted p-value']].copy()
    df.columns = ['gen_id', 'log2fc', 'pvalor_ajustado']
    
    # Manejo de p-values iguales a 0 (reemplazo con el valor más pequeño posible)
    min_p = df.loc[df['pvalor_ajustado'] > 0, 'pvalor_ajustado'].min()
    df['pvalor_ajustado'] = df['pvalor_ajustado'].replace(0, min_p/100)
    
    return df

def generar_volcano_plot(df, pval_umbral=1e-5, fc_umbral=1, output_file='volcano_plot.png'):
    """
    Genera volcano plot con justificación de umbrales
    
    Parámetros:
    - pval_umbral: 1e-5 (umbral estándar para estudios transcriptómicos)
    - fc_umbral: 1 (corresponde a un cambio de 2x en expresión)
    """
    plt.figure(figsize=(12, 8))
    
    # Cálculo de -log10(p-value)
    df['-log10_p'] = -np.log10(df['pvalor_ajustado'])
    
    # Identificación de genes significativos
    df['significativo'] = (abs(df['log2fc']) >= fc_umbral) & (df['pvalor_ajustado'] <= pval_umbral)
    
    # Filtrado de valores no finitos
    plot_data = df[np.isfinite(df['-log10_p']) & np.isfinite(df['log2fc'])]
    
    # Visualización
    plt.scatter(
        x=plot_data['log2fc'],
        y=plot_data['-log10_p'],
        c=np.where(plot_data['significativo'], '#e41a1c', '#999999'),
        alpha=0.6,
        s=25,
        edgecolor='none'
    )
    
    # Líneas de umbral con anotaciones explicativas
    plt.axvline(x=fc_umbral, color='#377eb8', linestyle='--', linewidth=1, alpha=0.7)
    plt.axvline(x=-fc_umbral, color='#377eb8', linestyle='--', linewidth=1, alpha=0.7)
    plt.axhline(y=-np.log10(pval_umbral), color='#377eb8', linestyle='--', linewidth=1, alpha=0.7)
    
    # Anotaciones explicativas de los umbrales
    plt.text(fc_umbral+0.1, plt.ylim()[1]*0.95, 
             f'FC ≥ {2**fc_umbral:.1f}x\n(Log2FC ≥ {fc_umbral})',
             ha='left', va='top', color='#377eb8')
    plt.text(-fc_umbral-0.1, plt.ylim()[1]*0.95, 
             f'FC ≤ {2**-fc_umbral:.2f}x\n(Log2FC ≤ -{fc_umbral})',
             ha='right', va='top', color='#377eb8')
    plt.text(plt.xlim()[1]*0.95, -np.log10(pval_umbral)+0.5, 
             f'p-ajustado < {pval_umbral:.0e}',
             ha='right', va='bottom', color='#377eb8')
    
    # Etiquetado de los top genes
    top_genes = pd.concat([
        df[df['significativo']].nlargest(5, 'log2fc'),
        df[df['significativo']].nsmallest(5, 'log2fc')
    ]).drop_duplicates()
    
    texts = []
    for _, row in top_genes.iterrows():
        texts.append(plt.text(
            x=row['log2fc'],
            y=row['-log10_p'],
            s=row['gen_id'],
            fontsize=8,
            ha='center',
            va='center'
        ))
    
    if texts:
        adjust_text(texts,
                   arrowprops=dict(arrowstyle='-', color='black', lw=0.5),
                   expand_points=(1.5, 1.5))
    
    # Configuración del gráfico
    plt.title('Volcano Plot: Expresión Diferencial en Fiebre de Norilsk', pad=15)
    plt.xlabel('log2(Fold Change)', labelpad=10)
    plt.ylabel('-log10(p-value ajustado)', labelpad=10)
    plt.grid(True, linestyle='--', alpha=0.2)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return df[df['significativo']]

def analisis_go(genes_sig, output_dir='resultados_GO'):
    """Realiza análisis de enriquecimiento de Gene Ontology"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Convertir ENSEMBL IDs a lista
    gene_list = genes_sig['gen_id'].tolist()
    
    # Análisis de enriquecimiento GO
    go_results = gp.enrichr(
        gene_list=gene_list,
        gene_sets=['GO_Biological_Process_2023', 'GO_Molecular_Function_2023', 'GO_Cellular_Component_2023'],
        organism='human',
        outdir=output_dir
    )
    
    # Procesar resultados
    if go_results is not None:
        # Guardar resultados completos
        go_results.results.to_csv(os.path.join(output_dir, 'resultados_GO_completos.csv'), index=False)
        
        # Filtrar los 5 términos más significativos por categoría
        top_go = go_results.results.sort_values('Adjusted P-value').groupby('Gene_set').head(5)
        top_go.to_csv(os.path.join(output_dir, 'top_GO_terms.csv'), index=False)
        
        return top_go
    else:
        print("No se encontraron términos GO significativos")
        return None

def main():
    # Configuración de rutas
    input_path = 'RNA-Seq-expression-Norilsk2019.csv'
    output_dir = 'resultados_RNAseq'
    os.makedirs(output_dir, exist_ok=True)
    
    # 2a. Volcano plot con justificación de umbrales
    print("Cargando datos...")
    df = cargar_datos(input_path)
    print(f"Datos cargados: {len(df)} genes")
    
    print("\nGenerando volcano plot...")
    genes_sig = generar_volcano_plot(
        df,
        pval_umbral=1e-5,
        fc_umbral=1,
        output_file=os.path.join(output_dir, 'volcano_plot.png')
    )
    print(f"Genes significativos identificados: {len(genes_sig)}")
    
    # 2b. Análisis de Gene Ontology
    print("\nRealizando análisis de Gene Ontology...")
    go_results = analisis_go(genes_sig[genes_sig['log2fc'] > 0],  # Solo genes sobreexpresados
                           os.path.join(output_dir, 'GO_analysis'))
    
    if go_results is not None:
        print("\nTérminos GO más enriquecidos:")
        print(go_results[['Gene_set', 'Term', 'Adjusted P-value']].head(10))
    
    print(f"\nResultados guardados en: {os.path.abspath(output_dir)}")

if __name__ == "__main__":
    main()
