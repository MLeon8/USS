import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
import gseapy as gp
from gseapy.plot import barplot
from mygene import MyGeneInfo

def cargar_datos(input_path):
    """Carga y prepara los datos con manejo de valores extremos"""
    df = pd.read_csv(input_path, sep=',')
    df = df[['Gene', 'log_2 fold change', 'Adjusted p-value']].copy()
    df.columns = ['gen_id', 'log2fc', 'pvalor_ajustado']
    
    # Manejo de p-values iguales a 0
    min_p = df.loc[df['pvalor_ajustado'] > 0, 'pvalor_ajustado'].min()
    df['pvalor_ajustado'] = df['pvalor_ajustado'].replace(0, min_p/100)
    
    return df

def generar_volcano_plot(df, pval_umbral=1e-5, fc_umbral=1, output_file='volcano_plot.png'):
    """Genera volcano plot con parámetros estándar"""
    plt.figure(figsize=(12, 8))
    
    df['-log10_p'] = -np.log10(df['pvalor_ajustado'])
    df['significativo'] = (abs(df['log2fc']) >= fc_umbral) & (df['pvalor_ajustado'] <= pval_umbral)
    
    # Scatter plot
    plt.scatter(
        x=df['log2fc'],
        y=df['-log10_p'],
        c=np.where(df['significativo'], '#e41a1c', '#999999'),
        alpha=0.6,
        s=25
    )
    
    # Líneas de umbral
    plt.axvline(x=fc_umbral, color='#377eb8', linestyle='--', linewidth=1)
    plt.axvline(x=-fc_umbral, color='#377eb8', linestyle='--', linewidth=1)
    plt.axhline(y=-np.log10(pval_umbral), color='#377eb8', linestyle='--', linewidth=1)
    
    # Etiquetado de top genes
    top_genes = pd.concat([
        df[df['significativo']].nlargest(5, 'log2fc'),
        df[df['significativo']].nsmallest(5, 'log2fc')
    ])
    
    texts = []
    for _, row in top_genes.iterrows():
        texts.append(plt.text(
            x=row['log2fc'],
            y=row['-log10_p'],
            s=row['gen_id'],
            fontsize=8,
            ha='center'
        ))
    
    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
    
    plt.title('Volcano Plot: Expresión Diferencial\n(FDR < 1e-5, |log2FC| ≥ 1)')
    plt.xlabel('log2(Fold Change)')
    plt.ylabel('-log10(p-value ajustado)')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return df[df['significativo']]

def convertir_a_symbols(ensembl_ids):
    """Convierte IDs ENSEMBL a símbolos de genes usando MyGene"""
    mg = MyGeneInfo()
    # Consultamos en lotes para evitar límites de la API
    batch_size = 1000
    results = []
    
    for i in range(0, len(ensembl_ids), batch_size):
        batch = ensembl_ids[i:i+batch_size]
        try:
            batch_results = mg.querymany(batch, scopes='ensembl.gene', fields='symbol', species='human')
            results.extend(batch_results)
        except Exception as e:
            print(f"Error en conversión de IDs: {str(e)}")
            continue
    
    # Procesar resultados
    symbol_map = {}
    for item in results:
        if 'symbol' in item:
            symbol_map[item['query']] = item['symbol']
    
    return symbol_map

def analisis_go(genes_sig, output_dir='GO_analysis'):
    """Realiza análisis de enriquecimiento GO con mejor manejo de IDs y parámetros optimizados"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Convertir ENSEMBL IDs a símbolos de genes
    print("Convirtiendo IDs ENSEMBL a símbolos de genes...")
    ensembl_ids = genes_sig['gen_id'].tolist()
    symbol_map = convertir_a_symbols(ensembl_ids)
    
    # Filtrar genes con símbolos válidos
    valid_genes = [symbol_map.get(gid, '') for gid in ensembl_ids]
    valid_genes = [g for g in valid_genes if g]  # Eliminar entradas vacías
    
    if not valid_genes:
        print("Error: No se pudieron convertir los IDs ENSEMBL a símbolos de genes")
        return None
    
    print(f"\nGenes convertidos a símbolos: {len(valid_genes)}/{len(ensembl_ids)}")
    
    # Bases de datos GO actuales en Enrichr
    go_databases = {
        'Biological_Process': 'GO_Biological_Process_2023',
        'Molecular_Function': 'GO_Molecular_Function_2023',
        'Cellular_Component': 'GO_Cellular_Component_2023',
        'KEGG_2021_Human': 'KEGG_2021_Human'  # Añadimos KEGG como alternativa
    }
    
    all_results = []
    
    for go_type, db_name in go_databases.items():
        try:
            print(f"\nAnalizando {go_type}...")
            
            # Análisis GO con parámetros optimizados
            enr = gp.enrichr(
                gene_list=valid_genes,
                gene_sets=[db_name],
                organism='human',
                outdir=None,
                cutoff=0.1  # Umbral inicial más relajado
            )
            
            if enr.results is not None and not enr.results.empty:
                # Filtrar por significancia con FDR < 0.25 inicialmente
                significant = enr.results[enr.results['Adjusted P-value'] < 0.25]
                
                if not significant.empty:
                    significant['GO_type'] = go_type
                    all_results.append(significant)
                    
                    print(f"Términos {go_type} encontrados: {len(significant)}")
                    print(significant[['Term', 'Adjusted P-value', 'Odds Ratio']].head(5))
                    
                    # Guardar resultados individuales
                    significant.to_csv(os.path.join(output_dir, f'GO_{go_type}.csv'), index=False)
                    
                    # Generar gráfico de barras para los términos más significativos
                    try:
                        plt.figure(figsize=(10, 6))
                        top_terms = significant.nsmallest(10, 'Adjusted P-value')
                        barplot(top_terms, title=f'Top {go_type} Terms')
                        plt.tight_layout()
                        plt.savefig(os.path.join(output_dir, f'GO_{go_type}_plot.png'), dpi=300)
                        plt.close()
                    except Exception as e:
                        print(f"Error generando gráfico para {go_type}: {str(e)}")
            
        except Exception as e:
            print(f"Error en análisis {go_type}: {str(e)}")
            continue
    
    # Combinar todos los resultados
    if all_results:
        final_results = pd.concat(all_results).sort_values('Adjusted P-value')
        final_results.to_csv(os.path.join(output_dir, 'all_GO_results.csv'), index=False)
        
        # Generar gráfico combinado de los mejores términos
        try:
            plt.figure(figsize=(12, 8))
            top_terms = final_results.nsmallest(15, 'Adjusted P-value')
            
            colors = {
                'Biological_Process': '#1f77b4',
                'Molecular_Function': '#ff7f0e',
                'Cellular_Component': '#2ca02c',
                'KEGG_2021_Human': '#d62728'
            }
            
            for i, (_, row) in enumerate(top_terms.iterrows()):
                plt.barh(
                    y=i,
                    width=-np.log10(row['Adjusted P-value']),
                    color=colors[row['GO_type']],
                    label=row['GO_type'] if i == 0 else ""
                )
                plt.text(
                    -np.log10(row['Adjusted P-value']) + 0.1,
                    i,
                    f"{row['Term']} ({row['GO_type']})",
                    va='center'
                )
            
            plt.yticks([])
            plt.xlabel('-log10(FDR)')
            plt.title('Términos más enriquecidos')
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'combined_enrichment_plot.png'), dpi=300)
            plt.close()
        except Exception as e:
            print(f"Error generando gráfico combinado: {str(e)}")
        
        return final_results
    else:
        print("\nNo se encontraron términos significativos en ninguna categoría")
        
        # Guardar lista de genes para análisis manual
        pd.DataFrame({'ENSEMBL_ID': ensembl_ids, 'Gene_Symbol': [symbol_map.get(gid, '') for gid in ensembl_ids]}).to_csv(
            os.path.join(output_dir, 'gene_symbol_mapping.csv'), index=False)
        
        return None

def main():
    # Configuración
    input_path = 'RNA-Seq-expression-Norilsk2019.csv'
    output_dir = 'resultados_RNAseq'
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Cargar datos
    print("Cargando datos...")
    df = cargar_datos(input_path)
    print(f"Datos cargados: {len(df)} genes")
    
    # 2. Generar Volcano Plot
    print("\nGenerando Volcano Plot...")
    genes_sig = generar_volcano_plot(
        df,
        pval_umbral=1e-5,
        fc_umbral=1,
        output_file=os.path.join(output_dir, 'volcano_plot.png')
    )
    print(f"Genes significativos identificados: {len(genes_sig)}")
    
    # 3. Análisis de enriquecimiento (GO y KEGG)
    print("\nRealizando análisis de enriquecimiento funcional...")
    
    # Analizar genes up y down por separado
    genes_up = genes_sig[genes_sig['log2fc'] > 0]
    genes_down = genes_sig[genes_sig['log2fc'] < 0]
    
    if len(genes_up) > 0:
        print("\nAnalizando genes sobreexpresados...")
        go_up = analisis_go(genes_up, os.path.join(output_dir, 'GO_upregulated'))
    else:
        print("\nNo hay genes sobreexpresados para analizar")
    
    if len(genes_down) > 0:
        print("\nAnalizando genes subexpresados...")
        go_down = analisis_go(genes_down, os.path.join(output_dir, 'GO_downregulated'))
    else:
        print("\nNo hay genes subexpresados para analizar")
    
    print(f"\nAnálisis completado. Resultados guardados en: {os.path.abspath(output_dir)}")

if __name__ == "__main__":
    main()
