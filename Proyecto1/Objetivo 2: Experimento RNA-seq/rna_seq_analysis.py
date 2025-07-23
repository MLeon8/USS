import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from gprofiler import GProfiler

# --- Cargar datos
df = pd.read_csv("RNA-Seq-expression-Norilsk2019.csv")
print(df.head())
print(df.columns)

# Renombrar columnas para facilitar el análisis (opcional pero recomendable)
df.rename(columns={
    'log_2 fold change': 'log2FoldChange',
    'Adjusted p-value': 'pvalue',
    'Gene': 'gene_id'
}, inplace=True)

# Parámetros de corte
logfc_cutoff = 1
pval_cutoff = 0.05

# Cálculo de -log10(pvalue)
df['-log10(pval)'] = -np.log10(df['pvalue'].replace(0, 1e-300))  # evitar log(0)

# Marcamos genes significativos
df['significant'] = (df['log2FoldChange'] > logfc_cutoff) & (df['pvalue'] < pval_cutoff)

# --- VOLCANO PLOT
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x='log2FoldChange', y='-log10(pval)', hue='significant',
                palette={True: 'red', False: 'grey'}, alpha=0.6)
plt.axhline(-np.log10(pval_cutoff), color='blue', linestyle='--')
plt.axvline(logfc_cutoff, color='green', linestyle='--')
plt.title('Volcano Plot - Genes diferencialmente expresados')
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(p-value)')
plt.legend(title='Significativo')
plt.tight_layout()
plt.savefig("volcano_plot.png")
plt.show()

# --- Extraer genes sobreexpresados
overexpressed_genes = df[(df['log2FoldChange'] > logfc_cutoff) & (df['pvalue'] < pval_cutoff)]
print(f"Genes sobreexpresados detectados: {len(overexpressed_genes)}")

# --- Análisis funcional con GO (Gene Ontology)
gp = GProfiler(return_dataframe=True)
result = gp.profile(organism='hsapiens', query=overexpressed_genes['gene_id'].tolist())

# Guardar los resultados
result.to_csv("GO_annotated_overexpressed_genes.csv", index=False)
print("Anotación GO completada y guardada.")
print(result[['native', 'name', 'p_value']].head())
