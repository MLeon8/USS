#!/bin/bash
# Script modificado para identificación de patógeno

# ---------------------------
# CONFIGURACIÓN INICIAL
# ---------------------------
WORKDIR="$PWD/identificacion_patogeno"
SEQUENCES_FILE="$PWD/sequences.fasta"
THREADS=8
EVALUE=1e-10  # Valor menos estricto inicialmente

# Verificar dependencias
check_dependency() {
    if ! command -v $1 &> /dev/null; then
        echo "Error: $1 no está instalado. Por favor instálelo primero."
        exit 1
    fi
}

check_dependency prokka
check_dependency blastn
check_dependency clustalo
check_dependency fasttree
check_dependency efetch
check_dependency pandoc

# Crear estructura de directorios
mkdir -p $WORKDIR/{anotacion_secuencias,analisis_genomico,reporte_final}

# ---------------------------
# 1. PREPARACIÓN DE DATOS
# ---------------------------
echo "=== PASO 1: Preparación de datos ==="

# Verificar si el archivo de secuencias existe
if [ ! -f "$SEQUENCES_FILE" ]; then
    echo "Error: No se encuentra el archivo $SEQUENCES_FILE"
    exit 1
fi

# Dividir el archivo multifasta en secuencias individuales
echo "Dividiendo archivo multifasta..."
awk '/^>/ {filename=sprintf("seq_%02d.fasta",++i); print > filename; next} {print >> filename}' $SEQUENCES_FILE
mv seq_*.fasta $WORKDIR/anotacion_secuencias/

# ---------------------------
# 2. ANOTACIÓN DE SECUENCIAS
# ---------------------------
echo -e "\n=== PASO 2: Anotación de secuencias ==="

for seq in $WORKDIR/anotacion_secuencias/seq_*.fasta; do
    SEQNAME=$(basename $seq .fasta)
    SEQDIR="$WORKDIR/anotacion_secuencias/$SEQNAME"
    mkdir -p $SEQDIR
    
    echo "Anotando $SEQNAME..."
    mv $seq $SEQDIR/
    
    # Anotación con Prokka para virus
    echo "Ejecutando Prokka..."
    prokka \
        --outdir $SEQDIR/prokka \
        --prefix $SEQNAME \
        --kingdom Viruses \
        --evalue 1e-05 \
        --coverage 70 \
        --cpus $THREADS \
        $SEQDIR/$SEQNAME.fasta || echo "Prokka falló para $SEQNAME, continuando..."
    
    # BLAST online como alternativa
    echo "Ejecutando BLAST online para $SEQNAME..."
    blastn \
        -query $SEQDIR/$SEQNAME.fasta \
        -remote \
        -db nt \
        -out $SEQDIR/blast_results.xml \
        -outfmt 5 \
        -evalue $EVALUE \
        -max_target_seqs 10 || echo "BLAST online falló para $SEQNAME, continuando..."
    
    # Convertir resultados BLAST a tabla si existe
    if [ -f "$SEQDIR/blast_results.xml" ]; then
        blast_formatter \
            -archive $SEQDIR/blast_results.xml \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
            -out $SEQDIR/blast_results.tsv
        
        # Filtrar top hits
        if [ -f "$SEQDIR/blast_results.tsv" ]; then
            sort -k12,12nr $SEQDIR/blast_results.tsv | head -n 5 > $SEQDIR/top_hits.tsv
            echo "Top hit para $SEQNAME: $(head -n1 $SEQDIR/top_hits.tsv | cut -f13)"
        fi
    fi
done

# ---------------------------
# 3. ALINEAMIENTOS MÚLTIPLES
# ---------------------------
echo -e "\n=== PASO 3: Alineamientos múltiples ==="

# Combinar todas las secuencias
find $WORKDIR/anotacion_secuencias/ -name "*.fasta" -exec cat {} + > $WORKDIR/todas_secuencias.fasta 2>/dev/null

# Verificar que hay secuencias
if [ ! -s "$WORKDIR/todas_secuencias.fasta" ]; then
    echo "Error: No se encontraron secuencias para alinear"
    exit 1
fi

# Alineamiento con CLUSTAL Omega
echo "Realizando alineamiento múltiple..."
clustalo \
    -i $WORKDIR/todas_secuencias.fasta \
    -o $WORKDIR/alineamiento_multiple.aln \
    --outfmt=clu \
    --iter=3 \
    --threads=$THREADS || echo "Clustalo falló, continuando..."

# Generar árbol filogenético si existe el alineamiento
if [ -f "$WORKDIR/alineamiento_multiple.aln" ]; then
    echo "Generando árbol filogenético..."
    fasttree \
        -nt \
        -gtr \
        $WORKDIR/alineamiento_multiple.aln > $WORKDIR/arbol_filogenetico.nwk 2>/dev/null
fi

# ---------------------------
# 4. ANÁLISIS GENÓMICO
# ---------------------------
echo -e "\n=== PASO 4: Análisis genómico ==="

# Intentar identificar el mejor hit global
BEST_HIT=$(find $WORKDIR/anotacion_secuencias/ -name "top_hits.tsv" -exec cat {} + | sort -k12,12nr | head -n1 | cut -f2)

if [ -z "$BEST_HIT" ]; then
    echo "No se pudo identificar un hit significativo, usando Norovirus GII.4 como referencia"
    BEST_HIT="MG893928.1"  # Norovirus GII.4 como fallback
fi

echo "Usando genoma de referencia: $BEST_HIT"

# Descargar genoma completo
echo "Descargando genoma $BEST_HIT..."
efetch -db nucleotide -id $BEST_HIT -format fasta > $WORKDIR/analisis_genomico/genoma_referencia.fna || \
(echo "Falló efetch, usando genoma alternativo"; \
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$BEST_HIT&rettype=fasta" > $WORKDIR/analisis_genomico/genoma_referencia.fna)

# Anotación del genoma de referencia si se descargó
if [ -s "$WORKDIR/analisis_genomico/genoma_referencia.fna" ]; then
    echo "Anotando genoma de referencia..."
    prokka \
        --outdir $WORKDIR/analisis_genomico/prokka_annotation \
        --prefix genoma_referencia \
        --kingdom Viruses \
        --evalue 1e-05 \
        --cpus $THREADS \
        $WORKDIR/analisis_genomico/genoma_referencia.fna || echo "Prokka falló para el genoma de referencia"
else
    echo "No se pudo obtener genoma de referencia"
fi

# ---------------------------
# 5. GENERACIÓN DE REPORTES
# ---------------------------
echo -e "\n=== PASO 5: Generación de reportes ==="

# Crear reporte consolidado
echo "Generando reporte final..."
{
echo "=== REPORTE DE IDENTIFICACIÓN DE PATÓGENO ==="
echo "Fecha: $(date)"
echo -e "\n[1] SECUENCIAS ANALIZADAS"
echo "------------------------"
echo "Número de secuencias: $(grep -c ">" $SEQUENCES_FILE)"
echo "Longitud promedio: $(awk '/^>/ {getline; print length}' $SEQUENCES_FILE | awk '{sum+=$1} END {print sum/NR}') bp"

echo -e "\n[2] RESULTADOS DE BLAST"
echo "------------------------"
find $WORKDIR/anotacion_secuencias/ -name "top_hits.tsv" -exec cat {} + | sort -k12,12nr | head -n10 | awk -F'\t' '{print $1"\t"$2"\t"$3"%\t"$12"\t"$13}'

echo -e "\n[3] ANÁLISIS FILOGENÉTICO"
echo "-------------------------"
if [ -f "$WORKDIR/arbol_filogenetico.nwk" ]; then
    echo "Árbol filogenético generado: $WORKDIR/arbol_filogenetico.nwk"
else
    echo "No se pudo generar el árbol filogenético"
fi

echo -e "\n[4] ANOTACIÓN GENÓMICA"
echo "----------------------"
if [ -f "$WORKDIR/analisis_genomico/prokka_annotation/genoma_referencia.gff" ]; then
    echo "Genes identificados en la referencia:"
    grep "CDS" $WORKDIR/analisis_genomico/prokka_annotation/genoma_referencia.gff | cut -f9 | cut -d';' -f1 | sort | uniq -c
else
    echo "No se pudo obtener la anotación genómica"
fi

echo -e "\n[5] CONCLUSIONES"
echo "----------------"
BEST_MATCH=$(find $WORKDIR/anotacion_secuencias/ -name "top_hits.tsv" -exec cat {} + | sort -k12,12nr | head -n1 | cut -f13)
echo "El análisis sugiere que el agente patógeno es compatible con:"
echo "$BEST_MATCH"
} > $WORKDIR/reporte_final/resumen_analisis.txt

# Convertir a PDF
if command -v pandoc &> /dev/null; then
    pandoc $WORKDIR/reporte_final/resumen_analisis.txt -o $WORKDIR/reporte_final/identificacion_patogeno.pdf --pdf-engine=wkhtmltopdf || \
    pandoc $WORKDIR/reporte_final/resumen_analisis.txt -o $WORKDIR/reporte_final/identificacion_patogeno.html
fi

echo -e "\n=== ANÁLISIS COMPLETADO ==="
echo "Resultados guardados en: $WORKDIR"
tree $WORKDIR
