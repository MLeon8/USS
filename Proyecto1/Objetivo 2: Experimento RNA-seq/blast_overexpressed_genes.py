#!/usr/bin/env python3
# Script: blast_overexpressed_genes.py
# Requiere: Biopython, BLAST+ instalado localmente

from Bio.Blast import NCBIWWW
from Bio import SeqIO
import time

fasta_file = "overexpressed_sequences.fasta"  # Asegúrate de tener este archivo con las secuencias
output_dir = "blast_results"

import os
os.makedirs(output_dir, exist_ok=True)

for record in SeqIO.parse(fasta_file, "fasta"):
    print(f"Enviando BLAST para: {record.id}")
    try:
        result_handle = NCBIWWW.qblast(program="blastn", database="nt", sequence=record.seq)
        output_path = os.path.join(output_dir, f"{record.id}_blast.xml")
        with open(output_path, "w") as out_handle:
            out_handle.write(result_handle.read())
        print(f"Resultado guardado en: {output_path}")
        time.sleep(3)  # evitar saturar servidor NCBI
    except Exception as e:
        print(f"Error con {record.id}: {e}")