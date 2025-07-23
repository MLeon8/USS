#!/usr/bin/env python3
"""
Extrae los mejores hits de los resultados BLAST (.xml) y los guarda en formato FASTA.
Requiere: Biopython
"""

from Bio.Blast import NCBIXML
import os

# Directorios de entrada y salida
blast_dir = "blast_results"
output_dir = "top_hits"
os.makedirs(output_dir, exist_ok=True)

for filename in os.listdir(blast_dir):
    if filename.endswith(".xml"):
        filepath = os.path.join(blast_dir, filename)
        with open(filepath) as result_handle:
            blast_record = NCBIXML.read(result_handle)
            output_path = os.path.join(output_dir, filename.replace(".xml", ".fasta"))
            with open(output_path, "w") as fasta_out:
                for i, alignment in enumerate(blast_record.alignments[:3]):  # Top 3 hits
                    for hsp in alignment.hsps:
                        fasta_out.write(f">hit_{i+1}_{alignment.hit_id}\n")
                        fasta_out.write(f"{hsp.sbjct}\n")
                        break  # Solo 1 HSP por hit
        print(f"[✓] {output_path} generado.")
