#!/usr/bin/env python3
from Bio.Blast import NCBIXML
from Bio import SeqIO
from pathlib import Path

blast_dir = Path("blast_results")
output_dir = Path("top_hits")
output_dir.mkdir(exist_ok=True)

out_fasta = output_dir / "unknown_blast.fasta"

with open(out_fasta, "w") as fasta_out:
    for blast_file in blast_dir.glob("*.xml"):
        with open(blast_file) as result_handle:
            blast_record = next(NCBIXML.parse(result_handle))
            top_hit = blast_record.alignments[0] if blast_record.alignments else None

            if top_hit:
                title = top_hit.hit_def.split(" >")[0]
                sequence = top_hit.hsps[0].sbjct
                header = f">{blast_file.stem}_{title[:50]}"
                fasta_out.write(f"{header}\n{sequence}\n")
