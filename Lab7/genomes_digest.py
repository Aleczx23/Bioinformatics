"""
Lab7 - Requirement 3 on two genomes

Take 2 genome sequences from NCBI and use 3 restriction enzymes on them.
For each genome: number of cleavages per enzyme, cleavage positions,
fragment lengths and one electrophoresis gel.

Genomes (two versions of the same virus):
- Monkeypox virus clade I,   strain Zaire_1979-005, NCBI DQ011155.1 (196,967 bp)
- Monkeypox virus clade IIb, 2022 outbreak,         NCBI NC_063383.1 (197,209 bp)
"""

from restriction_digest import run_digest

ENZYMES_USED = ["EcoRI", "BamHI", "HindIII"]
GENOMES = [
    ("MPXV_cladeI_DQ011155.1.fasta", "cladeI"),
    ("MPXV_cladeIIb_NC_063383.1.fasta", "cladeIIb"),
]


def main():
    for fasta, short in GENOMES:
        run_digest(fasta, ENZYMES_USED,
                   gel_file=f"gel_{short}.png",
                   results_file=f"results_{short}.txt",
                   title=f"Monkeypox virus {short} - EcoRI, BamHI, HindIII")


if __name__ == "__main__":
    main()
