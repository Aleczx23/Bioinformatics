"""
Lab7 - Homework

Digest 2 genomes from 2 versions of the same virus and show the digested
segments on separate columns of the same electrophoresis gel.

Virus: Monkeypox virus (mpox)
- version 1: clade I,   strain Zaire_1979-005 (NCBI DQ011155.1)
- version 2: clade IIb, 2022 outbreak        (NCBI NC_063383.1)
The script uses restriction_enzymes.py and the two FASTA files from the Lab7 folder.
Run it from anywhere:  py homework_virus_versions.py

Enzymes: EcoRI, BamHI, HindIII -> the gel has two columns per enzyme (I, IIb) side by side,
so the differences between the two versions are easy to see.
"""

import os
import sys

# the shared code and the genome files are in the Lab7 folder, one level up
HERE = os.path.dirname(os.path.abspath(__file__))
LAB7 = os.path.dirname(HERE)
sys.path.insert(0, LAB7)

from restriction_enzymes import read_fasta, get_enzyme, find_cleavages, fragment_lengths, draw_gel

ENZYMES_USED = ["EcoRI", "BamHI", "HindIII"]
VERSIONS = [
    (os.path.join(LAB7, "MPXV_cladeI_DQ011155.1.fasta"), "I", "orange"),
    (os.path.join(LAB7, "MPXV_cladeIIb_NC_063383.1.fasta"), "IIb", "limegreen"),
]


def different_bands(a, b, tolerance=0.03, min_size=500):
    """Bands (>= 500 bp) of one version with no band of similar size (within 3%) in the other.
    A gel cannot tell apart e.g. 13773 bp and 13772 bp, so those count as the same band."""
    a = sorted((x for x in a if x >= min_size), reverse=True)
    b = sorted((x for x in b if x >= min_size), reverse=True)
    only_a = []
    for x in a:
        match = next((y for y in b if abs(x - y) <= tolerance * x), None)
        if match is None:
            only_a.append(x)
        else:
            b.remove(match)
    return only_a, b


def main():
    genomes = [(read_fasta(f)[0], read_fasta(f)[1], short, color) for f, short, color in VERSIONS]
    lanes, colors, lines = [], [], []

    for header, sequence, short, _ in genomes:
        lines.append(f"Clade {short}: {header} ({len(sequence)} bp)")
    lines.append("")

    for enzyme in ENZYMES_USED:
        name, site = get_enzyme(enzyme)
        frags = {}
        for header, sequence, short, color in genomes:
            positions = find_cleavages(sequence, site)
            frags[short] = fragment_lengths(len(sequence), positions)
            lanes.append((f"{name}\n{short}", frags[short]))
            colors.append(color)
            lines.append(f"{name} ({site}) clade {short}: {len(positions)} cleavage(s), "
                         f"{len(frags[short])} fragments")
            lines.append(f"  cleavage positions: {positions}")
            lines.append(f"  fragments (bp): {frags[short]}")

        only_i, only_iib = different_bands(frags["I"], frags["IIb"])
        lines.append(f"  bands only in clade I:   {only_i}")
        lines.append(f"  bands only in clade IIb: {only_iib}")
        lines.append("")

    text = "\n".join(lines)
    print(text)
    with open(os.path.join(HERE, "homework_results.txt"), "w") as f:
        f.write(text)

    draw_gel(lanes, "Monkeypox virus clade I (orange) vs clade IIb (green)",
             os.path.join(HERE, "homework_gel.png"), colors=colors)


if __name__ == "__main__":
    main()
