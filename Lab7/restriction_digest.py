"""
Lab7 - Exercises 1, 2 and 3: restriction digest

1. Take an arbitrary DNA sequence from NCBI between 1000 and 3000 nucleotides.
   -> pUC19c cloning vector, NCBI L09137.2 (2,686 bp)
2. Use restriction enzymes (name, recognition sequence, cleavage site).
   -> table ENZYMES in restriction_enzymes.py
3. Input:  the recognition sequence of each restriction enzyme
           a DNA sequence to be digested (cut into pieces)
   Output: the number of cleavages for each restriction enzyme
           the cleavage positions and the length of the fragments
           a simulation of the electrophoresis gel based on the number of enzymes used

Run:
    py restriction_digest.py                                  (pUC19 + EcoRI, BamHI, HindIII)
    py restriction_digest.py pUC19c_L09137.2.fasta EcoRI PstI G^ANTC
Enzymes can be given by name (from the table) or as a recognition sequence with "^".
"""

import sys

from restriction_enzymes import (ENZYMES, read_fasta, get_enzyme, find_cleavages,
                                 fragment_lengths, digest, draw_gel)


def print_enzyme_table():
    print(f"{'Enzyme':<9}{'Recognition seq':<17}{'Cleavage site':<15}")
    for name, site in ENZYMES.items():
        print(f"{name:<9}{site.replace('^', ''):<17}{site:<15}")
    print()


def run_digest(fasta, enzyme_inputs, gel_file, results_file, title=None, show=True):
    header, sequence = read_fasta(fasta)
    enzymes = [get_enzyme(e) for e in enzyme_inputs]
    lines = [f"Sequence: {header}", f"Length: {len(sequence)} bp", ""]
    lanes = []

    for name, site in enzymes:
        positions = find_cleavages(sequence, site)
        fragments = fragment_lengths(len(sequence), positions)
        lanes.append((name, fragments))
        lines.append(f"{name} ({site}): {len(positions)} cleavage(s)")
        lines.append(f"  cleavage positions: {positions}")
        lines.append(f"  fragments (bp): {fragments}")
        lines.append("")

    if len(enzymes) > 1:
        positions, fragments = digest(sequence, [site for _, site in enzymes])
        label = "+".join(name for name, _ in enzymes)
        lanes.append(("all", fragments))
        lines.append(f"{label} together: {len(positions)} cleavage(s)")
        lines.append(f"  cleavage positions: {positions}")
        lines.append(f"  fragments (bp): {fragments}")
        lines.append("")

    text = "\n".join(lines)
    print(text)
    with open(results_file, "w") as f:
        f.write(text)

    draw_gel(lanes, title or f"Restriction digest - {header.split()[0]}", gel_file, show=show)
    return lanes


def main():
    print_enzyme_table()
    fasta = sys.argv[1] if len(sys.argv) > 1 else "pUC19c_L09137.2.fasta"
    enzymes = sys.argv[2:] if len(sys.argv) > 2 else ["EcoRI", "BamHI", "HindIII"]
    run_digest(fasta, enzymes, "digest_gel.png", "digest_results.txt")


if __name__ == "__main__":
    main()
