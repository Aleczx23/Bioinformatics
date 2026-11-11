"""
Lab7 - Restriction enzymes (shared code)

A restriction enzyme recognizes a short DNA pattern (the recognition sequence)
and cuts the DNA at a fixed place inside it (the cleavage site).
Here a cleavage site is written with "^", for example EcoRI = G^AATTC:
the enzyme cuts between G and A on the top strand.

This file holds:
- the table of enzymes (name, recognition sequence, cleavage site)
- reading a FASTA file
- finding the cleavage positions and the fragment lengths
- drawing an electrophoresis gel with one column per digest
"""

import math
import re

import matplotlib.pyplot as plt

# name: recognition sequence with "^" at the cleavage site (top strand)
ENZYMES = {
    "EcoRI":   "G^AATTC",
    "BamHI":   "G^GATCC",
    "HindIII": "A^AGCTT",
    "PstI":    "CTGCA^G",
    "SmaI":    "CCC^GGG",
    "KpnI":    "GGTAC^C",
    "SacI":    "GAGCT^C",
    "XbaI":    "T^CTAGA",
    "SalI":    "G^TCGAC",
    "EcoRV":   "GAT^ATC",
    "NotI":    "GC^GGCCGC",
    "HinfI":   "G^ANTC",      # free pattern: N = any base
}

# IUPAC letters, so a recognition sequence can be a "free pattern" (e.g. G^ANTC)
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def read_fasta(path):
    header = ""
    sequence = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            else:
                sequence.append(line.upper())
    return header, "".join(sequence)


def parse_site(site):
    """'G^AATTC' -> ('GAATTC', 1). Without '^' the cut is put in the middle."""
    site = site.upper()
    if "^" in site:
        cut = site.index("^")
        recseq = site.replace("^", "")
    else:
        recseq = site
        cut = len(recseq) // 2
    return recseq, cut


def get_enzyme(name_or_site):
    """Accepts an enzyme name from the table (EcoRI) or a recognition sequence (G^AATTC)."""
    for name, site in ENZYMES.items():
        if name.lower() == name_or_site.lower():
            return name, site
    return name_or_site.upper(), name_or_site.upper()


def find_cleavages(sequence, site):
    """Positions where the enzyme cuts (number of bases before the cut, top strand)."""
    recseq, cut = parse_site(site)
    pattern = "".join(IUPAC[b] for b in recseq)
    # (?=...) also finds overlapping matches
    return [m.start() + cut for m in re.finditer(f"(?=({pattern}))", sequence)]


def fragment_lengths(seq_length, positions):
    """Linear DNA: n cuts give n + 1 fragments."""
    edges = [0] + sorted(positions) + [seq_length]
    return [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]


def digest(sequence, sites):
    """Cut with several enzymes at once (double / triple digest)."""
    positions = set()
    for site in sites:
        positions.update(find_cleavages(sequence, site))
    positions = sorted(positions)
    return positions, fragment_lengths(len(sequence), positions)


# ---------------------------------------------------------------- gel

SMALL_LADDER = [3000, 2500, 2000, 1500, 1000, 750, 500, 250, 100]
LARGE_LADDER = [48500, 24000, 15000, 10000, 8000, 6000, 4000, 3000, 2000, 1500, 1000, 500]


def choose_ladder(lanes):
    biggest = max(max(lengths) for _, lengths in lanes)
    return SMALL_LADDER if biggest <= 3000 else LARGE_LADDER


def distance(length, ladder):
    """Short fragments travel far, long fragments stay near the well (log scale)."""
    top, bottom = max(ladder), min(ladder)
    pos = (math.log10(top) - math.log10(length)) / (math.log10(top) - math.log10(bottom))
    return 1 + 8 * pos


def draw_gel(lanes, title, filename, colors=None, show=True):
    """lanes = [(label, [fragment lengths]), ...]; one column per lane plus the ladder L."""
    ladder = choose_ladder(lanes)
    n = len(lanes)
    fig, ax = plt.subplots(figsize=(max(6, 1.3 * n + 3), 8))
    ax.set_facecolor("black")
    ax.set_xlim(-1.2, n + 0.6)
    ax.set_ylim(10.3, -0.3)
    ax.set_xticks(range(n + 1))
    ax.set_xticklabels(["L"] + [label for label, _ in lanes], fontsize=8)
    ax.set_ylabel("distance from the well (cm)")
    ax.set_title(title)

    for x in range(n + 1):                       # wells
        ax.plot([x - 0.35, x + 0.35], [0.4, 0.4], color="dimgray", linewidth=6)

    for size in ladder:
        y = distance(size, ladder)
        ax.plot([-0.3, 0.3], [y, y], color="deepskyblue", linewidth=3)
        ax.text(-0.42, y, f"{size}", color="white", fontsize=7, ha="right", va="center")

    for lane, (label, lengths) in enumerate(lanes, start=1):
        color = colors[lane - 1] if colors else "orange"
        biggest = max(lengths)
        for size in lengths:
            y = min(distance(max(size, 1), ladder), 10.1)   # very small pieces run to the bottom
            # brighter band = more DNA (longer fragment)
            alpha = 0.35 + 0.65 * math.sqrt(size / biggest)
            ax.plot([lane - 0.3, lane + 0.3], [y, y], color=color, linewidth=3, alpha=alpha)
        if len(lengths) <= 8:                     # label sizes only when the lane is readable
            tiny = sorted(s for s in set(lengths) if s < min(ladder))
            for size in set(lengths) - set(tiny):
                ax.text(lane + 0.33, distance(size, ladder), str(size),
                        color="white", fontsize=6, va="center")
            if tiny:                              # below the ladder: one label at the bottom
                ax.text(lane + 0.33, 10.1, ", ".join(map(str, tiny)),
                        color="white", fontsize=6, va="center")

    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    print(f"Gel saved as {filename}")
    if show:
        plt.show()
    plt.close(fig)
