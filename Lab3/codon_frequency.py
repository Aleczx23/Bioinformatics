
import gzip
import sys
from itertools import product
import matplotlib
matplotlib.use("Agg")  # save the chart to a file (no window needed)
import matplotlib.pyplot as plt

# Genetic code from the slide, used to show the amino acid of each codon
GENETIC_CODE = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}


def read_fasta(path):
    """Read the first record of a FASTA file -> (header, sequence as ARN)."""
    header, parts = "", []
    opener = gzip.open if path.endswith(".gz") else open   # .fna.gz from NCBI
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    break
                header = line[1:]
            else:
                parts.append(line)
    seq = "".join(parts).upper().replace("T", "U")   # DNA -> ARN
    return header, seq


def count_codons(seq):
    """Slide a 3-letter window one position at a time and count each codon."""
    counts = {"".join(c): 0 for c in product("UCAG", repeat=3)}   # all 64
    skipped = 0
    for i in range(len(seq) - 2):
        codon = seq[i:i + 3]
        if codon in counts:
            counts[codon] += 1
        else:
            skipped += 1        # windows with unknown letters (e.g. N)
    return counts, skipped


def make_chart(top10, total, title, out_file):
    labels = [f"{c}\n{GENETIC_CODE[c]}" for c, _ in top10]
    freqs = [n / total * 100 for _, n in top10]

    fig, ax = plt.subplots(figsize=(10, 5.5))
    bars = ax.bar(labels, freqs, color="#3a8f8a", width=0.65)
    for bar, f in zip(bars, freqs):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{f:.2f}%", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Frequency (% of all windows)")
    ax.set_title(f"Top 10 codons – {title}", fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(0, max(freqs) * 1.15)
    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    plt.close(fig)


def main(path):
    header, seq = read_fasta(path)
    counts, skipped = count_codons(seq)
    total = sum(counts.values())

    print(f"Genome      : {header}")
    print(f"Length      : {len(seq)} nucleotides")
    print(f"Windows read: {total}  (sliding window of 3, step 1)")
    if skipped:
        print(f"Skipped     : {skipped} windows with unknown letters")

    ranked = sorted(counts.items(), key=lambda x: x[1], reverse=True)

    print("\nAll 64 codons:")
    for i, (codon, n) in enumerate(ranked, 1):
        print(f"{i:2}. {codon} ({GENETIC_CODE[codon]:4}) {n:7}  {n / total * 100:5.2f}%")

    top10 = ranked[:10]
    print("\nTOP 10 codons:")
    for i, (codon, n) in enumerate(top10, 1):
        print(f"{i:2}. {codon} ({GENETIC_CODE[codon]:4}) {n:7}  {n / total * 100:5.2f}%")

    name = " ".join(header.split()[:4]) if header else path
    out_file = "top10_codons.png"
    make_chart(top10, total, name, out_file)
    print(f"\nChart saved to {out_file}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python codon_frequency.py genome.fasta")
        sys.exit(1)
    main(sys.argv[1])
