"""
Lab8 - DNA patterns in gene promoters (C+G% and Information Content)

Implement a software app that makes DNA patterns from DNA sequences in gene
promoters. Compute the C+G% and the information content values from each
sliding window.

Steps 1-7 (checking the algorithm on the test sequence):
1. Use the sequence S = CGGACTGATCTATCT(A...A)CGTAGCATCTATCGATCTAT CTAGCGATCTA
2. Use a sliding window of L = 30 bases
3. Function for the C+G content       (board: CG ~ 29.27)
4. Function for the information content (board: IC ~ 27.53)
5. Plot the pattern
6. Calculate the center of weight of the pattern
7. Take the center of each pattern and plot it on a 2nd chart
   (the 2nd chart with several promoters is made by promoter_genbank.py)

Run:
    py promoter_pattern.py
    py promoter_pattern.py ACGT...            (any other sequence)
"""

import os
import sys

import matplotlib.pyplot as plt

WINDOW = 30
# On the board the sequence continues on a second line (-> CTAGCGATCTA).
# With 22 A's in the middle this gives CG = 29.41 and IC = 27.55, the closest to the
# board values (~29.27 and ~27.53). Change NUMBER_OF_A to try other lengths.
NUMBER_OF_A = 22
S = "CGGACTGATCTATCT" + "A" * NUMBER_OF_A + "CGTAGCATCTATCGATCTAT" + "CTAGCGATCTA"


# ---------- Step 3: C+G content ----------
def cg_content(seq):
    """Percentage of C and G bases in the sequence."""
    seq = seq.upper()
    return (seq.count("C") + seq.count("G")) / len(seq) * 100


# ---------- Step 4: information content ----------
def information_content(seq):
    """
    Information content (Kappa index of coincidence), as in Gagniuc's algorithm:

        f(xi, xi+u) = 1 if xi == xi+u, else 0
        N  = |s| - 1
        IC = [ sum over u = 1..N of ( matches(u) / length(B) * 100 ) ] / N

    For every shift u, the sequence A is compared with B = A[u+1..N]
    (the copy shifted by u positions, without the last base); we take the %
    of positions where the bases match, then average these % over all shifts.
    Repetitive DNA (like AAAA...) gives a high IC, mixed DNA gives a low IC.
    """
    a = seq.upper()
    n = len(a) - 1
    total = 0.0
    for u in range(1, n + 1):
        b = a[u:n]                      # shifted copy, last base not included
        if not b:
            continue
        matches = sum(1 for i in range(len(b)) if a[i] == b[i])
        total += matches / len(b) * 100
    return total / n


# ---------- Step 2: sliding window ----------
def sliding_windows(seq, window=WINDOW):
    """All windows of length `window`, moved one base at a time."""
    return [seq[i:i + window] for i in range(len(seq) - window + 1)]


def make_pattern(seq, window=WINDOW):
    """The pattern = list of (C+G%, IC) points, one point per sliding window."""
    if len(seq) < window:
        raise ValueError(f"The sequence must have at least {window} bases.")
    return [(cg_content(w), information_content(w)) for w in sliding_windows(seq, window)]


# ---------- Step 6: center of weight ----------
def center_of_weight(pattern):
    """Average point of the pattern: (mean C+G%, mean IC)."""
    cg = sum(p[0] for p in pattern) / len(pattern)
    ic = sum(p[1] for p in pattern) / len(pattern)
    return cg, ic


# ---------- Step 5: plot the pattern ----------
def plot_pattern(pattern, title, filename, show=True):
    cg = [p[0] for p in pattern]
    ic = [p[1] for p in pattern]
    center = center_of_weight(pattern)
    positions = range(1, len(pattern) + 1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # left: the two signals along the sequence
    ax1.plot(positions, cg, color="tab:blue", label="C+G %")
    ax1.plot(positions, ic, color="tab:red", label="IC")
    ax1.set_xlabel(f"Window start position (window = {WINDOW} bases)")
    ax1.set_ylabel("Value (%)")
    ax1.set_title("C+G% and IC along the sequence")
    ax1.legend()
    ax1.grid(alpha=0.3)

    # right: the pattern as a cloud of points + its center of weight
    ax2.plot(cg, ic, color="gray", alpha=0.4, linewidth=0.8)
    ax2.scatter(cg, ic, c=list(positions), cmap="viridis", s=25, label="windows")
    ax2.scatter(*center, color="red", marker="X", s=200, label=f"center ({center[0]:.2f}, {center[1]:.2f})")
    ax2.set_xlabel("C+G %")
    ax2.set_ylabel("IC")
    ax2.set_title("Pattern (each point = one window)")
    ax2.legend()
    ax2.grid(alpha=0.3)

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    print(f"Saved {filename}")
    if show:
        plt.show()
    plt.close(fig)


# ---------- Step 7: centers on a 2nd chart ----------
def plot_centers(centers, filename, title="Centers of weight of the patterns", show=True):
    """centers = list of (name, (cg, ic))"""
    fig, ax = plt.subplots(figsize=(8, 6))
    for name, (cg, ic) in centers:
        ax.scatter(cg, ic, s=90)
        ax.annotate(name, (cg, ic), textcoords="offset points", xytext=(6, 6), fontsize=9)
    ax.set_xlabel("C+G % (center)")
    ax.set_ylabel("IC (center)")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    print(f"Saved {filename}")
    if show:
        plt.show()
    plt.close(fig)


def main():
    # save all output files next to this script, wherever it is run from
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    seq = sys.argv[1].upper() if len(sys.argv) > 1 else S

    print(f"Sequence ({len(seq)} bases): {seq}\n")
    print("Whole sequence:")
    print(f"  C+G % = {cg_content(seq):.2f}   (board: ~29.27)")
    print(f"  IC    = {information_content(seq):.2f}   (board: ~27.53)\n")

    pattern = make_pattern(seq)
    lines = [f"{'Win':>4} {'Start':>6}  {'C+G%':>6}  {'IC':>6}  Window"]
    for i, ((cg, ic), w) in enumerate(zip(pattern, sliding_windows(seq)), 1):
        lines.append(f"{i:>4} {i:>6}  {cg:6.2f}  {ic:6.2f}  {w}")
    center = center_of_weight(pattern)
    lines.append("")
    lines.append(f"Number of windows: {len(pattern)}")
    lines.append(f"Center of weight: C+G % = {center[0]:.2f}, IC = {center[1]:.2f}")
    text = "\n".join(lines)
    print(text)

    with open("pattern_results.txt", "w") as f:
        f.write(f"Sequence ({len(seq)} bases): {seq}\n")
        f.write(f"Whole sequence: C+G % = {cg_content(seq):.2f}, IC = {information_content(seq):.2f}\n\n")
        f.write(text + "\n")
    print("Saved pattern_results.txt")

    plot_pattern(pattern, "Pattern of the test sequence S", "pattern_S.png")
    plot_centers([("S", center)], "center_S.png", "Center of weight of the pattern of S")


if __name__ == "__main__":
    main()
