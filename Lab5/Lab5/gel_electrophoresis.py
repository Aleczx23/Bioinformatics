"""
Exercise 1

Gel electrophoresis is an analysis method used in all disciplines of the life sciences.
The results of gel electrophoresis indicate the relative sizes of the fragments, which is
useful for restriction mapping and for analyzing PCR fragments.

1. Take an arbitrary DNA sequence from NCBI (National Center for Biotechnology Information),
   between 1000 and 3000 nucleotides (letters).
2. Take 10 random samples from this sequence, between 100 and 3000 bases.
3. Store these samples in an array.
4. Simulate the migration of these DNA segments on the electrophoresis gel, based on their
   molecular weight (number of nucleotides).

Note: Short DNA fragments meet small friction forces and travel fast through the
electrophoresis gel, while long DNA fragments exhibit a high friction force and travel
slowly through the gel.

Sequence used: Cloning vector pUC19c, complete sequence (NCBI L09137.2, 2,686 bp).
"""

import math
import random

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

LADDER = [3000, 2500, 2000, 1500, 1000, 750, 500, 250, 100]


def read_fasta(path):
    sequence = ""
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip().upper()
    return sequence


def take_samples(sequence, n=10, min_len=100, max_len=3000):
    max_len = min(max_len, len(sequence))
    samples = []
    for _ in range(n):
        length = random.randint(min_len, max_len)
        start = random.randint(0, len(sequence) - length)
        samples.append(sequence[start:start + length])
    return samples


def distance(length, time=1.0):
    position = (math.log10(3000) - math.log10(length)) / (math.log10(3000) - math.log10(100))
    return (1 + 8 * position) * time


def draw_gel(ax, lengths, time=1.0):
    ax.clear()
    ax.set_facecolor("black")
    ax.set_xlim(-1, len(lengths) + 0.5)
    ax.set_ylim(10, -0.5)
    ax.set_xticks(range(len(lengths) + 1))
    ax.set_xticklabels(["L"] + [str(i) for i in range(1, len(lengths) + 1)])
    ax.set_ylabel("distance (cm)")
    ax.set_title("Gel electrophoresis")

    for size in LADDER:
        y = distance(size, time)
        ax.plot([-0.3, 0.3], [y, y], color="deepskyblue", linewidth=4)
        ax.text(-0.4, y, str(size), color="white", fontsize=7, ha="right", va="center")

    for lane, size in enumerate(lengths, start=1):
        y = distance(size, time)
        ax.plot([lane - 0.3, lane + 0.3], [y, y], color="orange", linewidth=4)
        ax.text(lane, y + 0.3, str(size), color="white", fontsize=7, ha="center")


def main():
    sequence = read_fasta("pUC19c_L09137.2.fasta")
    print("Sequence length:", len(sequence))

    samples = take_samples(sequence)
    lengths = [len(s) for s in samples]

    for i, length in enumerate(lengths, start=1):
        print(f"Sample {i}: {length} bp, distance {distance(length):.2f} cm")

    fig, ax = plt.subplots(figsize=(9, 7))

    for frame in range(1, 51):
        draw_gel(ax, lengths, frame / 50)
        plt.pause(0.05)

    fig.savefig("gel_electrophoresis.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    main()