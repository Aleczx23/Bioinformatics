"""
Exercise 2

Take an arbitrary DNA sequence from NCBI, between 1000 and 3000 nucleotides. Implement a
software application that detects repetitions between 6 and 10 bases in the DNA sequence.
The length of the repeated unit is between 3 and 5 bases.

Sequence used: the same Cloning vector pUC19c (NCBI L09137.2, 2,686 bp) from Exercise 1.
"""


def read_fasta(path):
    sequence = ""
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip().upper()
    return sequence


def is_primitive(unit):
    for p in range(1, len(unit)):
        if len(unit) % p == 0 and unit[:p] * (len(unit) // p) == unit:
            return False
    return True


def find_repeats(sequence):
    repeats = []
    for k in range(3, 6):
        i = 0
        while i + k < len(sequence):
            j = i
            while j + k < len(sequence) and sequence[j] == sequence[j + k]:
                j += 1
            copies = (j - i + k) // k
            unit = sequence[i:i + k]
            total = copies * k
            if copies >= 2 and 6 <= total <= 10 and is_primitive(unit):
                repeats.append((i + 1, unit, copies, sequence[i:i + total]))
            i = j + 1 if j > i else i + 1
    return sorted(repeats)


def main():
    sequence = read_fasta("pUC19c_L09137.2.fasta")
    repeats = find_repeats(sequence)

    for position, unit, copies, text in repeats:
        print(f"position {position}: {unit} x {copies} = {text}")

    print("Total repetitions found:", len(repeats))


if __name__ == "__main__":
    main()
