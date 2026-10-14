

ALPHABET = ("A", "C", "G", "T")


def validate_sequence(sequence: str) -> str:
    if not sequence:
        raise ValueError("The DNA sequence is empty.")
    invalid = sorted(set(sequence) - set(ALPHABET))
    if invalid:
        raise ValueError(f"Invalid characters in sequence: {', '.join(invalid)}")
    return sequence


def count_nucleotides(sequence: str) -> dict:
    """Return the absolute count of each nucleotide."""
    counts = {letter: 0 for letter in ALPHABET}
    for letter in sequence:
        counts[letter] += 1
    return counts


def relative_frequencies(counts: dict) -> dict:
    """Return the relative frequency (0..1) of each nucleotide."""
    total = sum(counts.values())
    return {letter: counts[letter] / total for letter in ALPHABET}


def print_report(sequence: str) -> None:
    sequence = validate_sequence(sequence)
    counts = count_nucleotides(sequence)
    freqs = relative_frequencies(counts)

    print(f"Sequence S = {sequence}")
    print(f"Length     = {len(sequence)}\n")
    print(f"{'Letter':<8}{'Count':>6}{'Frequency':>12}{'Percent':>10}")
    print("-" * 36)
    for letter in ALPHABET:
        print(f"{letter:<8}{counts[letter]:>6}{freqs[letter]:>12.2f}{freqs[letter] * 100:>9.1f}%")
    print("-" * 36)
    print(f"{'Total':<8}{len(sequence):>6}{1:>12.2f}{100:>9.1f}%")


if __name__ == "__main__":
    S = "ACGTTGCATTACGTTGCAGT"
    print_report(S)