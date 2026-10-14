
GENETIC_CODE = {
    # U _ _
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",
    # C _ _
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    # A _ _
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    # G _ _
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

START_CODON = "AUG"


def clean_sequence(seq):
    """Uppercase, remove spaces/newlines, and convert DNA (T) to ARN (U)."""
    seq = "".join(seq.split()).upper().replace("T", "U")
    invalid = set(seq) - set("UCAG")
    if invalid:
        raise ValueError(f"Invalid letters in sequence: {', '.join(sorted(invalid))}")
    return seq


def translate(arn):
    """
    Translate an ARN sequence into a protein.
    Returns (protein_list, codons_read, start_index, stopped_properly).
    """
    arn = clean_sequence(arn)

    start = arn.find(START_CODON)
    if start == -1:
        return [], [], -1, False  # no start codon -> no protein

    protein = []
    codons = []
    for i in range(start, len(arn) - 2, 3):
        codon = arn[i:i + 3]
        amino_acid = GENETIC_CODE[codon]
        codons.append(codon)
        if amino_acid == "Stop":
            return protein, codons, start, True
        protein.append(amino_acid)

    # reached the end of the sequence without a stop codon
    return protein, codons, start, False


def show_result(arn):
    protein, codons, start, stopped = translate(arn)
    print(f"\nARN sequence : {clean_sequence(arn)}")

    if start == -1:
        print("No START codon (AUG) found -> no protein is made.")
        return

    print(f"START codon found at position {start + 1}")
    print("Codons read  : " + " ".join(codons))
    print("Protein      : " + "-".join(protein))
    print(f"Length       : {len(protein)} amino acids")
    if stopped:
        print(f"Ended at STOP codon {codons[-1]}")
    else:
        print("Warning: no STOP codon found - the chain ended with the sequence.")


if __name__ == "__main__":
    # Example from the exercise style
    example = "GCAUGUUUCCAGGAUAAGCU"
    print("=== Example ===")
    show_result(example)

    print("\n=== Your turn ===")
    while True:
        user_seq = input("\nEnter an ARN sequence (or press Enter to quit): ")
        if not user_seq.strip():
            break
        try:
            show_result(user_seq)
        except ValueError as e:
            print(e)
