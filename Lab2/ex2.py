import math


def count_bases(sequence):
    return (sequence.count("A"), sequence.count("T"),
            sequence.count("C"), sequence.count("G"))


def validate(sequence):
    sequence = sequence.upper()
    if not sequence or any(base not in "ATCG" for base in sequence):
        raise ValueError("Sequence must contain only A, T, C, G")
    return sequence


def tm_basic(sequence):
    a, t, c, g = count_bases(sequence)
    return 4 * (c + g) + 2 * (a + t)


def gc_percent(sequence):
    a, t, c, g = count_bases(sequence)
    return (c + g) / len(sequence) * 100


def tm_salt(sequence, na_concentration):
    if na_concentration <= 0:
        raise ValueError("Na+ concentration must be greater than 0")
    gc = gc_percent(sequence)
    return 81.5 + 16.6 * math.log10(na_concentration) + 0.41 * gc - 600 / len(sequence)


def main():
    sequence = input("Enter DNA sequence (default ACGCGTGCCA): ").strip() or "ACGCGTGCCA"
    na_input = input("Enter Na+ concentration in M (default 0.05): ").strip() or "0.05"
    try:
        sequence = validate(sequence)
        na = float(na_input)
        a, t, c, g = count_bases(sequence)
        basic = tm_basic(sequence)
        salt = tm_salt(sequence, na)
    except ValueError as error:
        print(f"Error: {error}")
        return
    print(f"\nSequence: {sequence}")
    print(f"Length: {len(sequence)}")
    print(f"A = {a}, T = {t}, C = {c}, G = {g}")
    print(f"GC content: {gc_percent(sequence):.1f}%")
    print(f"\nFormula 1: Tm = 4*(C+G) + 2*(A+T)")
    print(f"Tm = 4*({c}+{g}) + 2*({a}+{t}) = {basic} °C")
    print(f"\nFormula 2: Tm = 81.5 + 16.6*log10[Na+] + 0.41*(%GC) - 600/length")
    print(f"Tm = 81.5 + 16.6*log10({na}) + 0.41*{gc_percent(sequence):.1f} - 600/{len(sequence)} = {salt:.2f} °C")


if __name__ == "__main__":
    main()
    