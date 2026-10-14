def calculate_tm(sequence):
    sequence = sequence.upper()
    valid = set("ATCG")
    if not sequence or any(base not in valid for base in sequence):
        raise ValueError("Sequence must contain only A, T, C, G")
    a = sequence.count("A")
    t = sequence.count("T")
    c = sequence.count("C")
    g = sequence.count("G")
    tm = 4 * (c + g) + 2 * (a + t)
    return a, t, c, g, tm


def main():
    sequence = input("Enter DNA sequence (default ATCGCGTA): ").strip() or "ATCGCGTA"
    try:
        a, t, c, g, tm = calculate_tm(sequence)
    except ValueError as error:
        print(error)
        return
    print(f"Sequence: {sequence.upper()}")
    print(f"Length: {len(sequence)}")
    print(f"A = {a}, T = {t}, C = {c}, G = {g}")
    gc_content = (c + g) / len(sequence) * 100
    print(f"GC content: {gc_content:.1f}%")
    print(f"Tm = 4*({c}+{g}) + 2*({a}+{t}) = {tm} °C")


if __name__ == "__main__":
    main()