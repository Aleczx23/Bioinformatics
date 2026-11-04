
import random

BASES = "ACGT"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

IR_MIN = 8
IR_MAX = 12
TSD_MIN = 4
TSD_MAX = 6
TE_MIN = 30
TE_MAX = 200


def random_dna(length):
    return "".join(random.choice(BASES) for _ in range(length))


def reverse_complement(seq):
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def make_transposon(ir_len, body_len):
    ir = random_dna(ir_len)
    while ir == reverse_complement(ir):
        ir = random_dna(ir_len)
    return ir + random_dna(body_len) + reverse_complement(ir)


def insert(seq, elements, name, te, pos, tsd_len=5):
    target = seq[pos:pos + tsd_len]
    new_seq = seq[:pos + tsd_len] + te + seq[pos:]
    shift = len(te) + tsd_len
    updated = []
    for n, s, e in elements:
        if s >= pos + tsd_len:
            s += shift
            e += shift
        elif e >= pos + tsd_len:
            e += shift
        updated.append((n, s, e))
    start = pos + tsd_len
    updated.append((name, start, start + len(te) - 1))
    return new_seq, updated, target


def overwrite(seq, elements, name, te, pos, tsd_len=5):
    target = random_dna(tsd_len)
    block = target + te + target
    new_seq = seq[:pos] + block + seq[pos + len(block):]
    start = pos + tsd_len
    elements = elements + [(name, start, start + len(te) - 1)]
    return new_seq, elements, target


def find_transposons(seq):
    hits = {}
    n = len(seq)
    for k in range(IR_MAX, IR_MIN - 1, -1):
        for i in range(TSD_MIN, n - k):
            left = seq[i:i + k]
            right = reverse_complement(left)
            if left == right:
                continue
            for j in range(i + TE_MIN - k, min(i + TE_MAX - k, n - k) + 1):
                if seq[j:j + k] != right:
                    continue
                end = j + k - 1
                tsd = 0
                for t in range(TSD_MAX, TSD_MIN - 1, -1):
                    if i - t >= 0 and end + 1 + t <= n and seq[i - t:i] == seq[end + 1:end + 1 + t]:
                        tsd = t
                        break
                if tsd and (i, end) not in hits:
                    hits[(i, end)] = (k, seq[i - tsd:i])
    return sorted((s, e, k, tsd) for (s, e), (k, tsd) in hits.items())


def build_separate():
    seq = random_dna(120)
    elements = []
    for name, pos in [("T1", 100), ("T2", 75), ("T3", 45), ("T4", 15)]:
        te = make_transposon(random.randint(IR_MIN, IR_MAX), random.randint(25, 40))
        seq, elements, _ = insert(seq, elements, name, te, pos)
    return seq, elements


def build_two_overlapping():
    seq = random_dna(260)
    elements = []
    t1 = make_transposon(10, 50)
    seq, elements, _ = insert(seq, elements, "T1", t1, 40)
    t2 = make_transposon(10, 40)
    t1_end = elements[0][2]
    seq, elements, _ = overwrite(seq, elements, "T2", t2, t1_end - 30)
    t3 = make_transposon(9, 30)
    seq, elements, _ = insert(seq, elements, "T3", t3, len(seq) - 20)
    return seq, elements


def build_nested():
    seq = random_dna(200)
    elements = []
    t1 = make_transposon(10, 60)
    seq, elements, _ = insert(seq, elements, "T1", t1, 50)
    t2 = make_transposon(9, 30)
    t1_start = elements[0][1]
    seq, elements, _ = insert(seq, elements, "T2", t2, t1_start + 30)
    t3 = make_transposon(10, 35)
    seq, elements, _ = insert(seq, elements, "T3", t3, 10)
    return seq, elements


def build_three_intersecting():
    seq = random_dna(200)
    elements = []
    t1 = make_transposon(11, 60)
    seq, elements, _ = insert(seq, elements, "T1", t1, 30)
    t1_start, t1_end = elements[0][1], elements[0][2]
    t2 = make_transposon(10, 45)
    seq, elements, _ = overwrite(seq, elements, "T2", t2, t1_end - 25)
    t2_start = elements[1][1]
    t3 = make_transposon(9, 25)
    seq, elements, _ = insert(seq, elements, "T3", t3, t2_start + 15)
    t4 = make_transposon(10, 30)
    seq, elements, _ = insert(seq, elements, "T4", t4, t1_start + 3)
    return seq, elements


def write_fasta(path, name, seq):
    with open(path, "w") as f:
        f.write(f">{name} length={len(seq)}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")


def report(title, seq, elements, filename):
    write_fasta(filename, title.replace(" ", "_"), seq)
    found = find_transposons(seq)
    print("=" * 60)
    print(title)
    print("Sequence length:", len(seq))
    print("Saved to:", filename)
    print()
    print("Inserted transposons:")
    for name, s, e in sorted(elements, key=lambda x: x[1]):
        print(f"  {name}: start {s + 1}, end {e + 1}")
    print()
    print("Detected transposons:")
    if not found:
        print("  none")
    for s, e, k, tsd in found:
        match = [n for n, ts, te in elements if ts == s and te == e]
        label = match[0] if match else "false positive"
        print(f"  start {s + 1}, end {e + 1}, IR {k} bp {seq[s:s + k]}, direct repeat {tsd} -> {label}")
    detected = {(s, e) for s, e, _, _ in found}
    missed = [n for n, s, e in elements if (s, e) not in detected]
    print()
    print("Missed:", ", ".join(missed) if missed else "none")
    print()


def main():
    random.seed(6)
    scenarios = [
        ("Case 1 separate transposons", build_separate, "case1_separate.fasta"),
        ("Case 2 two overlapping transposons", build_two_overlapping, "case2_overlap.fasta"),
        ("Case 3 nested transposons", build_nested, "case3_nested.fasta"),
        ("Case 4 three intersecting transposons", build_three_intersecting, "case4_three_intersecting.fasta"),
    ]
    for title, builder, filename in scenarios:
        seq, elements = builder()
        report(title, seq, elements, filename)


if __name__ == "__main__":
    main()
