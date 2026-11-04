import os
from bisect import bisect_left, bisect_right

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

IR_MIN = 10
IR_MAX = 30
BODY_MIN = 10
BODY_MAX = 3000
TSD_MIN = 4
TSD_MAX = 9

FOLDER = os.path.dirname(os.path.abspath(__file__))
GENOME_FILE = os.path.join(FOLDER, "NC_000915.1.fasta")
RESULTS_FILE = os.path.join(FOLDER, "genome_transposons_results.txt")


def read_fasta(path):
    header = ""
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
            else:
                parts.append(line.strip().upper())
    return header, "".join(parts)


def reverse_complement(seq):
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(seq))


def build_index(seq, k):
    index = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if "N" in kmer or any(b not in COMPLEMENT for b in kmer):
            continue
        index.setdefault(kmer, []).append(i)
    return index


def pairs(a, b):
    return COMPLEMENT.get(a) == b


def find_transposons(seq):
    n = len(seq)
    index = build_index(seq, IR_MIN)
    found = []
    for i in range(1, n - IR_MIN):
        left = seq[i:i + IR_MIN]
        if left not in index:
            continue
        right = reverse_complement(left)
        if right == left or right not in index:
            continue
        positions = index[right]
        lo = bisect_left(positions, i + IR_MIN + BODY_MIN)
        hi = bisect_right(positions, i + IR_MIN + BODY_MAX)
        for j in positions[lo:hi]:
            end = j + IR_MIN - 1
            k = IR_MIN
            while k < IR_MAX and (j - (k - IR_MIN)) - (i + k) - 2 >= BODY_MIN and pairs(seq[i + k], seq[j - (k - IR_MIN) - 1]):
                k += 1
            body = (j - (k - IR_MIN)) - (i + k)
            tsd = 0
            for t in range(TSD_MAX, TSD_MIN - 1, -1):
                if i - t >= 0 and end + 1 + t <= n and seq[i - t:i] == seq[end + 1:end + 1 + t]:
                    tsd = t
                    break
            if tsd:
                found.append((i, end, k, body, seq[i - tsd:i]))
    return merge_shifted(found)


def merge_shifted(found):
    kept = []
    for hit in sorted(found, key=lambda h: (h[0] + h[1], h[0])):
        if kept:
            last = kept[-1]
            shift = hit[0] - last[0]
            if last[0] + last[1] == hit[0] + hit[1] and abs(shift) <= TSD_MAX:
                if len(hit[4]) > len(last[4]):
                    kept[-1] = hit
                continue
        kept.append(hit)
    return sorted(kept)


def main():
    header, seq = read_fasta(GENOME_FILE)
    print("Genome:", header)
    print("Length:", len(seq))
    print(f"IR {IR_MIN}-{IR_MAX} bp, body {BODY_MIN}-{BODY_MAX} bp, direct repeats {TSD_MIN}-{TSD_MAX} bp")
    print()

    found = find_transposons(seq)

    with open(RESULTS_FILE, "w") as f:
        f.write(f"Genome: {header}\n")
        f.write(f"Length: {len(seq)}\n")
        f.write(f"IR {IR_MIN}-{IR_MAX} bp, body {BODY_MIN}-{BODY_MAX} bp, direct repeats {TSD_MIN}-{TSD_MAX} bp\n")
        f.write(f"Possible transposons found: {len(found)}\n\n")
        f.write("start\tend\tlength\tIR_length\tbody_length\tdirect_repeat\tleft_IR\tright_IR\n")
        for s, e, k, body, tsd in found:
            f.write(f"{s + 1}\t{e + 1}\t{e - s + 1}\t{k}\t{body}\t{tsd}\t{seq[s:s + k]}\t{seq[e - k + 1:e + 1]}\n")

    for s, e, k, body, tsd in found[:20]:
        print(f"start {s + 1}, end {e + 1}, length {e - s + 1}, IR {k} bp {seq[s:s + k]}, direct repeat {tsd}")
    if len(found) > 20:
        print(f"... {len(found) - 20} more in {os.path.basename(RESULTS_FILE)}")
    print()
    print("Possible transposons found:", len(found))

    nested = 0
    overlapping = 0
    ordered = sorted((s, e) for s, e, _, _, _ in found)
    for a in range(len(ordered)):
        for b in range(a + 1, len(ordered)):
            if ordered[b][0] > ordered[a][1]:
                break
            if ordered[b][1] <= ordered[a][1]:
                nested += 1
            else:
                overlapping += 1
    print("Nested pairs:", nested)
    print("Overlapping pairs:", overlapping)
    print("Results saved to:", os.path.basename(RESULTS_FILE))


if __name__ == "__main__":
    main()
