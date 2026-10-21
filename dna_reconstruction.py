import math
import random
import sys

NUM_SAMPLES = 2000
MIN_LEN, MAX_LEN = 100, 120
SEED = 42


def read_fasta(path):
    header, parts = "", []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            elif line:
                parts.append(line.upper())
    return header, "".join(parts)


def sample_reads(genome, n, min_len, max_len, rng):
    reads, positions = [], []
    for _ in range(n):
        length = rng.randint(min_len, max_len)
        start = rng.randint(0, len(genome) - length)
        reads.append(genome[start:start + length])
        positions.append(start)
    return reads, positions


def naive_rebuild(reads):
    return "".join(reads)


def similarity(a, b):
    n = min(len(a), len(b))
    same = sum(1 for i in range(n) if a[i] == b[i])
    return same / max(len(a), len(b))


def coverage_gaps(genome_len, reads, positions):
    covered = [0] * genome_len
    for read, start in zip(reads, positions):
        for i in range(start, start + len(read)):
            covered[i] += 1
    gaps, i = [], 0
    while i < genome_len:
        if covered[i] == 0:
            j = i
            while j < genome_len and covered[j] == 0:
                j += 1
            gaps.append((i, j))
            i = j
        else:
            i += 1
    uncovered = sum(1 for c in covered if c == 0)
    return covered, gaps, uncovered


def remove_contained(reads):
    unique = sorted(set(reads), key=len, reverse=True)
    kept = []
    for r in unique:
        if not any(r in k for k in kept):
            kept.append(r)
    return kept


def greedy_assemble(reads, min_overlap=20):
    reads = remove_contained(reads)
    k = min_overlap
    by_prefix = {}
    for idx, r in enumerate(reads):
        by_prefix.setdefault(r[:k], []).append(idx)

    overlaps = []
    for a, ra in enumerate(reads):
        for i in range(1, len(ra) - k + 1):
            for b in by_prefix.get(ra[i:i + k], []):
                if b != a and reads[b].startswith(ra[i:]):
                    overlaps.append((len(ra) - i, a, b))
    overlaps.sort(reverse=True)

    nxt, prv, olen = {}, {}, {}
    for length, a, b in overlaps:
        if a in nxt or b in prv:
            continue
        x = b
        while x in nxt:
            x = nxt[x]
        if x == a:
            continue
        nxt[a], prv[b], olen[a] = b, a, length

    contigs = []
    for start in range(len(reads)):
        if start in prv:
            continue
        seq, x = reads[start], start
        while x in nxt:
            seq += reads[nxt[x]][olen[x]:]
            x = nxt[x]
        contigs.append(seq)
    return sorted(contigs, key=len, reverse=True)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "genome.fasta"
    header, genome = read_fasta(path)
    G = len(genome)
    rng = random.Random(SEED)

    print("=" * 70)
    print("Genome:", header)
    print(f"Length: {G:,} bases")

    reads, positions = sample_reads(genome, NUM_SAMPLES, MIN_LEN, MAX_LEN, rng)
    total = sum(len(r) for r in reads)
    print(f"\nSampled {len(reads)} reads, stored in a list.")
    print(f"First read ({len(reads[0])} bp): {reads[0][:60]}...")
    print(f"Total sampled bases: {total:,}  ->  coverage = {total / G:.1f}x")

    naive = naive_rebuild(reads)
    print("\n Naive rebuild (concatenate reads in the order we have them) ")
    print(f"Rebuilt length : {len(naive):,} bases  (real genome: {G:,})")
    print(f"Length error   : {len(naive) / G:.1f}x too long")
    print(f"First {G:,} bases match the real genome at {similarity(naive[:G], genome):.1%}"
          " of positions  (pure chance gives ~25%)")
    print(f"Rebuilt == genome? {naive == genome}")

    print("\n Why it fails ")
    sorted_rebuild = naive_rebuild([r for _, r in sorted(zip(positions, reads))])
    print("1) We don't know where each read comes from (order is lost).")
    print(f"   Even sorting by the TRUE positions and just concatenating gives "
          f"{len(sorted_rebuild):,} bases,")
    print("   because neighbouring reads OVERLAP and the shared bases get repeated.")
    digits = math.lgamma(NUM_SAMPLES + 1) / math.log(10)
    print(f"2) Trying every possible order = {NUM_SAMPLES}! ~ 10^{digits:.0f} orderings"
          " -> impossible by brute force.")
    covered, gaps, uncovered = coverage_gaps(G, reads, positions)
    expected_gaps = NUM_SAMPLES * math.exp(-NUM_SAMPLES * (110 - 20) / G)
    print(f"3) Coverage holes: {uncovered} bases ({uncovered / G:.2%}) were never sampled,")
    print(f"   in {len(gaps)} gaps {gaps[:5]}{' ...' if len(gaps) > 5 else ''}")
    print(f"   -> no method can recover those bases from these reads.")
    print(f"   Lander-Waterman estimate of breaks (min overlap 20): ~{expected_gaps:.0f}")
    print(f"4) Overlap between reads is the only clue. Comparing every pair = "
          f"{NUM_SAMPLES * (NUM_SAMPLES - 1):,} comparisons.")

    print("\n For comparison: greedy overlap assembly ")
    contigs = greedy_assemble(reads, min_overlap=20)
    largest = contigs[0]
    print(f"Contigs produced: {len(contigs)}")
    print(f"Largest contig  : {len(largest):,} bases "
          f"({len(largest) / G:.1%} of the genome)")
    print(f"Largest contig is an exact piece of the real genome? {largest in genome}")
    print(f"Total length of contigs: {sum(map(len, contigs)):,} bases")
    print("=" * 70)


if __name__ == "__main__":
    main()