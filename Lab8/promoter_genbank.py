"""
Lab8 - Step 8: take real promoter sequences from GenBank and generate their patterns

For each gene:
  - find the gene in NCBI (Homo sapiens) and its position on the chromosome
  - download the promoter = the 1000 bases upstream of the transcription start,
    on the gene's own strand, and save it as a FASTA file next to this script
  - make the pattern (C+G% and IC for every 30-base window) and plot it
  - compute the center of weight of the pattern
Then all the centers (plus the test sequence S) are plotted on a 2nd chart.

Run (needs internet the first time; afterwards the saved FASTA files are used):
    py promoter_genbank.py                      (TP53, GAPDH, ACTB, INS, BRCA1)
    py promoter_genbank.py MYC HBB EGFR         (any human gene symbols)
    py promoter_genbank.py my_promoter.fasta    (a FASTA you downloaded yourself)
"""

import json
import os
import sys
import time
import urllib.parse
import urllib.request

from promoter_pattern import (S, WINDOW, cg_content, information_content, make_pattern,
                              center_of_weight, plot_pattern, plot_centers)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PROMOTER_LENGTH = 1000
DEFAULT_GENES = ["TP53", "GAPDH", "ACTB", "INS", "BRCA1"]


def ncbi(tool, **params):
    params["tool"] = "bioinformatics_lab8"
    url = EUTILS + tool + ".fcgi?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as r:
        data = r.read().decode()
    time.sleep(0.4)  # NCBI allows max 3 requests per second
    return data


def read_fasta(path):
    header, seq = "", []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            elif line:
                seq.append(line.upper())
    return header, "".join(seq)


def download_promoter(gene):
    """Download the 1000 bp upstream of the gene's start (human), save it as FASTA."""
    filename = f"{gene}_promoter.fasta"
    if os.path.exists(filename):
        return filename

    ids = json.loads(ncbi("esearch", db="gene", retmode="json",
                          term=f"{gene}[sym] AND Homo sapiens[orgn]"))["esearchresult"]["idlist"]
    if not ids:
        raise ValueError(f"Gene {gene} not found in NCBI Gene")
    summary = json.loads(ncbi("esummary", db="gene", id=ids[0], retmode="json"))["result"][ids[0]]
    info = summary["genomicinfo"][0]
    chrom, start, stop = info["chraccver"], int(info["chrstart"]), int(info["chrstop"])

    # NCBI gives 0-based coordinates; start is the transcription start of the gene
    if start <= stop:   # gene on the + strand -> promoter is to the left
        seq_start, seq_stop, strand = start + 1 - PROMOTER_LENGTH, start, 1
    else:               # gene on the - strand -> promoter is to the right, reverse complement
        seq_start, seq_stop, strand = start + 2, start + 1 + PROMOTER_LENGTH, 2

    fasta = ncbi("efetch", db="nuccore", id=chrom, rettype="fasta", retmode="text",
                 seq_start=seq_start, seq_stop=seq_stop, strand=strand)
    with open(filename, "w") as f:
        f.write(fasta)
    print(f"Downloaded {gene} promoter: {chrom}:{seq_start}-{seq_stop} "
          f"({'+' if strand == 1 else '-'} strand) -> {filename}")
    return filename


def main():
    # save all output files next to this script, wherever it is run from
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    inputs = sys.argv[1:] or DEFAULT_GENES
    centers = [("S (test)", center_of_weight(make_pattern(S)))]
    lines = [f"Window = {WINDOW} bases, promoter = {PROMOTER_LENGTH} bases upstream of the gene start", ""]

    for item in inputs:
        try:
            path = item if item.lower().endswith((".fasta", ".fa", ".fna")) else download_promoter(item)
        except Exception as e:
            print(f"Could not get {item}: {e}")
            continue
        name = os.path.splitext(os.path.basename(path))[0].replace("_promoter", "")
        header, seq = read_fasta(path)

        pattern = make_pattern(seq)
        center = center_of_weight(pattern)
        centers.append((name, center))

        lines.append(f"{name}: {header}")
        lines.append(f"  length: {len(seq)} bases, windows: {len(pattern)}")
        lines.append(f"  whole sequence: C+G % = {cg_content(seq):.2f}, IC = {information_content(seq):.2f}")
        lines.append(f"  C+G % range in windows: {min(p[0] for p in pattern):.2f} - {max(p[0] for p in pattern):.2f}")
        lines.append(f"  IC range in windows:    {min(p[1] for p in pattern):.2f} - {max(p[1] for p in pattern):.2f}")
        lines.append(f"  center of weight: C+G % = {center[0]:.2f}, IC = {center[1]:.2f}")
        lines.append("")

        plot_pattern(pattern, f"Pattern of the {name} promoter", f"pattern_{name}.png", show=False)

    text = "\n".join(lines)
    print(text)
    with open("promoter_results.txt", "w") as f:
        f.write(text + "\n")
    print("Saved promoter_results.txt")

    plot_centers(centers, "centers_promoters.png", "Centers of weight of the promoter patterns")


if __name__ == "__main__":
    main()
