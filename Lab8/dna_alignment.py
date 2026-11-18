
import os
import sys

import matplotlib.pyplot as plt

MATCH, MISMATCH, GAP = 1, -1, -2


def score_matrix(a, b, match=MATCH, mismatch=MISMATCH, gap=GAP, local=False):
    rows, cols = len(a) + 1, len(b) + 1
    m = [[0] * cols for _ in range(rows)]
    if not local:
        for i in range(1, rows):
            m[i][0] = i * gap
        for j in range(1, cols):
            m[0][j] = j * gap
    for i in range(1, rows):
        for j in range(1, cols):
            diag = m[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            up = m[i - 1][j] + gap
            left = m[i][j - 1] + gap
            m[i][j] = max(diag, up, left, 0) if local else max(diag, up, left)
    return m


def traceback(a, b, m, match=MATCH, mismatch=MISMATCH, gap=GAP, local=False):
    """Returns (aligned_a, aligned_b, path) where path = list of matrix cells visited."""
    if local:   # start from the highest score in the matrix
        i, j = max(((i, j) for i in range(len(m)) for j in range(len(m[0]))),
                   key=lambda c: m[c[0]][c[1]])
    else:       # start from the bottom-right corner
        i, j = len(a), len(b)

    out_a, out_b, path = [], [], [(i, j)]
    while i > 0 or j > 0:
        if local and m[i][j] == 0:
            break
        if i > 0 and j > 0 and m[i][j] == m[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch):
            out_a.append(a[i - 1]); out_b.append(b[j - 1]); i, j = i - 1, j - 1
        elif i > 0 and m[i][j] == m[i - 1][j] + gap:
            out_a.append(a[i - 1]); out_b.append("-"); i -= 1
        else:
            out_a.append("-"); out_b.append(b[j - 1]); j -= 1
        path.append((i, j))
    return "".join(reversed(out_a)), "".join(reversed(out_b)), path[::-1]


def align(a, b, match=MATCH, mismatch=MISMATCH, gap=GAP, local=False):
    a, b = clean(a), clean(b)
    m = score_matrix(a, b, match, mismatch, gap, local)
    al_a, al_b, path = traceback(a, b, m, match, mismatch, gap, local)
    matches = sum(1 for x, y in zip(al_a, al_b) if x == y)
    gaps = al_a.count("-") + al_b.count("-")
    score = max(max(r) for r in m) if local else m[-1][-1]
    return {
        "a": a, "b": b, "matrix": m, "path": path,
        "aligned_a": al_a, "aligned_b": al_b,
        "middle": "".join("|" if x == y else (" " if "-" in (x, y) else ".") for x, y in zip(al_a, al_b)),
        "score": score, "matches": matches, "gaps": gaps,
        "mismatches": len(al_a) - matches - gaps, "length": len(al_a),
        "identity": matches / len(al_a) * 100 if al_a else 0,
    }


def clean(seq):
    """Keep only A, C, G, T (removes spaces, new lines, FASTA header)."""
    lines = [l for l in seq.strip().splitlines() if not l.startswith(">")]
    return "".join(c for c in "".join(lines).upper() if c in "ACGT")


def report(r, local=False, width=60):
    lines = [f"{'Local (Smith-Waterman)' if local else 'Global (Needleman-Wunsch)'} alignment",
             f"Sequence 1: {len(r['a'])} bases   Sequence 2: {len(r['b'])} bases", ""]
    for k in range(0, r["length"], width):
        lines.append(f"Seq1  {r['aligned_a'][k:k + width]}")
        lines.append(f"      {r['middle'][k:k + width]}")
        lines.append(f"Seq2  {r['aligned_b'][k:k + width]}")
        lines.append("")
    lines += [f"Score: {r['score']}",
              f"Alignment length: {r['length']}",
              f"Matches: {r['matches']}   Mismatches: {r['mismatches']}   Gaps: {r['gaps']}",
              f"Identity: {r['identity']:.2f} %"]
    return "\n".join(lines)


def plot_matrix(r, filename="alignment_matrix.png", show=True):
    """Heat map of the score matrix with the traceback path."""
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.imshow(r["matrix"], cmap="viridis", aspect="auto")
    ys, xs = zip(*r["path"])
    ax.plot(xs, ys, color="red", linewidth=2, label="traceback path")
    if len(r["a"]) <= 40 and len(r["b"]) <= 40:
        ax.set_xticks(range(len(r["b"]) + 1)); ax.set_xticklabels(["-"] + list(r["b"]))
        ax.set_yticks(range(len(r["a"]) + 1)); ax.set_yticklabels(["-"] + list(r["a"]))
    ax.set_xlabel("Sequence 2"); ax.set_ylabel("Sequence 1")
    ax.set_title(f"Score matrix (score = {r['score']}, identity = {r['identity']:.1f}%)")
    ax.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    if show:
        plt.show()
    plt.close(fig)


# ---------------- the app window ----------------
def run_app():
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox

    root = tk.Tk()
    root.title("DNA sequence alignment")
    root.geometry("900x700")

    def seq_box(label, row):
        ttk.Label(root, text=label).grid(row=row, column=0, sticky="nw", padx=8, pady=4)
        box = tk.Text(root, height=4, width=80, font=("Consolas", 10))
        box.grid(row=row, column=1, columnspan=6, sticky="ew", padx=8, pady=4)
        ttk.Button(root, text="Load FASTA", command=lambda: load(box)).grid(row=row, column=7, padx=8)
        return box

    def load(box):
        path = filedialog.askopenfilename(filetypes=[("FASTA", "*.fasta *.fa *.fna *.txt"), ("All", "*.*")])
        if path:
            box.delete("1.0", "end")
            box.insert("1.0", open(path).read())

    box1 = seq_box("Sequence 1", 0)
    box2 = seq_box("Sequence 2", 1)
    box1.insert("1.0", "GATTACAGATTACA")
    box2.insert("1.0", "GCATGCAGATACA")

    params = {}
    for col, (name, value) in enumerate([("Match", MATCH), ("Mismatch", MISMATCH), ("Gap", GAP)]):
        ttk.Label(root, text=name).grid(row=2, column=2 * col, sticky="e", padx=4)
        var = tk.StringVar(value=str(value))
        ttk.Entry(root, textvariable=var, width=6).grid(row=2, column=2 * col + 1, sticky="w")
        params[name] = var
    mode = tk.StringVar(value="global")
    ttk.Radiobutton(root, text="Global", variable=mode, value="global").grid(row=2, column=6)
    ttk.Radiobutton(root, text="Local", variable=mode, value="local").grid(row=2, column=7)

    out = tk.Text(root, font=("Consolas", 10), wrap="none")
    out.grid(row=4, column=0, columnspan=8, sticky="nsew", padx=8, pady=8)
    root.grid_rowconfigure(4, weight=1)
    root.grid_columnconfigure(1, weight=1)
    last = {}

    def do_align():
        try:
            m, mm, g = (int(params[k].get()) for k in ("Match", "Mismatch", "Gap"))
        except ValueError:
            messagebox.showerror("Error", "Match, mismatch and gap must be whole numbers.")
            return
        a, b = box1.get("1.0", "end"), box2.get("1.0", "end")
        if not clean(a) or not clean(b):
            messagebox.showerror("Error", "Both sequences must contain A, C, G or T.")
            return
        if len(clean(a)) * len(clean(b)) > 25_000_000:
            messagebox.showerror("Error", "The sequences are too long (max ~5000 x 5000 bases).")
            return
        local = mode.get() == "local"
        r = align(a, b, m, mm, g, local)
        last["r"] = r
        text = report(r, local)
        out.delete("1.0", "end")
        out.insert("1.0", text)
        with open("alignment_result.txt", "w") as f:
            f.write(text + "\n")

    def show_matrix():
        if "r" not in last:
            do_align()
        if "r" in last:
            plot_matrix(last["r"])

    buttons = ttk.Frame(root)
    buttons.grid(row=3, column=0, columnspan=8, pady=4)
    ttk.Button(buttons, text="Align", command=do_align).pack(side="left", padx=6)
    ttk.Button(buttons, text="Show score matrix", command=show_matrix).pack(side="left", padx=6)

    root.mainloop()


def main():
    # save all output files next to this script, wherever it is run from
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    if len(sys.argv) >= 3:
        local = "--local" in sys.argv
        a, b = [x for x in sys.argv[1:] if x != "--local"][:2]
        if os.path.exists(a): a = open(a).read()
        if os.path.exists(b): b = open(b).read()
        r = align(a, b, local=local)
        text = report(r, local)
        print(text)
        with open("alignment_result.txt", "w") as f:
            f.write(text + "\n")
        plot_matrix(r)
    else:
        run_app()


if __name__ == "__main__":
    main()
