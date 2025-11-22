import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from collections import Counter

DATA_DIR = Path("data")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

fasta_file = DATA_DIR / "VANGL1_aligned.fasta"

# -------------------------
# Simple FASTA parser
# -------------------------
seqs = []
current = []
with open(fasta_file) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                seqs.append("".join(current))
                current = []
        else:
            current.append(line)
    if current:
        seqs.append("".join(current))

if not seqs:
    raise ValueError("No sequences found in FASTA")

L = len(seqs[0])
for s in seqs:
    if len(s) != L:
        raise ValueError("Sequences are not all the same length, alignment may be broken")

# -------------------------
# Take last ~40 alignment columns
# -------------------------
window = 40
start = max(0, L - window)
end = L

sub_seqs = [s[start:end] for s in seqs]
sub_L = end - start

# -------------------------
# Compute consensus and fraction-of-consensus per column
# -------------------------
consensus = []
frac_cons = []

for i in range(sub_L):
    col = [s[i] for s in sub_seqs]
    nongap = [aa for aa in col if aa not in ("-", ".")]
    if not nongap:
        consensus.append("-")
        frac_cons.append(0.0)
        continue
    from collections import Counter
    counts = Counter(nongap)
    aa, count = counts.most_common(1)[0]
    consensus.append(aa)
    frac_cons.append(count / len(nongap))

consensus = np.array(consensus)
frac_cons = np.array(frac_cons)

# Define PDZ motif as last 4 positions in this 40-aa window
pdz_len = 4
pdz_start_idx = max(0, sub_L - pdz_len)
pdz_mask = np.zeros(sub_L, dtype=bool)
pdz_mask[pdz_start_idx:] = True

# -------------------------
# Plot
# -------------------------
plt.style.use("default")
fig, ax = plt.subplots(figsize=(8,3))

x = np.arange(1, sub_L + 1)

ax.bar(x[~pdz_mask], frac_cons[~pdz_mask],
       color="tab:blue", alpha=0.8, label="C-terminal region")
ax.bar(x[pdz_mask], frac_cons[pdz_mask],
       color="tab:orange", alpha=0.9, label="PDZ-binding motif")

ax.set_xticks(x)
ax.set_xticklabels(consensus, fontsize=8)
ax.set_xlim(0.5, sub_L + 0.5)
ax.set_ylim(0, 1.05)

ax.set_ylabel("Fraction consensus")
ax.set_xlabel("Aligned C-terminal positions (last ~40 aa)")
ax.set_title("C-terminal VANGL1 alignment highlighting PDZ motif")
ax.legend(frameon=True, fontsize=8)

ax.grid(axis="y", linestyle="--", alpha=0.4)
plt.tight_layout()

outfile = OUTPUT_DIR / "Figure4_C_terminal_PDZ_region.png"
plt.savefig(outfile, dpi=300)
print("Saved:", outfile)
