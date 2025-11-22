import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import math

DATA_DIR = Path("data")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

cons_file = DATA_DIR / "VANGL1_conservation.tsv"
vars_file = DATA_DIR / "VANGL1_variants.tsv"
domains_file = DATA_DIR / "VANGL1_domains.tsv"

# -------------------------
# Load data
# -------------------------
cons_df = pd.read_csv(cons_file, sep="\t")
pos_col = cons_df.columns[0]
score_col = cons_df.columns[1]
cons_df = cons_df.rename(columns={pos_col: "pos", score_col: "score"})
cons_df["pos"] = cons_df["pos"].astype(int)

vars_df = pd.read_csv(vars_file, sep="\t")
vars_df["pos"] = vars_df["pos"].astype(int)

domains_df = pd.read_csv(domains_file, sep="\t")

sites = vars_df[["pos", "group"]].drop_duplicates()

# -------------------------
# Define categories
# -------------------------
q90 = cons_df["score"].quantile(0.90)
q95 = cons_df["score"].quantile(0.95)
q80 = cons_df["score"].quantile(0.80)

top10_pos = set(cons_df.loc[cons_df["score"] >= q90, "pos"])
top5_pos  = set(cons_df.loc[cons_df["score"] >= q95, "pos"])
top20_pos = set(cons_df.loc[cons_df["score"] >= q80, "pos"])

pdz_rows = domains_df[domains_df["name"].str.contains("PDZ", case=False, na=False)]
pdz_pos = set()
for _, r in pdz_rows.iterrows():
    pdz_pos.update(range(int(r["start"]), int(r["end"]) + 1))

categories = [
    ("Top 10% constraint", top10_pos),
    ("Top 5% constraint",  top5_pos),
    ("Top 20% constraint", top20_pos),
    ("PDZ motif",          pdz_pos),
]

def compute_or_ci(cat_pos):
    a = ((sites["group"] == "NTD") & (sites["pos"].isin(cat_pos))).sum()
    b = ((sites["group"] == "NTD") & (~sites["pos"].isin(cat_pos))).sum()
    c = ((sites["group"] == "benign") & (sites["pos"].isin(cat_pos))).sum()
    d = ((sites["group"] == "benign") & (~sites["pos"].isin(cat_pos))).sum()

    a1, b1, c1, d1 = a + 0.5, b + 0.5, c + 0.5, d + 0.5

    or_val = (a1 * d1) / (b1 * c1)
    log_or = math.log(or_val)
    se = math.sqrt(1/a1 + 1/b1 + 1/c1 + 1/d1)
    ci_low = math.exp(log_or - 1.96 * se)
    ci_high = math.exp(log_or + 1.96 * se)
    return or_val, ci_low, ci_high, (a, b, c, d)

results = []
for name, pos_set in categories:
    or_val, ci_low, ci_high, counts = compute_or_ci(pos_set)
    results.append((name, or_val, ci_low, ci_high, counts))

# -------------------------
# Forest plot
# -------------------------
plt.style.use("default")
fig, ax = plt.subplots(figsize=(6,4))

y = np.arange(len(results))[::-1]
or_vals = [r[1] for r in results]
ci_lows = [r[2] for r in results]
ci_highs = [r[3] for r in results]

xmin = min(ci_lows) * 0.6
xmax = max(ci_highs) * 1.6

colors = ["tab:blue", "tab:blue", "tab:blue", "tab:orange"]

for yi, or_v, lo, hi, col in zip(y, or_vals, ci_lows, ci_highs, colors):
    ax.plot([lo, hi], [yi, yi], color=col, linewidth=1.8)
    ax.scatter(or_v, yi, color=col, s=40, zorder=3)

ax.axvline(1.0, color="gray", linestyle="--", linewidth=1)

ax.set_yticks(y)
ax.set_yticklabels([r[0] for r in results])

ax.set_xscale("log")
ax.set_xlim(xmin, xmax)
ax.set_xlabel("Odds ratio (NTD vs benign)")
ax.set_title("NTD variant enrichment in key categories")

ax.grid(axis="x", linestyle="--", alpha=0.3)
plt.tight_layout(rect=[0.18, 0.08, 0.98, 0.95])

outfile = OUTPUT_DIR / "Figure3_enrichment_forest.png"
plt.savefig(outfile, dpi=300)
print("Saved:", outfile)

print("\nCategory\tOR\tCI_low\tCI_high\ta\tb\tc\td")
for name, or_v, lo, hi, (a,b,c,d) in results:
    print(f"{name}\t{or_v:.2f}\t{lo:.2f}\t{hi:.2f}\t{a}\t{b}\t{c}\t{d}")
