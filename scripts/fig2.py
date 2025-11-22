import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

DATA_DIR = Path("data")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

cons_file = DATA_DIR / "VANGL1_conservation.tsv"
vars_file = DATA_DIR / "VANGL1_variants.tsv"

# Load conservation
cons_df = pd.read_csv(cons_file, sep="\t")
pos_col = cons_df.columns[0]
score_col = cons_df.columns[1]
cons_df = cons_df.rename(columns={pos_col: "pos", score_col: "score"})
cons_df["pos"] = cons_df["pos"].astype(int)

# Load variants
vars_df = pd.read_csv(vars_file, sep="\t")
vars_df["pos"] = vars_df["pos"].astype(int)

# Merge
merged = pd.merge(vars_df, cons_df, on="pos", how="left").dropna(subset=["score"])

groups = ["NTD", "benign"]
data = [merged.loc[merged["group"] == g, "score"] for g in groups]

plt.style.use("default")
fig, ax = plt.subplots(figsize=(5,4))

# Boxplot: clean, white, no fliers
box = ax.boxplot(
    data,
    positions=[1, 2],
    widths=0.5,
    labels=["NTD variants", "Benign variants"],
    patch_artist=True,
    showfliers=False
)

for b in box["boxes"]:
    b.set(facecolor="white", edgecolor="black", linewidth=1.2)
for whisker in box["whiskers"]:
    whisker.set(color="black", linewidth=1.0)
for cap in box["caps"]:
    cap.set(color="black", linewidth=1.0)
for median in box["medians"]:
    median.set(color="black", linewidth=1.4)

# Jittered points: small, light gray
for i, g in enumerate(groups, start=1):
    y = merged.loc[merged["group"] == g, "score"]
    x = np.random.normal(i, 0.04, size=len(y))
    ax.scatter(x, y, color="gray", alpha=0.6, s=18, zorder=3)

ax.set_ylabel("Constraint score")
ax.set_title("Constraint at VANGL1 variant sites")
ax.set_xlim(0.5, 2.5)
ax.grid(axis="y", linestyle="--", alpha=0.4)
plt.tight_layout()

outfile = OUTPUT_DIR / "Figure2_variant_constraint_boxplot.png"
plt.savefig(outfile, dpi=300)
print("Saved:", outfile)
