import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

DATA = Path("data")
OUT = Path("figures")
OUT.mkdir(exist_ok=True)

# Load conservation
cons = pd.read_csv(DATA / "VANGL1_conservation.tsv", sep="\t")
cons.columns = ["pos", "score"]
cons["pos"] = cons["pos"].astype(int)

# Load domains
domains = pd.read_csv(DATA / "VANGL1_domains.tsv", sep="\t")

# Load variants
vars_df = pd.read_csv(DATA / "VANGL1_variants.tsv", sep="\t")
vars_df["pos"] = vars_df["pos"].astype(int)
ntd = vars_df[vars_df["group"] == "NTD"]
benign = vars_df[vars_df["group"] == "benign"]

plt.style.use("default")
fig, ax = plt.subplots(figsize=(10, 4))

# Plot conservation curve
ax.plot(cons["pos"], cons["score"], color="lightgray", lw=1.2, zorder=1)

# Shade TM helices and PDZ motif
for _, row in domains.iterrows():
    name = str(row["name"])
    if name.startswith("TM"):
        ax.axvspan(row["start"], row["end"], color="lightgray", alpha=0.3, zorder=0)
    if "PDZ" in name:
        ax.axvspan(row["start"], row["end"], color="orange", alpha=0.3, zorder=0)

# Plot benign variants (above the curve)
ax.scatter(
    benign["pos"],
    [1.05] * len(benign),
    color="steelblue",
    s=40,
    label="Benign variants",
    zorder=3,
)

# Plot NTD variants even higher, as stars
ax.scatter(
    ntd["pos"],
    [1.15] * len(ntd),
    color="crimson",
    s=80,
    marker="*",
    label="NTD-associated variants",
    zorder=4,
)

ax.set_xlim(cons["pos"].min() - 5, cons["pos"].max() + 5)
ax.set_ylim(0, 1.25)

ax.set_xlabel("Residue position")
ax.set_ylabel("Constraint / variant positions")
ax.set_title("Constraint landscape with NTD and benign variants mapped on VANGL1")

ax.legend(loc="upper left", frameon=False)
ax.grid(axis="y", linestyle="--", alpha=0.3)

plt.tight_layout()
outfile = OUT / "Figure5_variant_positions_on_constraint.png"
plt.savefig(outfile, dpi=300)
print("Saved Figure 5 to:", outfile)
