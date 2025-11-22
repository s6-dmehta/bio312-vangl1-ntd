import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

DATA_DIR = Path("data")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

cons_file = DATA_DIR / "VANGL1_conservation.tsv"
domains_file = DATA_DIR / "VANGL1_domains.tsv"

# -------------------------
# LOAD CONSERVATION DATA
# -------------------------
cons_df = pd.read_csv(cons_file, sep="\t")
pos_col = cons_df.columns[0]
score_col = cons_df.columns[1]

cons_df = cons_df.sort_values(pos_col)
cons_df["smooth"] = (
    cons_df[score_col]
    .rolling(window=7, center=True, min_periods=1)
    .mean()
)

# -------------------------
# LOAD DOMAIN DATA
# -------------------------
domains_df = pd.read_csv(domains_file, sep="\t")
start_col = "start"
end_col = "end"
name_col = "name"

# -------------------------
# PLOT
# -------------------------
plt.style.use("default")
fig, ax = plt.subplots(figsize=(10,4))

# Shade TM helices
tm_rows = domains_df[domains_df[name_col].str.contains("TM", case=False, na=False)]
for _, r in tm_rows.iterrows():
    ax.axvspan(r[start_col], r[end_col], color="lightgray", alpha=0.4)

# Shade PDZ-binding motif
pdz_rows = domains_df[domains_df[name_col].str.contains("PDZ", case=False, na=False)]
for _, r in pdz_rows.iterrows():
    ax.axvspan(r[start_col], r[end_col], color="orange", alpha=0.5)

# Plot conservation lines
ax.plot(cons_df[pos_col], cons_df[score_col], color="steelblue", lw=1, alpha=0.6)
ax.plot(cons_df[pos_col], cons_df["smooth"], color="navy", lw=2)

# Labels
ax.set_xlabel("Residue position")
ax.set_ylabel("Conservation score")
ax.set_title("VANGL1 conservation landscape")

ax.grid(axis="y", linestyle="--", alpha=0.5)
plt.tight_layout()

outfile = OUTPUT_DIR / "Figure1_conservation.png"
plt.savefig(outfile, dpi=300)
print("Saved:", outfile)
