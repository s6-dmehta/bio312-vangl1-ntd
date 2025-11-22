import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from pathlib import Path
import numpy as np

STRUCT = Path("structure")
OUT = Path("figures")
OUT.mkdir(exist_ok=True)

pdb_file = STRUCT / "VANGL1_AF2.pdb"

# ---- Parse C-alpha coordinates from PDB ----
records = []
with open(pdb_file) as f:
    for line in f:
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        try:
            resnum = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        records.append((resnum, x, y, z))

if not records:
    raise RuntimeError("No C-alpha atoms found in PDB")

resnums = np.array([r[0] for r in records])
xs = np.array([r[1] for r in records])
ys = np.array([r[2] for r in records])
zs = np.array([r[3] for r in records])

# Color from N-terminus to C-terminus
colors = (resnums - resnums.min()) / (resnums.max() - resnums.min())

# ---- Make 3D plot ----
plt.style.use("default")
fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111, projection="3d")

# backbone line + points
ax.plot(xs, ys, zs, linewidth=1.5, alpha=0.7)
sc = ax.scatter(xs, ys, zs, c=colors, cmap="viridis", s=10)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("VANGL1 AlphaFold structure (Cα trace, colored N→C)")

cb = fig.colorbar(sc, ax=ax, shrink=0.6)
cb.set_label("Residue index (N-terminus → C-terminus)")

plt.tight_layout()
out_file = OUT / "Figure5_VANGL1_structure.png"
plt.savefig(out_file, dpi=300)
print("Saved Figure 5 to:", out_file)
