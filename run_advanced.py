import json
from pathlib import Path
import subprocess
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem

BASE = Path(r"C:\dhfr-natlead")
vina_exe = BASE / "tools" / "vina" / "vina.exe"
raw_pdb_dir = BASE / "data" / "raw" / "pdb"
lig_dir = BASE / "data" / "processed" / "ligands"
out_dir = BASE / "results" / "docking"
adv_dir = BASE / "results" / "advanced"
adv_dir.mkdir(parents=True, exist_ok=True)


def convert_mol_to_pdbqt(mol, out_path):
  mol_h = Chem.AddHs(mol, addCoords=True)
  preparator = MoleculePreparation()
  mol_setups = preparator.prepare(mol_h)
  setup = mol_setups[0] if isinstance(mol_setups, list) else preparator.setup
  result = PDBQTWriterLegacy.write_string(setup)
  pdbqt_text = result[0] if isinstance(result, tuple) else result
  with open(out_path, "w") as f:
    f.write(pdbqt_text)


def prep_receptor(raw_pdb, out_pdbqt):
  with open(raw_pdb) as f:
    lines = f.readlines()
  pdbqt_lines = [
      l[:54]
      + "  1.00  0.00          "
      + (
          l[76:78].strip().rjust(2)
          if len(l) >= 78
          else l[12:14].strip().rjust(2)
      )
      + "\n"
      for l in lines
      if l.startswith(("ATOM", "HETATM"))
  ]
  with open(out_pdbqt, "w") as f:
    f.writelines(pdbqt_lines)


# ==========================================
# 1. HUMAN DHFR SELECTIVITY (PDB: 1U72)
# ==========================================
print("=== 1. PREPARING HUMAN DHFR (PDB: 1U72) ===")
hu_raw = raw_pdb_dir / "1U72.pdb"
hu_pdbqt = adv_dir / "1U72_human_dhfr.pdbqt"
prep_receptor(hu_raw, hu_pdbqt)

# Derive grid center from bound Methotrexate (MTX) in 1U72
hu_lines = open(hu_raw).readlines()
mtx_coords = [
    [float(l[30:38]), float(l[38:46]), float(l[46:54])]
    for l in hu_lines
    if l.startswith("HETATM") and "MTX" in l and len(l) >= 54
]
if mtx_coords:
  hu_center = [
      round(sum(c[i] for c in mtx_coords) / len(mtx_coords), 2)
      for i in range(3)
  ]
else:
  hu_center = [27.0, 30.0, 15.0]

print(f"Human DHFR Grid Center: {hu_center}")

# Dock key compounds against Human DHFR
test_compounds = ["Baicalein", "Berberine", "Curcumin", "Trimethoprim"]
human_docking = []

for comp in test_compounds:
  lig_pdbqt = lig_dir / f"{comp}.pdbqt"
  if not lig_pdbqt.exists():
    continue
  comp_out = adv_dir / f"{comp}_huDHFR.pdbqt"
  cmd = [
      str(vina_exe),
      "--receptor",
      str(hu_pdbqt),
      "--ligand",
      str(lig_pdbqt),
      "--center_x",
      str(hu_center[0]),
      "--center_y",
      str(hu_center[1]),
      "--center_z",
      str(hu_center[2]),
      "--size_x",
      "20.0",
      "--size_y",
      "20.0",
      "--size_z",
      "20.0",
      "--out",
      str(comp_out),
      "--exhaustiveness",
      "8",
      "--num_modes",
      "1",
  ]
  res = subprocess.run(cmd, capture_output=True, text=True)
  score = None
  for line in res.stdout.splitlines():
    parts = line.strip().split()
    if len(parts) >= 2 and parts[0] == "1":
      try:
        score = float(parts[1])
        break
      except ValueError:
        pass
  human_docking.append({"compound": comp, "human_affinity_kcal_mol": score})
  print(f"   -> Human DHFR vs {comp:<15}: {score} kcal/mol")

with open(adv_dir / "human_selectivity.json", "w") as f:
  json.dump(human_docking, f, indent=2)

# ==========================================
# 2. F98Y RESISTANCE (PDB: 3F0B WT vs 3F0U MUT)
# ==========================================
print("\n=== 2. RUNNING F98Y RESISTANCE ANALYSIS (3F0B vs 3F0U) ===")
wt_raw = raw_pdb_dir / "3F0B.pdb"
mut_raw = raw_pdb_dir / "3F0U.pdb"
wt_pdbqt = adv_dir / "3F0B_WT.pdbqt"
mut_pdbqt = adv_dir / "3F0U_F98Y.pdbqt"
prep_receptor(wt_raw, wt_pdbqt)
prep_receptor(mut_raw, mut_pdbqt)

# Grid center around NADPH/inhibitor pocket
wt_lines = open(wt_raw).readlines()
ndp_coords = [
    [float(l[30:38]), float(l[38:46]), float(l[46:54])]
    for l in wt_lines
    if l.startswith("HETATM") and "NDP" in l and len(l) >= 54
]
f98_center = (
    [round(sum(c[i] for c in ndp_coords) / len(ndp_coords), 2) for i in range(3)]
    if ndp_coords
    else [15.0, -10.0, 25.0]
)

resistance_results = []
for comp in ["Trimethoprim", "Baicalein", "Berberine"]:
  lig_pdbqt = lig_dir / f"{comp}.pdbqt"
  if not lig_pdbqt.exists():
    continue

  # WT Docking
  cmd_wt = [
      str(vina_exe),
      "--receptor",
      str(wt_pdbqt),
      "--ligand",
      str(lig_pdbqt),
      "--center_x",
      str(f98_center[0]),
      "--center_y",
      str(f98_center[1]),
      "--center_z",
      str(f98_center[2]),
      "--size_x",
      "22.0",
      "--size_y",
      "22.0",
      "--size_z",
      "22.0",
      "--out",
      str(adv_dir / f"{comp}_WT.pdbqt"),
      "--exhaustiveness",
      "8",
      "--num_modes",
      "1",
  ]
  res_wt = subprocess.run(cmd_wt, capture_output=True, text=True)
  s_wt = next(
      (
          float(l.split()[1])
          for l in res_wt.stdout.splitlines()
          if l.strip().startswith("1 ")
      ),
      None,
  )

  # F98Y Mutant Docking
  cmd_mut = [
      str(vina_exe),
      "--receptor",
      str(mut_pdbqt),
      "--ligand",
      str(lig_pdbqt),
      "--center_x",
      str(f98_center[0]),
      "--center_y",
      str(f98_center[1]),
      "--center_z",
      str(f98_center[2]),
      "--size_x",
      "22.0",
      "--size_y",
      "22.0",
      "--size_z",
      "22.0",
      "--out",
      str(adv_dir / f"{comp}_F98Y.pdbqt"),
      "--exhaustiveness",
      "8",
      "--num_modes",
      "1",
  ]
  res_mut = subprocess.run(cmd_mut, capture_output=True, text=True)
  s_mut = next(
      (
          float(l.split()[1])
          for l in res_mut.stdout.splitlines()
          if l.strip().startswith("1 ")
      ),
      None,
  )

  delta = round(s_mut - s_wt, 3) if (s_wt and s_mut) else None
  resistance_results.append({
      "compound": comp,
      "wildtype_kcal_mol": s_wt,
      "f98y_mutant_kcal_mol": s_mut,
      "affinity_loss_kcal_mol": delta,
  })
  print(
      f"   -> {comp:<15}: WT={s_wt} kcal/mol | F98Y={s_mut} kcal/mol |"
      f" Diff={delta} kcal/mol"
  )

with open(adv_dir / "f98y_resistance.json", "w") as f:
  json.dump(resistance_results, f, indent=2)

print("\n=== ADVANCED EXTENSION CALCULATIONS COMPLETE ===")