"""
=============================================================================
AFP / Cryoprotective Peptide Screening Pipeline — FreeSASA + BioPython
=============================================================================

DESCRIPTION:
    Automated analysis of multiple PDB files to identify antifreeze protein
    (AFP) candidates based on:
      1. Per-residue Solvent Accessible Surface Area (SASA)
      2. Hydrophobic patch detection and scoring
      3. IBS-typical residue exposure profiling
      4. Ice lattice complementarity (4.52 Å spacing matches)
      5. B-factor rigidity analysis
      6. Summary ranking table + CSV export + plots

INSTALLATION (run once):
    pip install freesasa biopython pandas scipy matplotlib seaborn

=============================================================================
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import cdist

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

PDB_DIR     = "./PDB"          # folder containing your .pdb files
OUTPUT_DIR  = "./afp_results"   # where results will be saved

# IBS-typical residues (known to ice-binding)
IBS_RESIDUES = {"THR", "ALA", "VAL", "SER", "GLY", "LEU", "ILE"}

# Ice lattice spacing targets (Å) — Ih ice a-axis and c-axis
ICE_LATTICE_A = 4.52   # basal / prism face
ICE_LATTICE_C = 7.36   # pyramidal / secondary prism face
LATTICE_TOLERANCE = 0.5  # Å — acceptable deviation

# SASA exposure thresholds for IBS candidates (%)
MIN_EXPOSURE = 10.0   
MAX_EXPOSURE = 120.0   

# B-factor threshold for rigidity
BFACTOR_RIGID = 30.0  # Å² — below this = rigid = good IBS (look for source)

# ─────────────────────────────────────────────────────────────────────────────
# SETUP
# ─────────────────────────────────────────────────────────────────────────────

os.makedirs(OUTPUT_DIR, exist_ok=True)

try:
    import freesasa
    HAS_FREESASA = True
except ImportError:
    HAS_FREESASA = False
    print("WARNING: freesasa not installed. Install with: pip install freesasa")
    print("         SASA analysis will be skipped.\n")

# adding error messages

try:
    from Bio.PDB import PDBParser, DSSP
    from Bio.PDB.DSSP import dssp_dict_from_pdb_file
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("WARNING: biopython not installed. Install with: pip install biopython")
    print("         B-factor and structure analysis will be limited.\n")


# ─────────────────────────────────────────────────────────────────────────────
# CORE ANALYSIS FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def run_freesasa(pdb_path):
    """
    Run FreeSASA on a PDB file.
    Returns a dict: {residue_key: {'abs_sasa': float, 'rel_sasa': float,
                                   'resname': str, 'chain': str, 'resid': int}}
    """
    if not HAS_FREESASA:
        return {}

    try:
        structure = freesasa.Structure(pdb_path)
        result    = freesasa.calc(structure)
        residue_areas = result.residueAreas()

        residues = {}
        for chain_id, chain_data in residue_areas.items():
            for res_id, area in chain_data.items():
                key = f"{area.residueType}{chain_id}{res_id}"
                residues[key] = {
                    "resname"  : area.residueType.strip(),
                    "chain"    : chain_id,
                    "resid"    : int(res_id),
                    "abs_sasa" : area.total,
                    "rel_sasa" : area.relativeTotal * 100  # convert to %
                }
        return residues

    except Exception as e:
        print(f"  [SASA ERROR] {os.path.basename(pdb_path)}: {e}")
        return {}


def get_bfactors(pdb_path):
    """
    Extract mean B-factor per residue from PDB ATOM records.
    Returns dict: {(chain, resid): mean_bfactor}
    """
    bfactors = {}
    current  = None
    vals     = []

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            try:
                chain = line[21].strip()
                resid = int(line[22:26].strip())
                bval  = float(line[60:66].strip())
                key   = (chain, resid)

                if key != current:
                    if current is not None and vals:
                        bfactors[current] = np.mean(vals)
                    current = key
                    vals = [bval]
                else:
                    vals.append(bval)
            except (ValueError, IndexError):
                continue

    if current is not None and vals:
        bfactors[current] = np.mean(vals)

    return bfactors


def get_cb_coords(pdb_path, residue_names=None):
    """
    Extract Cβ (or Cα for Gly) coordinates for specified residue types.
    Returns list of dicts with residue info and coordinates.
    """
    if residue_names is None:
        residue_names = IBS_RESIDUES

    coords = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            try:
                resname = line[17:20].strip()
                chain   = line[21].strip()
                resid   = int(line[22:26].strip())
                atom    = line[12:16].strip()
                x       = float(line[30:38].strip())
                y       = float(line[38:46].strip())
                z       = float(line[46:54].strip())

                # Use CB for most residues, CA for Gly
                target_atom = "CA" if resname == "GLY" else "CB"
                if atom == target_atom and resname in residue_names:
                    coords.append({
                        "resname": resname,
                        "chain"  : chain,
                        "resid"  : resid,
                        "coords" : np.array([x, y, z])
                    })
            except (ValueError, IndexError):
                continue

    return coords


def score_lattice_complementarity(cb_list):
    """
    Count pairs of IBS-type residues with Cβ spacing matching ice lattice.
    Returns dict with counts and a normalized score.
    """
    if len(cb_list) < 2:
        return {"a_axis_matches": 0, "c_axis_matches": 0,
                "total_matches": 0, "lattice_score": 0.0}

    coord_array = np.array([r["coords"] for r in cb_list])
    dmat = cdist(coord_array, coord_array)

    # Count matches to a-axis (4.52 Å) and c-axis (7.36 Å)
    a_matches = np.sum(
        (dmat > ICE_LATTICE_A - LATTICE_TOLERANCE) &
        (dmat < ICE_LATTICE_A + LATTICE_TOLERANCE) &
        (dmat > 0)
    ) // 2

    c_matches = np.sum(
        (dmat > ICE_LATTICE_C - LATTICE_TOLERANCE) &
        (dmat < ICE_LATTICE_C + LATTICE_TOLERANCE) &
        (dmat > 0)
    ) // 2

    total  = int(a_matches + c_matches)
    n_res  = len(cb_list)
    # Normalize by number of possible pairs
    n_pairs = max((n_res * (n_res - 1)) / 2, 1)
    score  = round(total / n_pairs, 4)

    return {
        "a_axis_matches" : int(a_matches),
        "c_axis_matches" : int(c_matches),
        "total_matches"  : total,
        "lattice_score"  : score
    }


def score_sasa(residue_sasa):
    """
    Score a peptide based on its SASA profile.
    Returns dict of AFP-relevant surface metrics.
    """
    if not residue_sasa:
        return {
            "n_ibs_candidates"   : 0,
            "mean_hydrophobic_sasa": 0.0,
            "hydrophobic_ratio"  : 0.0,
            "ibs_residue_list"   : ""
        }

    total_sasa     = sum(r["abs_sasa"] for r in residue_sasa.values())
    hydrophobic_res = {"ALA", "VAL", "ILE", "LEU", "PHE", "MET", "TRP", "PRO", "GLY"}

    hydrophobic_sasa = sum(
        r["abs_sasa"] for r in residue_sasa.values()
        if r["resname"] in hydrophobic_res
    )

    # IBS candidates: IBS-type residues with moderate exposure
    ibs_candidates = [
        r for r in residue_sasa.values()
        if r["resname"] in IBS_RESIDUES
        and MIN_EXPOSURE <= r["rel_sasa"] <= MAX_EXPOSURE
    ]

    ibs_list = ", ".join(
        f"{r['resname']}{r['resid']}" for r in
        sorted(ibs_candidates, key=lambda x: x["resid"])
    )

    hydro_ratio = round(hydrophobic_sasa / total_sasa, 3) if total_sasa > 0 else 0.0

    return {
        "total_sasa"           : round(total_sasa, 2),
        "hydrophobic_sasa"     : round(hydrophobic_sasa, 2),
        "hydrophobic_ratio"    : hydro_ratio,
        "n_ibs_candidates"     : len(ibs_candidates),
        "ibs_residue_list"     : ibs_list
    }


def score_bfactors(bfactor_dict, cb_list):
    """
    Score rigidity of IBS-candidate residues using B-factors.
    """
    if not bfactor_dict or not cb_list:
        return {"mean_ibs_bfactor": None, "rigid_ibs_fraction": None}

    ibs_bvals = []
    for res in cb_list:
        key = (res["chain"], res["resid"])
        if key in bfactor_dict:
            ibs_bvals.append(bfactor_dict[key])

    if not ibs_bvals:
        return {"mean_ibs_bfactor": None, "rigid_ibs_fraction": None}

    mean_b   = round(np.mean(ibs_bvals), 2)
    rigid_fr = round(sum(1 for b in ibs_bvals if b < BFACTOR_RIGID) / len(ibs_bvals), 3)

    return {
        "mean_ibs_bfactor"    : mean_b,
        "rigid_ibs_fraction"  : rigid_fr   # 1.0 = all rigid, 0.0 = all flexible
    }


# ─────────────────────────────────────────────────────────────────────────────
# SCORING SYSTEM
# ─────────────────────────────────────────────────────────────────────────────

def compute_afp_score(sasa_scores, lattice_scores, bfactor_scores):
    """
    Combine all metrics into a single AFP likelihood score (0–10).

    Score components:
      - Hydrophobic ratio        (0–2 pts)
      - N IBS candidates         (0–2 pts)
      - Lattice complementarity  (0–3 pts)
      - Rigidity (B-factor)      (0–3 pts)
    """
    score = 0.0
    breakdown = {}

    # 1. Hydrophobic ratio (0–2)
    hr = sasa_scores.get("hydrophobic_ratio", 0)
    if hr >= 0.45:
        pts = 2.0
    elif hr >= 0.30:
        pts = 1.0
    else:
        pts = 0.0
    score += pts
    breakdown["hydrophobic_ratio_pts"] = pts

    # 2. IBS candidate count (0–2)
    n = sasa_scores.get("n_ibs_candidates", 0)
    if n >= 5:
        pts = 2.0
    elif n >= 3:
        pts = 1.0
    else:
        pts = 0.0
    score += pts
    breakdown["ibs_candidates_pts"] = pts

    # 3. Lattice complementarity (0–3)
    ls = lattice_scores.get("lattice_score", 0)
    if ls >= 0.15:
        pts = 3.0
    elif ls >= 0.07:
        pts = 2.0
    elif ls >= 0.02:
        pts = 1.0
    else:
        pts = 0.0
    score += pts
    breakdown["lattice_score_pts"] = pts

    # 4. B-factor rigidity (0–3)
    rf = bfactor_scores.get("rigid_ibs_fraction", None)
    if rf is None:
        pts = 1.5  # neutral if no B-factors
    elif rf >= 0.8:
        pts = 3.0
    elif rf >= 0.5:
        pts = 2.0
    elif rf >= 0.3:
        pts = 1.0
    else:
        pts = 0.0
    score += pts
    breakdown["bfactor_rigidity_pts"] = pts

    breakdown["total_score"] = round(score, 2)
    return breakdown


# ─────────────────────────────────────────────────────────────────────────────
# MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def analyze_pdb(pdb_path):
    """Run full AFP analysis on a single PDB file."""
    name = os.path.splitext(os.path.basename(pdb_path))[0]
    print(f"  Analyzing: {name}")

    # 1. SASA
    residue_sasa = run_freesasa(pdb_path)
    sasa_scores  = score_sasa(residue_sasa)

    # 2. Ice lattice complementarity
    cb_list         = get_cb_coords(pdb_path)
    lattice_scores  = score_lattice_complementarity(cb_list)

    # 3. B-factor rigidity
    bfactor_dict   = get_bfactors(pdb_path)
    bfactor_scores = score_bfactors(bfactor_dict, cb_list)

    # 4. Combined AFP score
    afp_score = compute_afp_score(sasa_scores, lattice_scores, bfactor_scores)

    # 5. Save per-residue data
    if residue_sasa:
        res_df = pd.DataFrame(residue_sasa.values())
        res_df["is_IBS_type"]   = res_df["resname"].isin(IBS_RESIDUES)
        res_df["is_candidate"]  = (
            res_df["is_IBS_type"] &
            res_df["rel_sasa"].between(MIN_EXPOSURE, MAX_EXPOSURE)
        )
        res_df["bfactor"] = res_df.apply(
            lambda row: bfactor_dict.get((row["chain"], row["resid"]), None), axis=1
        )
        res_df = res_df.sort_values("resid")
        res_df.to_csv(f"{OUTPUT_DIR}/{name}_residues.csv", index=False)

    # Merge all results
    result = {"peptide": name}
    result.update(sasa_scores)
    result.update(lattice_scores)
    result.update(bfactor_scores)
    result.update(afp_score)

    return result


def run_pipeline():
    pdb_files = sorted(glob.glob(os.path.join(PDB_DIR, "*.pdb")))

    if not pdb_files:
        print(f"\nNo PDB files found in '{PDB_DIR}'")
        print("Please set PDB_DIR to the folder containing your .pdb files.")
        return

    print(f"\n{'='*60}")
    print(f"  AFP Screening Pipeline")
    print(f"  Found {len(pdb_files)} PDB file(s) in '{PDB_DIR}'")
    print(f"{'='*60}\n")

    results = []
    for pdb_path in pdb_files:
        result = analyze_pdb(pdb_path)
        results.append(result)

    df = pd.DataFrame(results)
    df = df.sort_values("total_score", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1

    # ── Save CSV ──────────────────────────────────────────────────────────────
    csv_path = f"{OUTPUT_DIR}/summary_ranking.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nCSV saved: {csv_path}")

    # ── Save human-readable report ────────────────────────────────────────────
    report_path = f"{OUTPUT_DIR}/summary_ranking.txt"
    with open(report_path, "w") as f:
        f.write("AFP CANDIDATE RANKING REPORT\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"{'Rank':<5} {'Peptide':<20} {'Score/10':<10} "
                f"{'Hydro%':<10} {'IBS_res':<10} "
                f"{'LatticeScore':<14} {'RigidFrac':<12}\n")
        f.write("-" * 80 + "\n")

        for _, row in df.iterrows():
            rigid = f"{row['rigid_ibs_fraction']:.2f}" \
                    if pd.notna(row.get('rigid_ibs_fraction')) else "N/A"
            f.write(
                f"{int(row['rank']):<5} "
                f"{row['peptide']:<20} "
                f"{row['total_score']:<10.2f} "
                f"{row.get('hydrophobic_ratio', 0):<10.3f} "
                f"{int(row.get('n_ibs_candidates', 0)):<10} "
                f"{row.get('lattice_score', 0):<14.4f} "
                f"{rigid:<12}\n"
            )

        f.write("\n\nDETAILED RESULTS\n")
        f.write("=" * 60 + "\n\n")
        for _, row in df.iterrows():
            f.write(f"[#{int(row['rank'])}] {row['peptide']}  —  "
                    f"Score: {row['total_score']}/10\n")
            f.write(f"  Hydrophobic ratio    : {row.get('hydrophobic_ratio', 'N/A')}\n")
            f.write(f"  IBS candidates (n)   : {row.get('n_ibs_candidates', 'N/A')}\n")
            f.write(f"  IBS residues         : {row.get('ibs_residue_list', 'N/A')}\n")
            f.write(f"  Lattice score        : {row.get('lattice_score', 'N/A')}\n")
            f.write(f"  a-axis matches (4.52Å): {row.get('a_axis_matches', 'N/A')}\n")
            f.write(f"  c-axis matches (7.36Å): {row.get('c_axis_matches', 'N/A')}\n")
            f.write(f"  Mean IBS B-factor    : {row.get('mean_ibs_bfactor', 'N/A')}\n")
            f.write(f"  Rigid IBS fraction   : {row.get('rigid_ibs_fraction', 'N/A')}\n")
            f.write(f"  Score breakdown:\n")
            f.write(f"    Hydrophobic ratio pts : {row.get('hydrophobic_ratio_pts', 'N/A')}/2\n")
            f.write(f"    IBS candidates pts    : {row.get('ibs_candidates_pts', 'N/A')}/2\n")
            f.write(f"    Lattice score pts     : {row.get('lattice_score_pts', 'N/A')}/3\n")
            f.write(f"    B-factor rigidity pts : {row.get('bfactor_rigidity_pts', 'N/A')}/3\n")
            f.write("\n")

    print(f"Report saved: {report_path}")

    # ── Plots ─────────────────────────────────────────────────────────────────
    plot_results(df)

    # ── Console summary ───────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("  RANKING SUMMARY")
    print(f"{'='*60}")
    print(f"{'Rank':<5} {'Peptide':<20} {'AFP Score':>10}")
    print("-" * 38)
    for _, row in df.iterrows():
        stars = "★" * int(round(row["total_score"] / 2))
        print(f"  {int(row['rank'])}    {row['peptide']:<20} "
              f"{row['total_score']:>5.1f}/10  {stars}")
    print(f"\nTop candidate: {df.iloc[0]['peptide']} "
          f"(score: {df.iloc[0]['total_score']}/10)")
    print(f"\nAll results saved to: {OUTPUT_DIR}/\n")


# ─────────────────────────────────────────────────────────────────────────────
# PLOTTING
# ─────────────────────────────────────────────────────────────────────────────

def plot_results(df):
    """Generate summary visualizations."""

    # ── 1. Score breakdown heatmap ────────────────────────────────────────────
    score_cols = [
        "hydrophobic_ratio_pts",
        "ibs_candidates_pts",
        "lattice_score_pts",
        "bfactor_rigidity_pts"
    ]
    score_cols = [c for c in score_cols if c in df.columns]

    if score_cols:
        fig, axes = plt.subplots(1, 2, figsize=(14, max(4, len(df) * 0.6 + 2)))

        # Heatmap
        hmap_data = df.set_index("peptide")[score_cols].astype(float)
        hmap_data.columns = ["Hydrophobic\nRatio", "IBS\nCandidates",
                             "Lattice\nScore", "B-factor\nRigidity"]
        sns.heatmap(
            hmap_data, ax=axes[0],
            cmap="YlOrRd", annot=True, fmt=".1f",
            linewidths=0.5, cbar_kws={"label": "Points"}
        )
        axes[0].set_title("Score Breakdown per Criterion", fontweight="bold")
        axes[0].set_ylabel("")

        # Total score bar chart
        colors = plt.cm.RdYlGn(df["total_score"] / 10)
        axes[1].barh(df["peptide"][::-1], df["total_score"][::-1], color=colors[::-1])
        axes[1].set_xlabel("Total AFP Score (out of 10)")
        axes[1].set_title("Overall AFP Likelihood Score", fontweight="bold")
        axes[1].axvline(x=5, color="gray", linestyle="--", alpha=0.5, label="Threshold")
        axes[1].legend()
        for i, (score, name) in enumerate(zip(df["total_score"][::-1],
                                               df["peptide"][::-1])):
            axes[1].text(score + 0.1, i, f"{score:.1f}", va="center", fontsize=9)

        plt.tight_layout()
        heatmap_path = f"{OUTPUT_DIR}/heatmap_scores.png"
        plt.savefig(heatmap_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Heatmap saved: {heatmap_path}")

    # ── 2. Per-residue SASA profiles ─────────────────────────────────────────
    res_files = glob.glob(f"{OUTPUT_DIR}/*_residues.csv")
    if res_files:
        n = len(res_files)
        fig, axes = plt.subplots(n, 1, figsize=(12, 3 * n), squeeze=False)

        for i, res_file in enumerate(sorted(res_files)):
            ax = axes[i][0]
            rdf = pd.read_csv(res_file)
            rdf = rdf.sort_values("resid")

            bar_colors = [
                "#e74c3c" if row["is_candidate"]
                else "#3498db" if row["is_IBS_type"]
                else "#bdc3c7"
                for _, row in rdf.iterrows()
            ]

            ax.bar(
                [f"{r['resname']}{r['resid']}" for _, r in rdf.iterrows()],
                rdf["rel_sasa"],
                color=bar_colors,
                edgecolor="white", linewidth=0.3
            )
            ax.axhline(y=MIN_EXPOSURE, color="orange", linestyle="--",
                       alpha=0.7, label=f"Min exposure ({MIN_EXPOSURE}%)")
            ax.axhline(y=MAX_EXPOSURE, color="red", linestyle="--",
                       alpha=0.7, label=f"Max exposure ({MAX_EXPOSURE}%)")
            ax.set_ylabel("Relative SASA (%)")
            peptide_name = os.path.basename(res_file).replace("_residues.csv", "")
            ax.set_title(f"{peptide_name} — Residue Exposure Profile",
                        fontweight="bold")
            ax.tick_params(axis="x", rotation=45, labelsize=7)
            ax.legend(fontsize=8)

            # Legend patches
            from matplotlib.patches import Patch
            legend_els = [
                Patch(facecolor="#e74c3c", label="IBS candidate"),
                Patch(facecolor="#3498db", label="IBS-type residue"),
                Patch(facecolor="#bdc3c7", label="Other residue"),
            ]
            ax.legend(handles=legend_els, fontsize=8, loc="upper right")

        plt.tight_layout()
        profile_path = f"{OUTPUT_DIR}/sasa_profiles.png"
        plt.savefig(profile_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"SASA profiles saved: {profile_path}")


# ─────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    run_pipeline()
