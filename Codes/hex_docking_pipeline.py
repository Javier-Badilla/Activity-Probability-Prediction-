"""
=============================================================================
HEX Docking Results Analysis Pipeline — Multi-Peptide, Per-Orientation  v1.2
=============================================================================
"""

import os
import glob
import shutil
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

warnings.filterwarnings("ignore")

# -----------------------------------------------------------------------------
# CONFIGURATION — edit these
# -----------------------------------------------------------------------------

# Parent folder containing one sub-folder per peptide, each holding
# the individual per-orientation PDB files from HEX.
# e.g.  results/pep42/pep42-1.pdb  results/pep42/pep42-2.pdb  ...
HEX_DIR    = "./results"
OUTPUT_DIR = "./hex_results"

# IBS residue numbers from your AFP pipeline output (ibs_residue_list column).
# If a peptide is not listed, ALL protein residues are used as fallback.
IBS_RESIDUE_MAP = {
    "pep153": [2, 6],
    "pep71" : [1, 3],
    "pep42" : [2, 3, 5, 8, 11, 12, 13, 14, 15],
    "pep40" : [2, 3, 4, 5, 6, 7, 11, 12, 14, 15],
    "pep11" : [3, 8, 9],
    "pep45" : [1, 2, 3, 5, 7, 8, 12, 13, 15],
    "pep24" : [1, 4, 7, 9],
    "pep32" : [5, 6, 8, 9]
}

# Ice/water residue names — "ICE" added for HEX slab output
KNOWN_ICE_RESNAMES = {
    "ICE", "TIP4", "TIP4P", "TP4", "WAT", "HOH",
    "SOL", "TIP3", "SPC", "SPCE", "TIP"
}

HBOND_MIN    = 2.4
HBOND_MAX    = 3.5
CONTACT_MAX  = 5.0
TOP_N_EXTRACT = 5
MIN_CONTACTS  = 3


# -----------------------------------------------------------------------------
# PDB PARSING HELPERS
# -----------------------------------------------------------------------------

def _parse_atom_line(line):
    """
    Parse one ATOM/HETATM line safely.
    Returns a dict or None if the line cannot be parsed.

    Handles blank chain ID (column 21) — common in HEX output for the ice slab.
    """
    try:
        record  = line[0:6].strip()
        name    = line[12:16].strip()
        resname = line[17:20].strip()
        chain   = line[21:22].strip()          # empty string when blank

        # Strip non-numeric chars so a missing chain does not corrupt resid
        raw_resid = line[22:26].strip()
        digits    = "".join(c for c in raw_resid if c.isdigit() or c == "-")
        resid     = int(digits) if digits else 0

        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())

        return dict(record=record, name=name, resname=resname,
                    chain=chain, resid=resid, x=x, y=y, z=z)
    except (ValueError, IndexError):
        return None


def _is_ice(atom):
    return atom["resname"] in KNOWN_ICE_RESNAMES


# -----------------------------------------------------------------------------
# UTILITY FUNCTIONS
# -----------------------------------------------------------------------------

def detect_ice_resnames(pdb_file):
    """Return sorted list of ice residue names found in a single PDB file."""
    found = set()
    with open(pdb_file) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom = _parse_atom_line(line)
            if atom and _is_ice(atom):
                found.add(atom["resname"])
    if not found:
        print(f"  WARNING: No ice residues detected in {os.path.basename(pdb_file)}")
        print(f"           Known names: {sorted(KNOWN_ICE_RESNAMES)}")
        print(f"           Check with: grep -E '^ATOM|^HETATM' <file> | cut -c18-20 | sort -u")
    return sorted(found)


def read_pdb_coords(pdb_file):
    """
    Read all ATOM/HETATM records from a single-orientation PDB file.
    Returns {"protein": [...], "ice": [...]}
    Each entry: (resid, resname, atom_name, x, y, z)
    """
    protein_atoms = []
    ice_atoms     = []

    with open(pdb_file) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom = _parse_atom_line(line)
            if atom is None:
                continue
            entry = (atom["resid"], atom["resname"], atom["name"],
                     atom["x"], atom["y"], atom["z"])
            if _is_ice(atom):
                ice_atoms.append(entry)
            else:
                protein_atoms.append(entry)

    return {"protein": protein_atoms, "ice": ice_atoms}


# -----------------------------------------------------------------------------
# CONTACT CALCULATION
# -----------------------------------------------------------------------------

def compute_contacts(coords_dict, ibs_residues, contact_cutoff=CONTACT_MAX):
    """
    Count IBS-residue atoms within H-bond and vdW distance of ice oxygens.
    Uses only oxygen atoms from the ice slab (atom name starts with O).
    """
    protein = coords_dict["protein"]
    if ibs_residues:
        prot_coords = np.array([
            [a[3], a[4], a[5]] for a in protein if a[0] in ibs_residues
        ])
    else:
        prot_coords = np.array([[a[3], a[4], a[5]] for a in protein])

    # Prefer oxygen atoms for accurate H-bond distances
    ice_coords = np.array([
        [a[3], a[4], a[5]] for a in coords_dict["ice"] if a[2].startswith("O")
    ])
    if len(ice_coords) == 0:  # fall back to all ice atoms
        ice_coords = np.array([[a[3], a[4], a[5]] for a in coords_dict["ice"]])

    if len(prot_coords) == 0 or len(ice_coords) == 0:
        return dict(n_contacts=0, hbond_contacts=0, vdw_contacts=0,
                    min_dist_A=999.0, passes=False)

    min_dist = 999.0
    n_hbond  = 0
    n_vdw    = 0
    chunk    = 500

    for i in range(0, len(prot_coords), chunk):
        blk  = prot_coords[i : i + chunk]
        diff = blk[:, np.newaxis, :] - ice_coords[np.newaxis, :, :]
        d    = np.sqrt((diff ** 2).sum(axis=2))
        chunk_min = d.min()
        if chunk_min < min_dist:
            min_dist = chunk_min
        n_hbond += int(np.sum((d >= HBOND_MIN) & (d <= HBOND_MAX)))
        n_vdw   += int(np.sum((d >  HBOND_MAX) & (d <= contact_cutoff)))

    n_total = n_hbond + n_vdw
    passes  = (n_total >= MIN_CONTACTS) and (HBOND_MIN < min_dist < contact_cutoff)

    return dict(n_contacts=n_total, hbond_contacts=n_hbond, vdw_contacts=n_vdw,
                min_dist_A=round(float(min_dist), 3), passes=passes)


# -----------------------------------------------------------------------------
# PER-PEPTIDE ANALYSIS
# -----------------------------------------------------------------------------

def analyze_peptide(peptide_dir, peptide_name, ibs_residues, peptide_output_dir):
    """
    Analyse all per-orientation PDB files in peptide_dir.
    Each file is one docking orientation produced by HEX.
    """
    os.makedirs(peptide_output_dir, exist_ok=True)

    pdb_files = sorted(glob.glob(os.path.join(peptide_dir, "*.pdb")))
    n_models  = len(pdb_files)
    if n_models == 0:
        print(f"  WARNING: No PDB files found in {peptide_dir}")
        return None

    # Detect ice residue names from the first file
    ice_resnames = detect_ice_resnames(pdb_files[0])

    print(f"\n  Peptide   : {peptide_name}")
    print(f"  Models    : {n_models}")
    print(f"  Ice found : {ice_resnames if ice_resnames else 'NONE DETECTED'}")
    print(f"  IBS resids: {ibs_residues if ibs_residues else 'all protein residues'}")

    results = []
    for idx, pdb_file in enumerate(pdb_files, 1):
        # Use the filename (without extension) as the model label
        model_label = os.path.splitext(os.path.basename(pdb_file))[0]
        coords  = read_pdb_coords(pdb_file)
        metrics = compute_contacts(coords, set(ibs_residues))
        metrics["model"]      = idx           # numeric rank for plotting
        metrics["model_file"] = model_label   # original filename stem
        metrics["pdb_path"]   = pdb_file
        results.append(metrics)
        if idx % 10 == 0:
            print(f"    Processed {idx}/{n_models} orientations...", end="\r")
    print(f"    Processed {n_models}/{n_models} orientations.   ")

    df = pd.DataFrame(results)
    df = df.sort_values(
        ["hbond_contacts", "n_contacts", "min_dist_A"],
        ascending=[False, False, True]
    ).reset_index(drop=True)
    df["rank"] = df.index + 1

    df.to_csv(os.path.join(peptide_output_dir, "orientation_ranking.csv"), index=False)

    # Copy top N PDB files directly (no extraction needed)
    top_models_dir = os.path.join(peptide_output_dir, "top_models")
    os.makedirs(top_models_dir, exist_ok=True)
    top_rows = df.head(TOP_N_EXTRACT)
    top_labels = top_rows["model_file"].tolist()
    for _, row in top_rows.iterrows():
        shutil.copy(row["pdb_path"], os.path.join(top_models_dir,
                                                   os.path.basename(row["pdb_path"])))

    plot_contact_profile(df, peptide_name, peptide_output_dir, top_rows["model"].tolist())

    best    = df.iloc[0]
    passing = int(df["passes"].sum())

    print(f"\n  RESULTS FOR {peptide_name}")
    print(f"  {'─'*52}")
    print(f"  Best orientation : {best['model_file']}")
    print(f"  H-bond contacts  : {int(best['hbond_contacts'])}")
    print(f"  Total contacts   : {int(best['n_contacts'])}")
    print(f"  Min distance     : {best['min_dist_A']:.3f} Å")
    print(f"  Orientations passing filter: {passing}/{n_models}")
    print(f"\n  {'':2} {'Rank':<6} {'File':<20} {'H-bonds':<10} {'Contacts':<10} {'Min dist':>10}")
    print(f"  {'─'*60}")
    for _, row in df.head(10).iterrows():
        flag = "✓" if row["passes"] else "✗"
        print(f"  {flag}  #{int(row['rank']):<5} {row['model_file']:<20} "
              f"{int(row['hbond_contacts']):<10} {int(row['n_contacts']):<10} "
              f"{row['min_dist_A']:>8.3f} Å")

    return dict(
        peptide              = peptide_name,
        n_models             = n_models,
        best_model           = best["model_file"],
        best_hbond_contacts  = int(best["hbond_contacts"]),
        best_total_contacts  = int(best["n_contacts"]),
        best_min_dist_A      = best["min_dist_A"],
        orientations_passing = passing,
        pass_rate_pct        = round(passing / n_models * 100, 1),
        top_models           = top_labels,
        ice_resnames_found   = ", ".join(ice_resnames) if ice_resnames else "NONE",
    )


# -----------------------------------------------------------------------------
# PLOTTING
# -----------------------------------------------------------------------------

def plot_contact_profile(df, peptide_name, output_dir, top_models):
    df_sorted = df.sort_values("model")
    fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)

    colors_hb = ["#e74c3c" if p else "#aec6e8" for p in df_sorted["passes"]]
    axes[0].bar(df_sorted["model"], df_sorted["hbond_contacts"],
                color=colors_hb, edgecolor="none", width=0.8)
    axes[0].set_ylabel("H-bond contacts\n(2.4 – 3.5 Å)", fontsize=10)
    axes[0].set_title(f"{peptide_name} — Contact Profile Across All Orientations",
                      fontsize=12, fontweight="bold")
    axes[0].axhline(y=MIN_CONTACTS, color="gray", linestyle="--", alpha=0.6)

    colors_ct = ["#2ecc71" if p else "#d5e8d4" for p in df_sorted["passes"]]
    axes[1].bar(df_sorted["model"], df_sorted["n_contacts"],
                color=colors_ct, edgecolor="none", width=0.8)
    axes[1].set_ylabel("Total contacts\n(< 5.0 Å)", fontsize=10)
    axes[1].set_xlabel("HEX Model Number", fontsize=10)

    for ax in axes:
        for m in top_models:
            ax.axvline(x=m, color="orange", linestyle=":", alpha=0.8, linewidth=1.5)

    from matplotlib.lines import Line2D
    legend_els = [
        Line2D([0], [0], color="orange", linestyle=":", linewidth=1.5,
               label=f"Top {TOP_N_EXTRACT}: {top_models}"),
        plt.Rectangle((0, 0), 1, 1, color="#e74c3c", label="Passes filter"),
        plt.Rectangle((0, 0), 1, 1, color="#aec6e8", label="Below threshold"),
    ]
    axes[0].legend(handles=legend_els, fontsize=8, loc="upper right")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "contact_profile.png"), dpi=150, bbox_inches="tight")
    plt.close()


def plot_cross_peptide_comparison(summary_df, output_dir):
    if len(summary_df) < 2:
        return

    fig    = plt.figure(figsize=(14, 10))
    gs     = gridspec.GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)
    bar_kw = dict(edgecolor="white", linewidth=0.5)
    peptides = summary_df["peptide"].tolist()
    x        = np.arange(len(peptides))

    ax1    = fig.add_subplot(gs[0, 0])
    max_hb = max(summary_df["best_hbond_contacts"].max(), 1)
    bars   = ax1.bar(x, summary_df["best_hbond_contacts"],
                     color=plt.cm.RdYlGn(summary_df["best_hbond_contacts"] / max_hb), **bar_kw)
    ax1.set_xticks(x); ax1.set_xticklabels(peptides, rotation=30, ha="right", fontsize=9)
    ax1.set_ylabel("H-bond contacts")
    ax1.set_title("Best H-bond Contacts\n(IBS residues ↔ ice)", fontweight="bold")
    for bar, val in zip(bars, summary_df["best_hbond_contacts"]):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height()+0.1,
                 str(int(val)), ha="center", va="bottom", fontsize=9)

    ax2      = fig.add_subplot(gs[0, 1])
    inv_dist = 1 / (summary_df["best_min_dist_A"] + 0.001)
    bars2    = ax2.bar(x, summary_df["best_min_dist_A"],
                       color=plt.cm.RdYlGn(inv_dist / inv_dist.max()), **bar_kw)
    ax2.axhline(y=HBOND_MAX, color="red",   linestyle="--", alpha=0.6, label=f"H-bond max ({HBOND_MAX} Å)")
    ax2.axhline(y=HBOND_MIN, color="green", linestyle="--", alpha=0.6, label=f"H-bond min ({HBOND_MIN} Å)")
    ax2.set_xticks(x); ax2.set_xticklabels(peptides, rotation=30, ha="right", fontsize=9)
    ax2.set_ylabel("Minimum distance (Å)")
    ax2.set_title("Closest IBS–Ice Contact\n(lower = tighter binding)", fontweight="bold")
    ax2.legend(fontsize=8)
    for bar, val in zip(bars2, summary_df["best_min_dist_A"]):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height()+0.02,
                 f"{val:.2f}", ha="center", va="bottom", fontsize=8)

    ax3   = fig.add_subplot(gs[1, 0])
    bars3 = ax3.bar(x, summary_df["pass_rate_pct"],
                    color=plt.cm.RdYlGn(summary_df["pass_rate_pct"]/100), **bar_kw)
    ax3.set_xticks(x); ax3.set_xticklabels(peptides, rotation=30, ha="right", fontsize=9)
    ax3.set_ylabel("% of orientations passing")
    ax3.set_title(f"Orientation Pass Rate\n(% with ≥{MIN_CONTACTS} IBS contacts)", fontweight="bold")
    ax3.set_ylim(0, 105)
    for bar, val in zip(bars3, summary_df["pass_rate_pct"]):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height()+0.5,
                 f"{val:.0f}%", ha="center", va="bottom", fontsize=9)

    ax4     = fig.add_subplot(gs[1, 1])
    metrics = summary_df[["peptide","best_hbond_contacts","best_total_contacts","pass_rate_pct"]].copy()
    metrics = metrics.set_index("peptide")
    metrics.columns = ["H-bond\nContacts", "Total\nContacts", "Pass\nRate %"]
    norm    = (metrics - metrics.min()) / (metrics.max() - metrics.min() + 1e-9)
    sns.heatmap(norm, ax=ax4, cmap="YlOrRd", annot=metrics.round(1),
                fmt="g", linewidths=0.5, cbar_kws={"label": "Normalized score"})
    ax4.set_title("Normalized Comparison\n(darker = better)", fontweight="bold")
    ax4.set_ylabel("")

    fig.suptitle("Cross-Peptide Docking Comparison — HEX Ice Slab Results",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.savefig(os.path.join(output_dir, "comparison_plot.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nComparison plot saved: {os.path.join(output_dir, 'comparison_plot.png')}")


# -----------------------------------------------------------------------------
# REPORT WRITER
# -----------------------------------------------------------------------------

def write_summary_report(summary_df, output_dir):
    report_path = os.path.join(output_dir, "summary_all_peptides.txt")
    ranked = summary_df.sort_values(
        ["best_hbond_contacts", "best_total_contacts"], ascending=False
    ).reset_index(drop=True)

    with open(report_path, "w") as fh:
        fh.write("=" * 70 + "\n")
        fh.write("  HEX DOCKING ANALYSIS — MULTI-PEPTIDE SUMMARY REPORT\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"{'Rank':<6} {'Peptide':<20} {'Best Model':<12} "
                 f"{'H-bonds':<10} {'Contacts':<10} {'Min dist':>10}  {'Pass rate'}\n")
        fh.write("-" * 75 + "\n")
        for i, row in ranked.iterrows():
            fh.write(f"  {i+1:<4} {row['peptide']:<20} {str(row['best_model']):<13} "
                     f"{int(row['best_hbond_contacts']):<10} {int(row['best_total_contacts']):<10} "
                     f"{row['best_min_dist_A']:>8.3f} Å  {row['pass_rate_pct']:.0f}%\n")
        fh.write("\n\nDETAILED RESULTS PER PEPTIDE\n" + "=" * 70 + "\n\n")
        for i, row in ranked.iterrows():
            fh.write(f"[#{i+1}] {row['peptide']}\n")
            fh.write(f"  Best orientation    : {row['best_model']}\n")
            fh.write(f"  H-bond contacts     : {int(row['best_hbond_contacts'])} (2.4–3.5 Å)\n")
            fh.write(f"  Total contacts      : {int(row['best_total_contacts'])} (< 5.0 Å)\n")
            fh.write(f"  Minimum distance    : {row['best_min_dist_A']:.3f} Å\n")
            fh.write(f"  Orientations passing: {int(row['orientations_passing'])}/{int(row['n_models'])} "
                     f"({row['pass_rate_pct']:.0f}%)\n")
            fh.write(f"  Top {TOP_N_EXTRACT} models extracted: {row['top_models']}\n")
            fh.write(f"  Ice residues found  : {row['ice_resnames_found']}\n\n")
        fh.write("\nINTERPRETATION GUIDE\n" + "=" * 70 + "\n")
        fh.write("  H-bond contacts >= 5   : Strong ice binding candidate\n")
        fh.write("  H-bond contacts  3-4   : Moderate — investigate top models\n")
        fh.write("  H-bond contacts  < 3   : Weak — deprioritize\n")
        fh.write("  Min distance < 2.4 Å   : Steric clash — check in PyMOL\n")
        fh.write("  Min distance 2.4-3.5 Å : Direct H-bond contact (ideal)\n")
        fh.write("  Min distance 3.5-5.0 Å : Van der Waals contact only\n")
        fh.write("  Min distance > 5.0 Å   : No real contact — discard\n")
        fh.write("  Pass rate > 30%        : Many favorable orientations\n")

    print(f"Report saved: {report_path}")
    return report_path


# -----------------------------------------------------------------------------
# MAIN PIPELINE
# -----------------------------------------------------------------------------

def run_pipeline():
    # Each sub-folder of HEX_DIR is one peptide
    peptide_dirs = sorted([
        d for d in glob.glob(os.path.join(HEX_DIR, "*"))
        if os.path.isdir(d)
    ])
    if not peptide_dirs:
        print(f"\nNo peptide sub-folders found in '{HEX_DIR}'")
        print("Expected layout:  results/pep42/pep42-1.pdb  results/pep42/pep42-2.pdb ...")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("\n" + "=" * 60)
    print("  HEX Docking Analysis Pipeline  (v1.2)")
    print(f"  Found {len(peptide_dirs)} peptide folder(s) in '{HEX_DIR}'")
    print("=" * 60)

    all_summaries = []
    for peptide_dir in peptide_dirs:
        peptide_name   = os.path.basename(peptide_dir)
        ibs_residues   = IBS_RESIDUE_MAP.get(peptide_name, [])
        peptide_outdir = os.path.join(OUTPUT_DIR, peptide_name)
        summary = analyze_peptide(peptide_dir, peptide_name, ibs_residues, peptide_outdir)
        if summary:
            all_summaries.append(summary)

    if not all_summaries:
        print("No results to summarise.")
        return

    summary_df = pd.DataFrame(all_summaries)
    summary_df.to_csv(os.path.join(OUTPUT_DIR, "summary_all_peptides.csv"), index=False)
    write_summary_report(summary_df, OUTPUT_DIR)
    plot_cross_peptide_comparison(summary_df, OUTPUT_DIR)

    ranked = summary_df.sort_values(
        ["best_hbond_contacts", "best_total_contacts"], ascending=False
    ).reset_index(drop=True)

    print("\n" + "=" * 60)
    print("  FINAL RANKING — ALL PEPTIDES")
    print("=" * 60)
    print(f"  {'Rank':<5} {'Peptide':<20} {'H-bonds':>8}  {'Min dist':>10}  Score")
    print("  " + "-" * 55)
    for i, row in ranked.iterrows():
        stars = "★" * min(int(row["best_hbond_contacts"]), 5)
        print(f"   {i+1:<4} {row['peptide']:<20} {int(row['best_hbond_contacts']):>8}  "
              f"{row['best_min_dist_A']:>8.3f} Å  {stars}")

    print(f"\n  Best candidate : {ranked.iloc[0]['peptide']}")
    print(f"  All results    : {OUTPUT_DIR}/\n")


if __name__ == "__main__":
    run_pipeline()
