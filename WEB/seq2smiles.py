"""
seq2smiles.py
Reads an Excel file (col 1 = ID, col 2 = protein sequence, one-letter code),
builds each peptide/protein as a single molecule (peptide bonds included)
and writes canonical SMILES to a new Excel file.

Requirements:
    pip install rdkit pandas openpyxl
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def clean_sequence(raw) -> str:
    """Uppercase, strip whitespace/digits/FASTA header, keep letters only."""
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return ""
    seq = str(raw).strip()
    if seq.startswith(">"):
        seq = "".join(seq.splitlines()[1:])
    seq = re.sub(r"[\s\d\-\*]", "", seq)
    return seq.upper()


def sequence_to_smiles(seq: str) -> str:
    """Convert a one-letter protein sequence to a canonical SMILES (whole molecule)."""
    if not seq:
        raise ValueError("empty sequence")

    invalid = sorted(set(seq) - VALID_AA)
    if invalid:
        raise ValueError(f"unsupported residue(s): {', '.join(invalid)}")

    # flavor=0 -> L-amino acids, protein; peptide bonds built by RDKit
    mol = Chem.MolFromSequence(seq, flavor=0)
    if mol is None:
        raise ValueError("RDKit could not build the molecule")

    return Chem.MolToSmiles(mol, canonical=True)


def process_file(input_path: str, output_path: str | None = None, progress_callback=None) -> str:
    in_path = Path(input_path)
    if not in_path.exists():
        raise FileNotFoundError(f"File not found: {in_path}")

    df = pd.read_excel(in_path, header=0, dtype=str)
    if df.shape[1] < 2:
        raise ValueError("The Excel file needs at least two columns (ID, sequence).")

    ids = df.iloc[:, 0]
    seqs = df.iloc[:, 1]

    smiles_out, status_out = [], []
    total = len(seqs)
    for i, raw in enumerate(seqs, start=1):
        try:
            smiles_out.append(sequence_to_smiles(clean_sequence(raw)))
            status_out.append("OK")
        except Exception as e:
            smiles_out.append("")
            status_out.append(f"ERROR: {e}")
        if progress_callback:
            progress_callback(i, total)

    result = pd.DataFrame(
        {
            df.columns[0]: ids,
            df.columns[1]: seqs,
            "Canonical_SMILES": smiles_out,
            "Status": status_out,
        }
    )

    out_path = Path(output_path) if output_path else in_path.with_name(in_path.stem + "_smiles.xlsx")
    result.to_excel(out_path, index=False)
    return str(out_path)


def main():
    parser = argparse.ArgumentParser(description="Protein sequences (Excel) -> canonical SMILES (Excel)")
    parser.add_argument("input", help="Input .xlsx file")
    parser.add_argument("-o", "--output", help="Output .xlsx file (default: <input>_smiles.xlsx)")
    args = parser.parse_args()

    try:
        out = process_file(args.input, args.output)
        print(f"Done: {out}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
