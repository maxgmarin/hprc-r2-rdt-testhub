#!/usr/bin/env python3
"""
gen_RDT_hub_v1.py
===============
Generates a UCSC trackhub directory structure for a user-supplied set of
HPRC R2 haplotype assemblies.

PRIMARY INPUT: a TSV with one row per haplotype you want in the hub.
Required columns:
    GenHapID    — e.g. HG02965_Hap2
    bb_url      — full HTTPS URL to the hosted .bb file
                  e.g. https://genome.ucsc.edu/hubspace/24/mgmarinds/.../PBKT.HG02965.H2.RDTs.bb

ASSEMBLY MAPPING: the standard HPRC ID mapping TSV
    (260507_HPRC_R2_HapAsm_IDMapping_V1.tsv)
    Used only to look up the GCA accession for each GenHapID.
    Assemblies not present in your input TSV are ignored.

OUTPUT structure (no .bb files — those are already hosted at their URLs):
    <hub_dir>/
    ├── hub.txt
    ├── genomes.txt
    ├── GCA_042031845.1/
    │   └── trackDb.txt    (bigDataUrl = full URL from your input TSV)
    ├── GCA_042032535.1/
    │   └── trackDb.txt
    ...

Usage:
    python3 gen_RDT_hub_v1.py \\
        --input    my_haplotypes.tsv \\
        --mapping  260507_HPRC_R2_HapAsm_IDMapping_V1.tsv \\
        --hub-dir  HPRC_R2_RDTs_Hub \\
        --email    marin@ds.dfci.harvard.edu

Example input TSV (my_haplotypes.tsv):
    GenHapID        bb_url
    HG02965_Hap2    https://genome.ucsc.edu/hubspace/24/mgmarinds/.../PBKT.HG02965.H2.RDTs.bb
    HG02965_Hap1    https://genome.ucsc.edu/hubspace/24/mgmarinds/.../PBKT.HG02965.H1.RDTs.bb
    HG01361_Hap1    https://genome.ucsc.edu/hubspace/24/mgmarinds/.../PBKT.HG01361.H1.RDTs.bb
"""

import argparse
import os
import sys
import pandas as pd


# ---------------------------------------------------------------------------
# Templates
# ---------------------------------------------------------------------------

HUB_TXT = """\
hub HPRC_R2_RDTs
shortLabel HPRC R2 RDTs
longLabel HPRC Release 2 Reference-Divergent Transcripts (PacBio Kinnex)
genomesFile genomes.txt
email {email}
"""

GENOME_STANZA = """\
genome {accession}
trackDb {accession}/trackDb.txt

"""

TRACKDB_STANZA = """\
track RDT_{genome_id}_{hap_short}
bigDataUrl {bb_url}
type bigBed 12 +
shortLabel {genome_id} {hap_short} RDTs
longLabel Reference-Divergent Transcripts — {genome_id} Haplotype {hap_num} (HPRC R2 Kinnex)
visibility pack
itemRgb off
searchIndex name,pbGeneId,associatedGene,cmpRef
mouseOver $associatedGene ($geneBiotype) | $classCode | FLNC reads: $flncReads | ORFanage: $orfanageStatus
url https://www.ncbi.nlm.nih.gov/nuccore/$cmpRef
urlLabel View reference transcript on NCBI RefSeq
skipEmptyFields on
"""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def hap_to_short(hap: str) -> str:
    """'Hap1' -> 'H1', 'Hap2' -> 'H2'"""
    return hap.replace("Hap", "H")


def hap_to_num(hap: str) -> str:
    """'Hap1' -> '1', 'Hap2' -> '2'"""
    return hap.replace("Hap", "")


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_inputs(input_df: pd.DataFrame, mapping_lookup: dict) -> set:
    """
    Cross-check input TSV against the assembly mapping.
    Returns the set of GenHapIDs that cannot be resolved and should be skipped.
    Exits immediately if required columns are missing.
    """
    errors = []

    for col in ("GenHapID", "bb_url"):
        if col not in input_df.columns:
            errors.append(f"Input TSV is missing required column: '{col}'")

    if errors:
        for e in errors:
            print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    missing_accession = []
    bad_url = []

    for _, row in input_df.iterrows():
        ghid = row["GenHapID"]
        url  = str(row["bb_url"]).strip()

        if ghid not in mapping_lookup:
            missing_accession.append(ghid)

        if not url.startswith("http"):
            bad_url.append(f"  {ghid}: '{url}'")

    if missing_accession:
        print(
            f"WARNING: {len(missing_accession)} GenHapID(s) not found in the assembly "
            f"mapping and will be skipped:\n  " + "\n  ".join(missing_accession),
            file=sys.stderr
        )
    if bad_url:
        print(
            f"WARNING: {len(bad_url)} bb_url value(s) don't look like HTTP(S) URLs "
            f"and will be skipped:\n" + "\n".join(bad_url),
            file=sys.stderr
        )

    skip = set(missing_accession) | {r["GenHapID"] for _, r in input_df.iterrows()
                                      if not str(r["bb_url"]).strip().startswith("http")}
    return skip


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--input", required=True,
        help="Input TSV with columns: GenHapID, bb_url"
    )
    parser.add_argument(
        "--mapping", required=True,
        help="HPRC ID mapping TSV (260507_HPRC_R2_HapAsm_IDMapping_V1.tsv)"
    )
    parser.add_argument(
        "--hub-dir", default="HPRC_R2_RDTs_Hub",
        help="Output hub directory name (default: HPRC_R2_RDTs_Hub)"
    )
    parser.add_argument(
        "--email", default="marin@ds.dfci.harvard.edu",
        help="Contact email for hub.txt"
    )
    args = parser.parse_args()

    # ── Load files ────────────────────────────────────────────────────────────
    input_df   = pd.read_csv(args.input,   sep="\t")
    mapping_df = pd.read_csv(args.mapping, sep="\t")

    for col in ("assembly_accession", "GenHapID", "haplotype", "GenomeID"):
        if col not in mapping_df.columns:
            print(f"ERROR: Assembly mapping is missing required column '{col}'",
                  file=sys.stderr)
            sys.exit(1)

    # ── Build mapping lookup: GenHapID -> full metadata row ──────────────────
    mapping_lookup = mapping_df.set_index("GenHapID")[
        ["assembly_accession", "GenomeID", "haplotype"]
    ].to_dict("index")

    # ── Validate ──────────────────────────────────────────────────────────────
    skip_ids = validate_inputs(input_df, mapping_lookup)
    input_df = input_df[~input_df["GenHapID"].isin(skip_ids)].copy()

    if input_df.empty:
        print("ERROR: No valid haplotypes remain after validation. Exiting.",
              file=sys.stderr)
        sys.exit(1)

    print(f"Building hub for {len(input_df)} haplotypes → {args.hub_dir}/")

    # ── Create hub root ───────────────────────────────────────────────────────
    os.makedirs(args.hub_dir, exist_ok=True)

    # ── hub.txt ───────────────────────────────────────────────────────────────
    with open(os.path.join(args.hub_dir, "hub.txt"), "w") as f:
        f.write(HUB_TXT.format(email=args.email))
    print("  wrote hub.txt")

    # ── genomes.txt + per-assembly trackDb.txt ────────────────────────────────
    genomes_lines = []
    n_written = 0

    for _, row in input_df.iterrows():
        genhapid = row["GenHapID"]
        bb_url   = str(row["bb_url"]).strip()

        meta      = mapping_lookup[genhapid]
        accession = meta["assembly_accession"]
        genome_id = meta["GenomeID"]
        hap       = meta["haplotype"]        # e.g. "Hap2"
        hap_short = hap_to_short(hap)        # "H2"
        hap_num   = hap_to_num(hap)          # "2"

        # Append genomes.txt stanza for this assembly
        genomes_lines.append(GENOME_STANZA.format(accession=accession))

        # Create assembly subdirectory and write trackDb.txt
        asm_dir = os.path.join(args.hub_dir, accession)
        os.makedirs(asm_dir, exist_ok=True)

        with open(os.path.join(asm_dir, "trackDb.txt"), "w") as f:
            f.write(TRACKDB_STANZA.format(
                genome_id=genome_id,
                hap_short=hap_short,
                hap_num=hap_num,
                bb_url=bb_url,
            ))
        n_written += 1

    # ── genomes.txt ───────────────────────────────────────────────────────────
    with open(os.path.join(args.hub_dir, "genomes.txt"), "w") as f:
        f.writelines(genomes_lines)

    print(f"  wrote genomes.txt  ({n_written} assemblies)")
    print(f"  wrote trackDb.txt  ({n_written} assembly subdirectories)")
    print(f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Hub structure created: {args.hub_dir}/

The directory contains only configuration files — no .bb files.
Each trackDb.txt references the full hosted URL you provided.

Upload {args.hub_dir}/ to HubSpace, then connect via:
    My Hubs → paste URL to hub.txt
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
""")


if __name__ == "__main__":
    main()
