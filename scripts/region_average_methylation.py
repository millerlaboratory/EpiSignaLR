#!/usr/bin/env python3
"""
Create region-averaged methylation matrix from sample methylation pileups.
The BED file lists individual CpG positions (1bp intervals), with column 4
giving the region (e.g. DMR) each CpG belongs to. Multiple CpGs can share the
same region name. For each region, methylation is aggregated across only the
CpG positions listed for that region, producing a matrix with regions as rows
and samples as columns.

Aggregation modes (per region, per sample), selected via --weighted:
    weighted: pooled read-level fraction
        region_methylation = sum(Nmod) / sum(Nvalid_cov)   over qualifying CpGs
        -- each CpG contributes in proportion to its coverage; this is the
           coverage-weighted estimator and is read directly from the count
           columns (no round-trip through the rounded percent_modified field).
    unweighted (default): simple mean of per-CpG fractions (Nmod / Nvalid_cov).

Two input pileup formats are supported, selected via --platform:

    ont (default): modkit bedmethyl format
        https://github.com/nanoporetech/modkit#bedmethyl-column-descriptions
        Relevant 0-based column indices used here:
            0  chrom
            1  start (0-based)
            9  Nvalid_cov     (valid coverage)
            11 Nmod           (count of modified calls) <-- used for aggregation

    pacbio: pb-CpG-tools bedmethyl format
        https://github.com/PacificBiosciences/pb-CpG-tools
        Relevant 0-based column indices used here:
            0  chrom
            1  begin (0-based)
            5  cov             (coverage)          -> treated as Nvalid_cov
            6  est_mod_count   (est. modified count) <-- treated as Nmod

    In both cases, per-CpG coverage and modified-call count are read directly
    from the count columns; the output matrix format is identical regardless
    of which platform the input came from.

Usage:
    python region_average_methylation.py --samples samples.txt --bed regions.bed --output region_matrix --platform ont
    python region_average_methylation.py --samples samples.txt --bed regions.bed --output region_matrix --platform pacbio

Input files:
    - samples.txt: Text file with one filepath per line to bedmethyl pileup files
      (optionally followed by a tab and the sample name to use)
    - regions.bed: BED file with one row per CpG position (chr, start, end, region name),
      e.g. "chr1  1439639  1439640  DMR_1"
    - Bedmethyl files: modkit (ONT) or pb-CpG-tools (PacBio) output files with methylation data

Output:
    - CSV/TSV matrix with regions as rows and samples as columns.
      Row labels are taken from column 4 of the BED file if present,
      otherwise formatted as chr:start-end.
"""

import pandas as pd
import numpy as np
import os
import argparse
import gzip
from pathlib import Path
import sys
from collections import defaultdict


def read_bed_regions(bed_file):
    """
    Read BED file and return a DataFrame with chr, start, end, region_id.
    Region labels come from column 4 if present, else chr:start-end.
    """
    try:
        try:
            bed_df = pd.read_csv(bed_file, sep='\t', header=None,
                                 usecols=[0, 1, 2, 3],
                                 names=['chr', 'start', 'end', 'name'],
                                 dtype={0: str, 1: int, 2: int, 3: str})
        except Exception:
            bed_df = pd.read_csv(bed_file, sep='\t', header=None,
                                 usecols=[0, 1, 2],
                                 names=['chr', 'start', 'end'],
                                 dtype={0: str, 1: int, 2: int})
            bed_df['name'] = None

        if bed_df['name'].notna().all():
            bed_df['region_id'] = bed_df['name'].astype(str)
        else:
            bed_df['region_id'] = (bed_df['chr'].astype(str) + ':' +
                                   bed_df['start'].astype(str) + '-' +
                                   bed_df['end'].astype(str))

        n_regions = bed_df['region_id'].nunique()
        print(f"Loaded {len(bed_df)} CpG positions covering {n_regions} region(s) from BED file")
        print(f"Chromosomes: {sorted(bed_df['chr'].unique())}")
        return bed_df

    except Exception as e:
        print(f"Error reading BED file {bed_file}: {e}")
        sys.exit(1)


def build_region_index(bed_df):
    """
    Build a per-chromosome dict mapping 0-based position -> list of region_ids,
    one entry per BED row (each row is a single CpG position assigned to a region).
    """
    index = defaultdict(lambda: defaultdict(list))
    for row in bed_df.itertuples(index=False):
        for pos in range(row.start, row.end):
            index[row.chr][pos].append(row.region_id)
    return {chrom: dict(positions) for chrom, positions in index.items()}


def find_regions(chrom, pos, region_index):
    """
    Return the region_id(s) whose BED entry covers this exact 0-based position.
    pos is the 0-based start from the bedmethyl file.
    """
    return region_index.get(chrom, {}).get(pos, [])


def read_sample_paths(samples_file):
    try:
        with open(samples_file, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
        # Each line is "<path>" or "<path>\t<sample name>"; without a name, the
        # sample name is derived from the filename.
        valid_paths = []
        for line in lines:
            path, _, name = line.partition('\t')
            path, name = path.strip(), name.strip() or None
            if os.path.exists(path):
                valid_paths.append((path, name))
            else:
                print(f"Warning: File not found: {path}")
        print(f"Found {len(valid_paths)} valid sample files out of {len(lines)} listed")
        return valid_paths
    except Exception as e:
        print(f"Error reading samples file {samples_file}: {e}")
        sys.exit(1)


# Per-platform column layout: 0-based indices into the tab-split line for
# chrom, start, coverage, and modified-call count, plus the minimum number
# of columns a line must have for those indices to be valid.
PLATFORM_COLUMNS = {
    'ont': {'chrom': 0, 'start': 1, 'cov': 9, 'mod': 11, 'min_cols': 12},
    'pacbio': {'chrom': 0, 'start': 1, 'cov': 5, 'mod': 6, 'min_cols': 9},
}


def read_methylation_by_region(file_path, region_index, min_coverage=5, min_cpgs=1,
                               weighted=False, platform='ont', sample_name=None):
    """
    Read a bedmethyl pileup file and compute per-region methylation.

    For each qualifying CpG we read the modified-call count and the valid
    coverage directly from the count columns; the column layout depends on
    the source platform (see PLATFORM_COLUMNS):
        ont:    Nvalid_cov (col 9), Nmod (col 11) -- modkit bedmethyl output.
                percent_modified (col 10) is NOT used, avoiding a round-trip
                through its rounded value.
        pacbio: cov (col 5), est_mod_count (col 6) -- pb-CpG-tools output.

    Parameters
    ----------
    file_path : str
    region_index : dict
        Output of build_region_index.
    min_coverage : int
        Minimum coverage to include a CpG site.
    min_cpgs : int
        Minimum number of qualifying CpG sites required to report a region value.
    weighted : bool
        If True, pool counts across CpGs: sum(Nmod) / sum(Nvalid_cov).
        If False, take the simple mean of per-CpG fractions (Nmod / Nvalid_cov).
    platform : str
        'ont' (modkit bedmethyl) or 'pacbio' (pb-CpG-tools bedmethyl).
    sample_name : str or None
        Sample name to report; if None, derived from the filename.

    Returns
    -------
    dict : region_id -> methylation fraction (0-1), or absent if < min_cpgs sites
    str  : sample name
    """
    if sample_name is None:
        sample_stem = Path(file_path).stem.split('.')[0]
        sample_name = sample_stem.split('_')[0]

    cols = PLATFORM_COLUMNS[platform]

    # region_id -> list of (nmod, nvalid_cov)
    region_data = defaultdict(list)

    try:
        open_func = gzip.open if file_path.endswith('.gz') else open
        mode = 'rt' if file_path.endswith('.gz') else 'r'

        with open_func(file_path, mode) as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or not line.strip():
                    continue
                try:
                    parts = line.strip().split('\t')
                    if len(parts) < cols['min_cols']:
                        continue
                    chrom = parts[cols['chrom']]
                    start_pos = int(parts[cols['start']])
                    nvalid_cov = float(parts[cols['cov']])
                    nmod = float(parts[cols['mod']])

                    if nvalid_cov < min_coverage:
                        continue

                    for region_id in find_regions(chrom, start_pos, region_index):
                        region_data[region_id].append((nmod, nvalid_cov))

                except (ValueError, IndexError) as e:
                    print(f"Warning: Skipping malformed line {line_num} in {file_path}: {e}")
                    continue

    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return {}, sample_name

    region_averages = {}
    for region_id, values in region_data.items():
        if len(values) < min_cpgs:
            continue
        nmods = np.array([v[0] for v in values], dtype=float)
        covs = np.array([v[1] for v in values], dtype=float)
        if weighted:
            # Coverage-weighted pooled fraction: total modified / total coverage.
            region_averages[region_id] = float(nmods.sum() / covs.sum())
        else:
            # Simple mean of per-CpG fractions.
            region_averages[region_id] = float(np.mean(nmods / covs))

    print(f"  -> {sample_name}: {len(region_averages)} regions with >= {min_cpgs} qualifying CpG(s)")
    return region_averages, sample_name


def create_region_matrix(sample_paths, bed_df, region_index,
                         min_coverage=5, min_cpgs=1, weighted=False, platform='ont'):
    """
    Build a DataFrame: rows = regions (in BED order), columns = samples.
    """
    all_region_ids = list(dict.fromkeys(bed_df['region_id']))
    sample_data = {}
    sample_names = []

    print(f"Processing {len(sample_paths)} sample files...")
    for i, (file_path, name) in enumerate(sample_paths, 1):
        print(f"Processing file {i}/{len(sample_paths)}: {os.path.basename(file_path)}")
        region_avgs, sample_name = read_methylation_by_region(
            file_path, region_index, sample_name=name,
            min_coverage=min_coverage, min_cpgs=min_cpgs, weighted=weighted, platform=platform)

        if not region_avgs:
            print(f"  -> Warning: No region data found for {sample_name}, skipping")
            continue

        unique_name = sample_name
        dup_i = 1
        while unique_name in sample_data:
            dup_i += 1
            unique_name = f"{sample_name}_{dup_i}"

        sample_data[unique_name] = region_avgs
        sample_names.append(unique_name)

    if not sample_data:
        print("Error: No methylation data found in any samples")
        sys.exit(1)

    print(f"\nSuccessfully processed {len(sample_data)} samples")
    print("Building region methylation matrix...")

    matrix_data = [
        [sample_data[s].get(r, np.nan) for s in sample_names]
        for r in all_region_ids
    ]
    matrix = pd.DataFrame(matrix_data, index=all_region_ids, columns=sample_names)

    print(f"Matrix dimensions: {matrix.shape}")
    print(f"Missing values: {matrix.isna().sum().sum()}")
    print(f"Data coverage: {(1 - matrix.isna().sum().sum() / matrix.size) * 100:.1f}%")

    return matrix, sample_names


def main():
    parser = argparse.ArgumentParser(
        description='Create region-averaged methylation matrix from a bedmethyl pileup file')
    parser.add_argument('--samples', '-s', required=True,
                        help='Text file with paths to bedmethyl pileup files (one per line)')
    parser.add_argument('--bed', '-b', required=True,
                        help='BED file with one row per CpG position (chr, start, end, region name)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output file path (without extension)')
    parser.add_argument('--format', '-f', choices=['csv', 'tsv'], default='csv',
                        help='Output format (default: csv)')
    parser.add_argument('--platform', choices=['ont', 'pacbio'], default='ont',
                        help='Pileup source platform: "ont" for modkit bedmethyl output '
                             '(Nvalid_cov/Nmod columns), or "pacbio" for pb-CpG-tools bedmethyl '
                             'output (cov/est_mod_count columns). Default: ont')
    parser.add_argument('--min-coverage', type=int, default=10,
                        help='Minimum read coverage per CpG site (default: 10)')
    parser.add_argument('--min-cpgs', type=int, default=1,
                        help='Minimum qualifying CpG sites within a region to report a value (default: 1)')
    parser.add_argument('--weighted', action='store_true',
                        help='Use coverage-weighted pooled fraction (sum Nmod / sum Nvalid_cov) '
                             'instead of a simple mean of per-CpG fractions')

    args = parser.parse_args()

    for path, label in [(args.samples, 'Samples file'), (args.bed, 'BED file')]:
        if not os.path.exists(path):
            print(f"Error: {label} not found: {path}")
            sys.exit(1)

    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("Region Methylation Matrix Creator")
    print("=" * 35)
    print(f"Samples file:   {args.samples}")
    print(f"BED file:       {args.bed}")
    print(f"Output:         {args.output}")
    print(f"Format:         {args.format}")
    print(f"Platform:       {args.platform}")
    print(f"Min coverage:   {args.min_coverage}")
    print(f"Min CpGs:       {args.min_cpgs}")
    print(f"Aggregation:    {'coverage-weighted pooled fraction (sum Nmod / sum Nvalid_cov)' if args.weighted else 'simple mean of per-CpG fractions'}")
    print()

    bed_df = read_bed_regions(args.bed)
    region_index = build_region_index(bed_df)
    sample_paths = read_sample_paths(args.samples)

    if not sample_paths:
        print("Error: No valid sample files found")
        sys.exit(1)

    matrix, sample_names = create_region_matrix(
        sample_paths, bed_df, region_index,
        min_coverage=args.min_coverage,
        min_cpgs=args.min_cpgs,
        weighted=args.weighted,
        platform=args.platform)

    sep = '\t' if args.format == 'tsv' else ','
    out_file = f"{args.output}.{args.format}"
    matrix.to_csv(out_file, sep=sep)
    print(f"\nSaved region methylation matrix to {out_file}")

    summary_file = f"{args.output}_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("Region Methylation Matrix Summary\n")
        f.write("=================================\n\n")
        f.write(f"Matrix dimensions: {matrix.shape}\n")
        f.write(f"Number of regions: {len(matrix.index)}\n")
        f.write(f"Number of samples: {len(sample_names)}\n")
        f.write(f"Platform:          {args.platform}\n")
        f.write(f"Min CpG coverage:  {args.min_coverage}\n")
        f.write(f"Min CpGs/region:   {args.min_cpgs}\n")
        f.write(f"Aggregation:       {'coverage-weighted pooled fraction (sum Nmod / sum Nvalid_cov)' if args.weighted else 'simple mean of per-CpG fractions'}\n")
        f.write(f"Missing values: {matrix.isna().sum().sum()}\n")
        f.write(f"Data coverage: {(1 - matrix.isna().sum().sum() / matrix.size) * 100:.1f}%\n\n")
        f.write("Samples:\n")
        for i, s in enumerate(sample_names, 1):
            f.write(f"{i:3d}. {s}\n")
    print(f"Saved summary to {summary_file}")

    print(f"\nProcessing complete!")
    print(f"Final matrix: {matrix.shape[0]} regions x {matrix.shape[1]} samples")


if __name__ == "__main__":
    main()
