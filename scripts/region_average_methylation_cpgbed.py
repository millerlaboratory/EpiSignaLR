#!/usr/bin/env python3
"""
Create region-averaged methylation matrix from per-CpG BED files produced by
gvcf_cpg_to_bed.py, instead of reading the gVCF files directly.

This produces the same output as region_average_methylation_gvcf.py. The gVCF
block rules (FILTER == "PASS", QUAL == ".") were already applied when the
per-CpG BED was created by gvcf_cpg_to_bed.py, so here each row is simply one
qualifying CpG:
    chr   pos (0-based)   pos+1   M5mC (fraction modified)   DPM5mC (depth)

From these:
    Nmod        = M5mC fraction * DPM5mC
    Nvalid_cov  = DPM5mC

Aggregation modes (per region, per sample), selected via --weighted:
    weighted: pooled read-level fraction
        region_methylation = sum(Nmod) / sum(Nvalid_cov)   over qualifying CpGs
    unweighted (default): simple mean of per-CpG fractions (Nmod / Nvalid_cov).

Note: --min-dpm5mc is applied here as well; if the CpG BED was created with a
stricter --min-dpm5mc, the stricter of the two effectively applies.

Sample names are taken from the optional second column of samples.txt, or
otherwise from the CpG BED filename with ".gz" and ".bed" extensions removed (the gVCF version uses the #CHROM header sample ID, which
is not stored in the CpG BED).

Usage:
    python region_average_methylation_cpgbed.py --samples samples.txt --bed regions.bed --output region_matrix
    python region_average_methylation_cpgbed.py --samples samples.txt --bed regions.bed --output region_matrix --min-dpm5mc 10 --weighted

Input files:
    - samples.txt: Text file with one filepath per line to per-CpG BED files
      (optionally followed by a tab and the sample name to use)
      (.bed or .bed.gz) from gvcf_cpg_to_bed.py
    - regions.bed: BED file with one row per CpG position (chr, start, end, region name),
      e.g. "chr1  1439639  1439640  DMR_1"

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


def sample_name_from_path(file_path):
    """Derive a sample name from a CpG BED filename by stripping .gz/.bed extensions."""
    name = Path(file_path).name
    for ext in ('.gz', '.bed'):
        if name.endswith(ext):
            name = name[:-len(ext)]
    return name


def read_methylation_by_region(file_path, region_index, min_dpm5mc=5, min_cpgs=1,
                               weighted=False, sample_name=None):
    """
    Read a per-CpG BED file (from gvcf_cpg_to_bed.py) and compute per-region
    methylation.

    Parameters
    ----------
    file_path : str
    region_index : dict
        Output of build_region_index.
    min_dpm5mc : int
        Minimum informative read depth (DPM5mC, column 5) to include a CpG site.
    min_cpgs : int
        Minimum number of qualifying CpG sites required to report a region value.
    weighted : bool
        If True, pool counts across CpGs: sum(Nmod) / sum(Nvalid_cov).
        If False, take the simple mean of per-CpG fractions (Nmod / Nvalid_cov).
    sample_name : str or None
        Sample name to report; if None, derived from the filename.

    Returns
    -------
    dict : region_id -> methylation fraction (0-1), or absent if < min_cpgs sites
    str  : sample name
    """
    if sample_name is None:
        sample_name = sample_name_from_path(file_path)

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
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 5:
                        continue

                    chrom = parts[0]
                    start = int(parts[1])
                    depth = float(parts[4])
                    if depth < min_dpm5mc or depth <= 0:
                        continue
                    fraction = float(parts[3])
                    nmod = fraction * depth

                    for region_id in find_regions(chrom, start, region_index):
                        region_data[region_id].append((nmod, depth))

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
                         min_dpm5mc=5, min_cpgs=1, weighted=False):
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
            min_dpm5mc=min_dpm5mc, min_cpgs=min_cpgs, weighted=weighted)

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
        description='Create region-averaged methylation matrix from per-CpG BED files '
                    '(output of gvcf_cpg_to_bed.py)')
    parser.add_argument('--samples', '-s', required=True,
                        help='Text file with paths to per-CpG BED files (one per line)')
    parser.add_argument('--bed', '-b', required=True,
                        help='BED file with one row per CpG position (chr, start, end, region name)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output file path (without extension)')
    parser.add_argument('--format', '-f', choices=['csv', 'tsv'], default='csv',
                        help='Output format (default: csv)')
    parser.add_argument('--min-dpm5mc', type=int, default=10,
                        help='Minimum informative read depth (DPM5mC) per CpG site (default: 10)')
    parser.add_argument('--min-cpgs', type=int, default=1,
                        help='Minimum qualifying CpG sites within a region to report a value (default: 1)')
    parser.add_argument('--weighted', action='store_true',
                        help='Use coverage-weighted pooled fraction (sum Nmod / sum DPM5mC) '
                             'instead of a simple mean of per-CpG fractions')

    args = parser.parse_args()

    for path, label in [(args.samples, 'Samples file'), (args.bed, 'BED file')]:
        if not os.path.exists(path):
            print(f"Error: {label} not found: {path}")
            sys.exit(1)

    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("Region Methylation Matrix Creator (CpG BED from gVCF)")
    print("=" * 40)
    print(f"Samples file:   {args.samples}")
    print(f"BED file:       {args.bed}")
    print(f"Output:         {args.output}")
    print(f"Format:         {args.format}")
    print(f"Min DPM5mC:     {args.min_dpm5mc}")
    print(f"Min CpGs:       {args.min_cpgs}")
    print(f"Aggregation:    {'coverage-weighted pooled fraction (sum Nmod / sum DPM5mC)' if args.weighted else 'simple mean of per-CpG fractions'}")
    print()

    bed_df = read_bed_regions(args.bed)
    region_index = build_region_index(bed_df)
    sample_paths = read_sample_paths(args.samples)

    if not sample_paths:
        print("Error: No valid sample files found")
        sys.exit(1)

    matrix, sample_names = create_region_matrix(
        sample_paths, bed_df, region_index,
        min_dpm5mc=args.min_dpm5mc,
        min_cpgs=args.min_cpgs,
        weighted=args.weighted)

    sep = '\t' if args.format == 'tsv' else ','
    out_file = f"{args.output}.{args.format}"
    matrix.to_csv(out_file, sep=sep)
    print(f"\nSaved region methylation matrix to {out_file}")

    summary_file = f"{args.output}_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("Region Methylation Matrix Summary (CpG BED from gVCF)\n")
        f.write("=====================================================\n\n")
        f.write(f"Matrix dimensions: {matrix.shape}\n")
        f.write(f"Number of regions: {len(matrix.index)}\n")
        f.write(f"Number of samples: {len(sample_names)}\n")
        f.write(f"Min DPM5mC:        {args.min_dpm5mc}\n")
        f.write(f"Min CpGs/region:   {args.min_cpgs}\n")
        f.write(f"Aggregation:       {'coverage-weighted pooled fraction (sum Nmod / sum DPM5mC)' if args.weighted else 'simple mean of per-CpG fractions'}\n")
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
