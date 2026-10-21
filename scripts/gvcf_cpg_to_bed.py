#!/usr/bin/env python3
"""
Extract per-CpG methylation calls from a single gVCF into a BED file.

Uses the same gVCF block rules as region_average_methylation_gvcf.py:
    - Only records with FILTER == "PASS" are used (other gVCF filter blocks,
      e.g. LowGQ, are skipped).
    - Only records with QUAL == "." are used (records with a real QUAL score
      are variant calls, not reference-confidence blocks/gblocks, and are
      skipped).

Within each qualifying block (POS..INFO/END), the INFO/M5mC field gives a
per-position context string for the reference allele: uppercase "Z" marks a
CpG dinucleotide with methylation reported (aggregated across both strands).
Each "Z" corresponds, in order, to one entry in the sample's FORMAT/M5mC
(fraction modified, 0-1) and FORMAT/DPM5mC (informative read depth) fields.

Output is a gzip-compressed 5-column BED file, one row per CpG:
    chr   pos (0-based)   pos+1   M5mC (fraction modified)   DPM5mC (depth)

Usage:
    python gvcf_cpg_to_bed.py --gvcf sample.gvcf --output sample_cpgs.bed.gz
    python gvcf_cpg_to_bed.py --gvcf sample.gvcf.gz --output sample_cpgs.bed.gz --min-dpm5mc 10
"""

import argparse
import gzip
import os
import sys


def parse_info_field(info_str):
    """Parse a VCF INFO field into a dict of key -> value (value is None for flags)."""
    info = {}
    for entry in info_str.split(';'):
        if '=' in entry:
            key, _, value = entry.partition('=')
            info[key] = value
        else:
            info[entry] = None
    return info


def iter_gvcf_cpgs(file_path, min_dpm5mc=0):
    """
    Yield (chrom, pos_0based, pos_1based_end, m5mc_fraction, dpm5mc_depth)
    for each qualifying CpG ("Z") in a gVCF file.

    Only PASS, QUAL=="." reference blocks are considered; within those, only
    CpGs with DPM5mC >= min_dpm5mc are yielded.
    """
    open_func = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'

    with open_func(file_path, mode) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#') or not line.strip():
                continue
            try:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 10:
                    continue

                chrom = parts[0]
                pos = int(parts[1])
                qual = parts[5]
                filt = parts[6]

                if filt != 'PASS' or qual != '.':
                    continue

                info = parse_info_field(parts[7])
                if 'M5mC' not in info or info['M5mC'] is None:
                    continue

                end = int(info['END']) if 'END' in info and info['END'] is not None else pos
                context = info['M5mC'].split(',')[0]

                z_offsets = [i for i, c in enumerate(context) if c == 'Z']
                if not z_offsets:
                    continue

                if len(context) != end - pos + 1:
                    print(f"Warning: Skipping line {line_num} in {file_path}: "
                          f"M5mC context length ({len(context)}) does not match "
                          f"block length ({end - pos + 1})", file=sys.stderr)
                    continue

                format_keys = parts[8].split(':')
                sample_values = parts[9].split(':')
                if 'M5mC' not in format_keys or 'DPM5mC' not in format_keys:
                    continue
                m5mc_str = sample_values[format_keys.index('M5mC')]
                dpm5mc_str = sample_values[format_keys.index('DPM5mC')]
                if m5mc_str in ('.', '') or dpm5mc_str in ('.', ''):
                    continue

                fractions = m5mc_str.split(',')
                depths = dpm5mc_str.split(',')
                if len(fractions) != len(z_offsets) or len(depths) != len(z_offsets):
                    print(f"Warning: Skipping line {line_num} in {file_path}: "
                          f"M5mC/DPM5mC value count does not match number of "
                          f"CpGs marked in the INFO/M5mC context string", file=sys.stderr)
                    continue

                for offset, frac_str, depth_str in zip(z_offsets, fractions, depths):
                    if frac_str == '.' or depth_str == '.':
                        continue
                    depth = float(depth_str)
                    if depth < min_dpm5mc:
                        continue
                    fraction = float(frac_str)

                    genomic_pos_1based = pos + offset
                    bed_start_0based = genomic_pos_1based - 1

                    yield chrom, bed_start_0based, bed_start_0based + 1, fraction, depth

            except (ValueError, IndexError) as e:
                print(f"Warning: Skipping malformed line {line_num} in {file_path}: {e}",
                      file=sys.stderr)
                continue


def main():
    parser = argparse.ArgumentParser(
        description='Extract per-CpG methylation calls from a gVCF into a BED file')
    parser.add_argument('--gvcf', '-g', required=True,
                        help='Input gVCF file (.gvcf or .gvcf.gz)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output BED file path (gzip-compressed; ".gz" appended if missing)')
    parser.add_argument('--min-dpm5mc', type=int, default=0,
                        help='Minimum informative read depth (FORMAT/DPM5mC) to include a '
                             'CpG site (default: 0, no filtering)')
    parser.add_argument('--no-header', action='store_true',
                        help='Do not write the "#chr..." header comment line')

    args = parser.parse_args()

    if not os.path.exists(args.gvcf):
        print(f"Error: gVCF file not found: {args.gvcf}")
        sys.exit(1)

    output_path = args.output if args.output.endswith('.gz') else args.output + '.gz'
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    n_written = 0
    with gzip.open(output_path, 'wt') as out:
        if not args.no_header:
            out.write("#chr\tpos\tpos+1\tM5mC\tDPM5mC\n")
        for chrom, start, end, fraction, depth in iter_gvcf_cpgs(args.gvcf, min_dpm5mc=args.min_dpm5mc):
            out.write(f"{chrom}\t{start}\t{end}\t{fraction:g}\t{depth:g}\n")
            n_written += 1

    print(f"Wrote {n_written} CpG position(s) to {output_path}")


if __name__ == "__main__":
    main()
