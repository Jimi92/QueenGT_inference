#!/usr/bin/python

"""
Uses haploid drones' genotypes to recreate queen genotypes.

The script reads the VCF that contain all the drones' genotypes (no index file required), and an accompanying pedigree file, 
to differentiate drones of the same colonies, and treat them as reads to call queen genotypes. Drones VCF should be filtered 
prior to using this script with GATK VariantFiltration to ensure quality of each positions.  

Args:
@parser@

Results:
    vcf.gz
        For each colony, a vcf.gz with the queen's genotype is output. 
        
"""

import argparse
import gzip
import subprocess
import pandas as pd
import pysam
from collections import defaultdict
from tqdm import tqdm
import os
import sys
from concurrent.futures import ProcessPoolExecutor

def parse_args():
    parser = argparse.ArgumentParser(description="Build queen's genotype based on drones' genotypes.")
    parser.add_argument('--vcf',
                    type=str, required=True, help='Address of the input vcf.gz file.')
    parser.add_argument('--list',
                    type=str, required=True, help='Address of the input list of queens and drones.')
    parser.add_argument('--threads',
                    type=int, default=1, help='Number of threads to use (default: 1)')
    parser.add_argument('--out',
                    type=str, required=True, help='Address of the output folder.')
    parser.add_argument('--min_drones', 
                    type=int, default=2, help='Minimum number of drone genotypes required per site (default: 2)')
    parser.add_argument('--het_thres', 
                    type=float, default=0.125, help='Threshold for heterozygosity (default: 0.125)')
    parser.add_argument('--hap', 
                    dest='haploid', action='store_true', help='Invoke if input genotypes are haploid')
    
    return parser.parse_args()

def open_maybe_gz(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode=mode, encoding='utf-8')
    else:
        return open(filename, mode=mode, encoding='utf-8')

def extract_vital_header_lines(vcf_file):
    vital_headers = []
    with open_maybe_gz(vcf_file) as file:
        for line in file:
            if line.startswith('##FILTER') or line.startswith('##ALT') or line.startswith('##FORMAT') or \
               line.startswith('##INFO') or line.startswith('##CONTIG'):
                vital_headers.append(line)
            if line.startswith('#CHROM'):
                vital_headers.append(line)  # Include the #CHROM line
                break  # Stop reading after the header is complete
    return vital_headers

def get_queen_genotypes(input_vcf, queen_name, drones, threshold, min_drones, haploid=False):
    genotypes = {}  # variant_key -> queen_gt
    info_map = {}   # variant_key -> info_cols

    with open_maybe_gz(input_vcf) as f:
        header_line = None
        for line in f:
            if line.startswith('#CHROM'):
                header_line = line.strip()
                break

        if not header_line:
            raise ValueError("No #CHROM line found in VCF")

        cols = header_line.lstrip('#').split('\t')
        sample_names = cols[9:]

        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            info_cols = parts[:9]
            sample_data = dict(zip(sample_names, parts[9:]))
            variant_key = '\t'.join(parts[:5])  # CHROM, POS, ID, REF, ALT

            drone_gts = []
            for d in drones:
                gt = sample_data.get(d, './.')
                gt = gt.split(':')[0] if gt != '.' else './.'
                drone_gts.append(gt)

            queen_gt = infer_queen_gt_from_list(drone_gts, threshold, min_drones, haploid=haploid)
            genotypes[variant_key] = queen_gt
            info_map[variant_key] = info_cols

    return queen_name, genotypes, info_map

def infer_queen_gt_from_list(drone_gts, threshold, min_drones, haploid=False):
    # Filter valid drone genotypes (ignore missing)
    if haploid:
        valid_gts = [gt for gt in drone_gts if gt not in ['.', '',]]
    else:
        valid_gts = [gt for gt in drone_gts if gt not in ['./.', '',]]
    total = len(valid_gts)
    if total < min_drones:
        return './.'  # Not enough data

    # Count how many have alternate allele (anything other than '0')
    if haploid:
        alt_count = sum(1 for gt in valid_gts if gt != '0')
    else:
        alt_count = sum(1 for gt in valid_gts if gt != '0/0')

    ratio = alt_count / total

    # Infer queen diploid genotype based on ratio
    if ratio < threshold:
        return '0/0'    # Mostly reference
    elif ratio > (1 - threshold):
        return '1/1'    # Mostly alternate
    else:
        return '0/1'    # Mixed alleles (heterozygous)

def extract_samples_from_header(vcf_file):
    with open_maybe_gz(vcf_file) as f:
        for line in f:
            if line.startswith('#CHROM'):
                header_cols = line.lstrip('#').strip().split('\t')
                # Samples start from the 10th column (index 9)
                return header_cols[9:]
    return []

def variant_key_sorter(key):
    chrom, pos, *_ = key.split('\t')
    pos = int(pos)
    return (chrom, pos)

def write_combined_vcf(output_path, input_vcf, queen_names, all_genotypes, info_map):
    temp_vcf = output_path.rstrip('.gz')
    variant_keys = sorted(all_genotypes.keys(), key=variant_key_sorter)

    with open(temp_vcf, 'w') as out:
        # Write headers
        with open_maybe_gz(input_vcf) as f:
            for line in f:
                if line.startswith('##'):
                    out.write(line)
                elif line.startswith('#CHROM'):
                    fixed_cols = line.strip().split('\t')[:9]
                    out.write('\t'.join(fixed_cols + queen_names) + '\n')
                    break

        for variant_key in variant_keys:
            info_cols = info_map[variant_key]
            queen_gts = [all_genotypes[variant_key].get(q, './.') for q in queen_names]
            out.write('\t'.join(info_cols + queen_gts) + '\n')

    subprocess.run(['bgzip', '-f', temp_vcf], check=True)
    subprocess.run(['tabix', '-p', 'vcf', output_path], check=True)

def main():
    if '--version' in sys.argv:
        print("Haplodiploid Queen Genotype Inference v1.0.0\n"
              "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
              "This is free software: you are free to change and redistribute it.\n"
              "There is NO WARRANTY, to the extent permitted by law.")
        return
    
    args = parse_args()

    queen_map = defaultdict(list)
    with open(args.list, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                queen, drone = parts[0], parts[1]
                queen_map[queen].append(drone)
            else:
                print(f"Bad line in list file: {line.strip()}")

    vcf_samples = extract_samples_from_header(args.vcf)

    for queen, drones in queen_map.items():
        missing = [d for d in drones if d not in vcf_samples]
        if missing:
            print(f"Warning: For queen {queen}, missing drones in VCF: {', '.join(missing)}")

    os.makedirs(args.out, exist_ok=True)

    all_genotypes = defaultdict(dict)  # variant_key -> {queen1: gt1, ...}
    shared_info_map = {}

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for queen, drones in queen_map.items():
            futures.append(executor.submit(get_queen_genotypes, args.vcf, queen, drones, args.het_thres, args.min_drones, args.haploid))

        for future in tqdm(futures, desc="Processing queens"):
            queen_name, queen_gts, info_map = future.result()
            for variant_key, gt in queen_gts.items():
                all_genotypes[variant_key][queen_name] = gt
                shared_info_map[variant_key] = info_map[variant_key]

    output_vcf_path = os.path.join(args.out, "all_queens.vcf.gz")
    write_combined_vcf(output_vcf_path, args.vcf, list(queen_map.keys()), all_genotypes, shared_info_map)

if __name__ == "__main__":
    main()
