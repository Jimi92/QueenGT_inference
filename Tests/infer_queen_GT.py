import pandas as pd
import argparse
from tqdm import tqdm
from collections import defaultdict

# Argument Parser
parser = argparse.ArgumentParser(description='Infers queen genotype using the genotype of drone sons.')
parser.add_argument('-v', '--vcf', type=str, help='Input VCF file', required=True)
parser.add_argument('-l', '--list', type=str, help='List of individual ID and queen ID', required=True)
parser.add_argument('-r', '--rownum', type=int, help='Number of header rows in the VCF file minus 1', required=True)
parser.add_argument('-ht', '--heterozygot_thres', type=float, default=0.125, help='Threshold for heterozygosity (default=0.125)')
parser.add_argument('--haploid', action='store_true', help='Use this flag if the input VCF contains haploid genotypes')

args = parser.parse_args()

# ---------------------- PARSE INPUT FILES ---------------------------------
# Read header lines including the Sample line
with open(args.vcf, 'r') as vcf_file:
    lines = vcf_file.readlines()

header_lines = lines[:args.rownum + 1]
column_line = header_lines[-1].strip().lstrip('#')
column_names = column_line.split('\t')

if len(column_names) != len(set(column_names)):
    duplicates = [x for x in column_names if column_names.count(x) > 1]
    print(f"Warning: Duplicate column names found in VCF header: {set(duplicates)}")

# Read VCF data
vcf = pd.read_csv(args.vcf, sep='\t', skiprows=args.rownum + 1, dtype=str)

# Read haploid individual list
droneID = []
queenID = []
with open(args.list, 'r') as file:
    for line in file:
        parts = line.strip().split()
        if len(parts) >= 2:
            droneID.append(parts[0].strip())
            queenID.append(parts[1].strip())
        else:
            print('Oops: Something is wrong with the list file.')

# Check samples in VCF match the list
vcf_samples = list(vcf.columns[9:])

print("From list:", droneID)
print("From VCF :", vcf_samples)

if droneID != vcf_samples:
    print("Oops: Individuals in the list do not match the VCF column names.")
else:
    print(f"Individuals matched between list and VCF: {', '.join(droneID)}")

# ----------------- FUNCTION TO INFER QUEEN GENOTYPES ----------------------
