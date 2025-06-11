# QueenGT_inference

### A tool for inferring bee queen's genotype based on the genotype of her haploid sons (aka the drones) 

-----------------------------------------------------------------------------------
## Required input

 * A VCF file

 * A map of droneID and the corresponding queenID. One individual per line.
  
 * The number of header lines in the vcf file.
  
## Usage

[usage]: python infer_queen_GT2_latest.py -v your.vcf.vcf -l ind2fam_map.txt -r 231 

options:
  -v VCF, --vcf VCF     Input VCF file
  
  -l LIST, --list LIST  List of individual ID and queen ID
  
  -r ROWNUM, --rownum ROWNUM Number of header rows in the VCF file minus 1
  
  -ht HETEROZYGOT_THRES, --heterozygot_thres HETEROZYGOT_THRES Threshold for heterozygosity (default=0.125)
  
--------------------------------------------------------------------------------
 
## Notes

By default the script considers that the queen is heterozygous if at least 0.125 of her ofspring carry the minor allele. This threshold can change using the -ht flag

If there is information by only three or less drones, the script sets queens genotype to missing on that possition

The script can handle up to two alleles.
