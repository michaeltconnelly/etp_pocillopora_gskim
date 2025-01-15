# python script to summarise blast hits 
# modified for loci as an argument by Mike Connelly, 1/15/25

import pandas as pd
import os
import argparse
from pathlib import Path
import shutil

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Summarise BLAST hits for a given locus.")
parser.add_argument("locus", type=str, help="Locus name (e.g., mtcob, ITS2, cp23S, psbA).")
args = parser.parse_args()

# Use the provided locus
locus = args.locus

# Directory that contains BLAST hits across samples
txt_dir = f"/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/noncoral_spades/{locus}_blast_hits"

# Output text path
out_txt = f"/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim/outputs/noncoral_spades/{locus}_blast_hits_summary.csv"

# Change to the specified directory
os.chdir(txt_dir)

marker_hit_num = []
seq_hit_list = []

for file in os.listdir():
    with open(file) as f:
        sample_list = []
        seq_hit = []
        if str(file).endswith("_contigs"):
            pass # Skip files ending with _contigs
        else: #looking for contigs.sorted
            print(str(file))
            sample_name = file.split("_contigs")[0]
            print(sample_name)
            line = f.readline().strip('\n')
            items = line.split('	')
            sample_list.append(sample_name)

            if len(items) > 1:
                seq_hit.append(sample_name)
                seq_hit.append(items[1])
                seq_hit_list.append(seq_hit)

            sample_list = sample_list + items
            marker_hit_num.append(sample_list)

# Convert into dataframe
df = pd.DataFrame(data=marker_hit_num)
print(df)

# Convert into excel
df.to_csv(out_txt, index=False)



