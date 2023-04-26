import pandas as pd
import os
import shutil
from tqdm import tqdm
pd.options.display.width = 0
sra_run_info = 'SraRunTable.txt'
kallisto_idx_path = 'Spenn-v2-cds-annot.fa.idx'
sra_data = pd.read_csv(sra_run_info)
sra_data.set_index('Run', inplace=True)
print(sra_data.head())
print('--------------------------------------------------------------------------------------------------\n')
print(f'Number of SRA runs is {sra_data.shape[0]}')
print('Starting processing')
print('--------------------------------------------------------------------------------------------------\n')

if not os.path.exists('Spenn'):
    os.system('mkdir Spenn')

for sra_run in tqdm(sra_data.index):
    if not os.path.exists(f'{sra_run}.fastq'):
        print(f"Run {sra_run} : Layout {sra_data.loc[sra_run, 'LibraryLayout']}")
        if not os.path.exists(f'{sra_run}'):
            os.system(f"fasterq-dump {sra_run} -e 10")

            # For SINGLE END Layouts
            if sra_data.loc[sra_run, 'LibraryLayout'] == 'SINGLE':
                os.system(f"sickle se -f {sra_run}.fastq -t sanger -o {sra_run}_trimm.fastq")
                os.system(f"kallisto quant -i {kallisto_idx_path} -o {sra_run} --single -l 200 -s 20 {sra_run}_trimm.fastq")
            # For PAIRED END Layouts
            elif sra_data.loc[sra_run, 'LibraryLayout'] == 'PAIRED':
                os.system(f"sickle pe -f {sra_run}_1.fastq -r {sra_run}_2.fastq -t sanger -o {sra_run}_1_trimm.fastq -p {sra_run}_2_trimm.fastq -s {sra_run}_singles.fastq")
                os.system(f"kallisto quant -i {kallisto_idx_path} -o {sra_run} {sra_run}_1_trimm.fastq {sra_run}_2_trimm.fastq")

            os.system(f"rm -rf *.fastq")


for file in os.listdir():
    if file.startswith('SRR'):
        if os.path.exists(f'{file}/abundance.tsv'):
            shutil.move(file, 'Spenn')
        else:
            os.system(f'rm -rf {file}')
