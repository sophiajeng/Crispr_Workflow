## Crispr Simplified Protocol
Library annotation files:
Using the Bioconductor package Biostrings, pattern_match.R goes through each sgRNA sequence and identifies the number of other sgRNA sequences that differ by 1 nucleotide.  pattern_match.R takes in four arguments
--lib_file: the csv file with sgRNA name, sequence, and gene
--fasta_file: fasta file for Crispr library
--num_mismatch: maximum number of mismatches 
--lib_name: library name
**NOTE** The resulting file will be named <lib_name>_sequence_sim_le_<num_mismatch>.txt
**NOTE** This needs to be performed just ONCE per LIBRARY.
Example:  
```bash
nohup srun -A mcweeney_lab --job-name=brunello --time=1440 -p exacloud --mail-type=END --mail-user=jengs@ohsu.edu --mem=12G R CMD BATCH --vanilla '--args --lib_file=/home/exacloud/gscratch/mcweeney_lab/resources/u54/crispr/libraries/brunello_library1.csv --fasta_file=/home/exacloud/gscratch/mcweeney_lab/resources/u54/crispr/libraries/brunello_library.fasta --num_mismatch=1 --lib_name=brunello' pattern_match.R > pattern_match.out &
```
This script has been run for the yusa, brunello, mousa, and eil20/eil21 libraries.

### Crispr Pipeline
1. Fastq trimming using cutadapt. It will take in as input the 5 prime and 3 prime adapter sequences found in the library_parameters.txt file
2. Use python script,crispr_counts.py,  to count number of occurrences of each sgRNA.
3. basicRPython.R. This produces the current_cell_line_crispr_counts.txt and current_cell_line_crispr_stats.txt files.

```bash
source /home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/activate crispr
```

To create a crispr conda environment:
```bash
/home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/conda create -n crispr
source /home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/activate crispr

#software from previous Crispr implementation
#/home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/conda config --add channels defaults
#/home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/conda config --add channels bioconda
/home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/conda config --add channels conda-forge
#bowtie2
#/home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/conda install -c bioconda bowtie2
#bamtools
#/home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/conda install -c bioconda bamtools
#samtools
#/home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/conda install -c bioconda samtools
#cutadapt
/home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/conda install -c bioconda cutadapt
#mageck
#/home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/conda install -c bioconda mageck
```

To display version of software installed
```bash
/home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/conda list crispr
```

The python_crispr.wdl takes in as input crispr_python_inputs.json:
```bash
{
    
    "crispr_python_workflow.sample_name_file": "sample_name.txt",
    "crispr_python_workflow.library_file": "library.txt",
    "crispr_python_workflow.library_config_file": "library_parameters.txt"
    "crispr_python_workflow.annot_file": "/home/groups/mcweeney_lab/jengs/crispr_library_scores/eil_20_21_sequence_sim_le_1.txt",
    "crispr_python_workflow.script_dir": "/home/groups/mcweeney_lab/jengs/Crispr"

}   
```

1. 'library.txt' is a single column file with the library that is used for this run
2. 'library_parameters.txt' is a tab-delimited file containing library file (csv file), 5 and 3 prime adapters
3. 'sample_name.txt' contains two columns with the first column being location of the fastq file and the second column the sample name
4. 'annot_file' is the annotation file produced by the pattern_match.R script
5. 'script_dir' is the directory where the scripts are (basicRPython.R and crispr_counts.py)

```bash
sbatch -A mcweeney_lab --mail-type END --mail-user jengs@ohsu.edu -o crispr.out -e crispr.err -t 2160 -p exacloud -c 2 --mem=5G  \
--wrap "
source /home/exacloud/gscratch/mcweeney_lab/resources/programs/miniconda2/bin/activate crispr

java -jar -Djava.io.tmpdir=tmp_dir \
-Dconfig.file=/home/exacloud/gscratch/mcweeney_lab/resources/programs/exacloud.cromwell.conf /home/exacloud/gscratch/mcweeney_lab/resources/programs/cromwell-48.jar run \
python_crispr.wdl \
-i crispr_python_inputs.json \
-m crispr.json
"
```

