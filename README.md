# A python script to generate RNA splicing site motifs from a GTF and a genome fasta files

## How to use
$ python generate_splicing_sites_logo.py -h  
usage: generate_splicing_sites_logo.py [-h] [--gtf GTF] [--genome_fasta GENOME_FASTA]  

generate consensus splicing site logos from a GTF file and a genome fasta file  

options:  
  -h, --help        &emsp;  &emsp; show this help message and exit  
  --gtf             &emsp;  &emsp; an uncompressed or gzipped gtf file from Ensembl or Gencode (default: None)  
  --genome_fasta    &emsp;  &emsp; an uncompressed or gzipped genome fasta file (default: None)  
  --save_filename   &emsp;  &emsp; an output filename for the sequence logo (default: logo.pdf)  


