# A python script to generate typical motif of RNA splicing_sites

## How to use
$ python generate_splicing_sites_logo.py -h
usage: generate_splicing_sites_logo.py [-h] [--gtf GTF] [--genome-fasta GENOME_FASTA]

generate consensus splicing site logos from a GTF file and a genome fasta file

options:
  -h, --help        show this help message and exit
  --gtf             an uncompressed or gzipped gtf file (default: None)
  --genome_fasta    an uncompressed or gzipped genome fasta file (default: None)
  --save_filename   an output filename for the sequence logo (default: logo_donor_site.pdf and logo_acceptor_site.pdf), the format is determined by the extension  (default: logo.pdf)


