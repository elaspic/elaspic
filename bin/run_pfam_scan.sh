#!/bin/bash

#~ python2.7 code_supl/helper/submit_job.py -c /home/kimlab1/strokach/working/pipeline/bin/run_pfam_scan.sh -i /home/kimlab1/strokach/working/databases/uniprot/filenames_01

#~ perl -MCPAN -e shell
#~ install Class::Load::XS

#~ PERL5LIB=/home/kimlab1/strokach/working/pipeline/bin/PfamScan/:$PERL5LIB /home/kimlab1/strokach/working/pipeline/bin/PfamScan/pfam_scan.pl -fasta $1 -dir /home/kimlab1/strokach/working/databases/pfam.janelia.org -outfile $1.pfamscan -clan_overlap -align -e_dom 0.0001 -e_seq 0.0001 -cpu 7


PERL5LIB=/home/kimlab1/strokach/working/pipeline/bin/PfamScan/:$PERL5LIB /home/kimlab1/strokach/working/pipeline/bin/PfamScan/pfam_scan.pl -fasta $1 -dir /home/kimlab1/strokach/working/databases/pfam -outfile $1.pfamscan -clan_overlap -align -cpu 7

