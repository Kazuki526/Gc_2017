#!/usr/bin/sh

#install nkf & BLAST & R & mafft

mkdir CS_cds
perl divide_ABD.pl Chinese_Spring/Annotation/iwgsc_refseqv1.0_HighConf_REPR_CDS_2017Apr03.fa.zip Chinese_Spring/Annotation/iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa.zip CS_cds/CS_HC
perl divide_ABD.pl Chinese_Spring/Annotation/iwgsc_refseqv1.0_LowConf_REPR_CDS_2017Apr03.fa.zip Chinese_Spring/Annotation/iwgsc_refseqv1.0_LowConf_CDS_2017Mar13.fa.zip CS_cds/CS_LC

mkdir CS_cds/DB
mkdir CS_cds/DB/HLmix
cd CS_cds/

echo "make blast db of each genomes CDS"
cat CS_HC_repr_cdsA.fa CS_LC_repr_cdsA.fa >DB/HLmix/CS_HLC_repr_cdsA.fa
cat CS_HC_repr_cdsB.fa CS_LC_repr_cdsB.fa >DB/HLmix/CS_HLC_repr_cdsB.fa
cat CS_HC_repr_cdsD.fa CS_LC_repr_cdsD.fa >DB/HLmix/CS_HLC_repr_cdsD.fa
makeblastdb -in DB/HLmix/CS_HLC_repr_cdsA.fa -out DB/HLmix/CS_HLC_repr_cdsA -dbtype nucl -parse_seqids
makeblastdb -in DB/HLmix/CS_HLC_repr_cdsB.fa -out DB/HLmix/CS_HLC_repr_cdsB -dbtype nucl -parse_seqids
makeblastdb -in DB/HLmix/CS_HLC_repr_cdsD.fa -out DB/HLmix/CS_HLC_repr_cdsD -dbtype nucl -parse_seqids

###### High Conf blast #####
#blastn -db DB/HLmix/CS_HLC_repr_cdsA -query CS_HC_repr_cdsA.fa -out blastAH2AHL.tsv -outfmt 6 -evalue 1e-6
#blastn -db DB/HLmix/CS_HLC_repr_cdsB -query CS_HC_repr_cdsB.fa -out blastBH2BHL.tsv -outfmt 6 -evalue 1e-6
#blastn -db DB/HLmix/CS_HLC_repr_cdsD -query CS_HC_repr_cdsD.fa -out blastDH2DHL.tsv -outfmt 6 -evalue 1e-6

##### High & Low Conf blast #####
echo "run blast to remove duplicated CDS"
blastn -db DB/HLmix/CS_HLC_repr_cdsA -query DB/HLmix/CS_HLC_repr_cdsA.fa -out blastAHL2AHL.tsv -outfmt 6 -evalue 1e-6
blastn -db DB/HLmix/CS_HLC_repr_cdsB -query DB/HLmix/CS_HLC_repr_cdsB.fa -out blastBHL2BHL.tsv -outfmt 6 -evalue 1e-6
blastn -db DB/HLmix/CS_HLC_repr_cdsD -query DB/HLmix/CS_HLC_repr_cdsD.fa -out blastDHL2DHL.tsv -outfmt 6 -evalue 1e-6

cd ..
mkdir singlecopy_CS_cds
##### install tidyverse package!!! #####
Rscript pick_single_copy.R

echo "make blast db of single copy CDSs & run blast between ABD each genomes"
perl blast_with_singlecopy_cds_id.pl perfectly
perl blast_with_singlecopy_cds_id.pl mostly

echo "extract ABD all matched sequence"
Rscript extract_all_aligned_cds.R perfectly
Rscript extract_all_aligned_cds.R mostly

perl check_all_aligned_cds_out_distinct_region.pl singlecopy_CS_cds/perfectly_singlecopy_all_aligned.tsv
perl check_all_aligned_cds_out_same_region.pl singlecopy_CS_cds/perfectly_singlecopy_all_aligned.tsv



