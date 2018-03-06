#!/user/bin/perl
use strict;
use warnings;

my $replicate_ref = "454-RNA_sequencing/T_aes_Cluster_AllTranscripts.fasta";
my $ident_ref = "ABD_identical_fa/ABD_identical.fa";
if((!-e $replicate_ref)||(!-e $ident_ref)){die "ref_file isnot exist!!\n";}
my ($seqdir,$replicatedir,$ABDidenticaldir)=("2017_RNA-Seq_Gc2","mizuno_replication","ABD_identical");
mkdir $replicatedir;
mkdir $ABDidenticaldir;
#system("bwa index $replicate_ref");
#system("bwa index $ident_ref");
open(IN,"$seqdir/file_list.txt");
while(<IN>){
		chomp;
		my$sample=$_;
		my($fastaq1,$fastaq2)=("$seqdir/$sample"."1.fastq.gz","$seqdir/$sample"."2.fastq.gz");
		my($pairfaq1,$pairfaq2)=("$seqdir/$sample"."pair_1.fastq","$seqdir/$sample"."pair_2.fastq");
		system("trimmomatic PE -phred33 $fastaq1 $fastaq2 $pairfaq1 /dev/null $pairfaq2 /dev/null ILLUMINACLIP:$seqdir/adaptar.fa:2:20:10 SLIDINGWINDOW:4:25 LEADING:20 TRAILING:20 MINLEN:100");
		my($replibam,$identbam,$identbam_filter)=("$replicatedir/$sample"."to_T_aes_Cluster_AllTranscripts.bam",
												  "$ABDidenticaldir/$sample"."to_ABD_identical.bam",
												  "$ABDidenticaldir/$sample"."to_ABD_identical_filter.bam");
		#system("bwa mem -M $replicate_ref $pairfaq1 $pairfaq2 |samtools sort >$replibam");
		system("bwa mem -M $ident_ref $pairfaq1 $pairfaq2 |samtools sort >$identbam");
		system("perl ~/git/Gc_2017/filter_bam.pl $identbam $identbam_filter");
}
close IN;
exit;
