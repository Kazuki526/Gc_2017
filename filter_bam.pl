#!/usr/bin/perl
use warnings;
use strict;

#read ABD_identical length 
my %abd_identical_length = ();
open(ABD,"ABD_identical_remove_same_seq.tsv") or die "ERROR::cant open ABD_identical_remove_same_seq.tsv";
<ABD>; # delete heeader
while(<ABD>){
		chomp;
		my @line = split(/\t/,);
		$abd_identical_length{$line[0]} = $line[3];
		$abd_identical_length{$line[1]} = $line[4];
		$abd_identical_length{$line[2]} = $line[5];
}
close ABD;

open(IN, "samtools view -h $ARGV[0]|");
open(OUT,"|samtools view -b >$ARGV[1]");
while(<IN>){
		if($_ =~ /^\@/){print OUT $_;next;}
		my @line = split(/\t/,);
		my($cigar,$length,$ref_lenght,$flag)=($line[5],length($line[9]),$abd_identical_length{$line[2]},$line[1]);
		if(($length < 100) || ($flag % 16 >= 8) || ($flag % 8 >= 4)){next;}
		if($cigar !~ /[SH]/){
				print OUT "$_";
		}elsif($line[3] == 1 ){
				my $cigar_delete_head_clip="";
				if($cigar =~/^\d+[SH](.+)$/){$cigar_delete_head_clip = $1;
				}else{$cigar_delete_head_clip ="S";}
				if($cigar_delete_head_clip !~ /[SH]/){print OUT "$_";}
		}elsif($line[3] > $ref_lenght - 100){
				my $cigar_delete_head_clip="";
				if($cigar =~/^(.+)\d+[SH]$/){$cigar_delete_head_clip = $1;
				}else{$cigar_delete_head_clip="S";}
				if($cigar_delete_head_clip !~ /[SH]/){print OUT "$_";}
		}
}
close IN;
close OUT;





