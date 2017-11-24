#!/usr/bin/perl
use warnings;
use strict;

#extracting representative gene UTR+CDS by High and Low conf gff

#directri check
my @ls=`ls`;chomp @ls;
if(! grep{$_ eq "Chinese_Spring"}@ls){die "ERROR::starting on wrong directri\n";}

my $ann_dir="Chinese_Spring/Annotation/";
my $ge_dir ="Chinese_Spring/Genome/";
#high conf gff and representative cds and all cds
my $hc_gff="iwgsc_refseqv1.0_HighConf_2017Mar13.gff3.zip";
my $hc_cds ="iwgsc_refseqv1.0_HighConf_REPR_CDS_2017Apr03.fa.zip";
my $hc_cds_all = "iwgsc_refseqv1.0_HighConf_CDS_2017Apr03.fa.zip";
#low conf gff and representative cds and all cds
my $lc_gff="iwgsc_refseqv1.0_LowConf_2017Mar13.gff3.zip";
my $lc_cds="iwgsc_refseqv1.0_LowConf_REPR_CDS_2017Apr03.fa.zip";
my $lc_cds_all="iwgsc_refseqv1.0_LowConf_CDS_2017Apr03.fa.zip";

my $gene_name_forcal=0;
#high confidence
my %repr=&pick_representative_trid("$ann_dir$hc_cds");
my %hc_gff=&extract_cdsutr_region_from_gff("$ann_dir$hc_gff");

#low confidence
%repr=&pick_representative_trid("$ann_dir$lc_cds");
my %lc_gff=&extract_cdsutr_region_from_gff("$ann_dir$lc_gff");


@ls=`ls $ge_dir`;chomp @ls;
my @chr_fa=grep{$_=~/^iwgsc_refseqv1.0_chr\d[ABD]\.fsa.zip$/}@ls;
my %hc_cds=&fasta2hash("$ann_dir$hc_cds_all");
my %lc_cds=&fasta2hash("$ann_dir$lc_cds_all");

my $out_dir="CS_utrcds";
open(OUTHA,"|gzip -c >$out_dir/CSA_hc_repr.fa.gz");
open(OUTHB,"|gzip -c >$out_dir/CSB_hc_repr.fa.gz");
open(OUTHD,"|gzip -c >$out_dir/CSD_hc_repr.fa.gz");
open(OUTHAO,"|gzip -c >$out_dir/CSA_hc_notrepr.fa.gz");
open(OUTHBO,"|gzip -c >$out_dir/CSB_hc_notrepr.fa.gz");
open(OUTHDO,"|gzip -c >$out_dir/CSD_hc_notrepr.fa.gz");
open(OUTLA,"|gzip -c >$out_dir/CSA_lc_repr.fa.gz");
open(OUTLB,"|gzip -c >$out_dir/CSB_lc_repr.fa.gz");
open(OUTLD,"|gzip -c >$out_dir/CSD_lc_repr.fa.gz");
open(OUTLAO,"|gzip -c >$out_dir/CSA_lc_notrepr.fa.gz");
open(OUTLBO,"|gzip -c >$out_dir/CSB_lc_notrepr.fa.gz");
open(OUTLDO,"|gzip -c >$out_dir/CSD_lc_notrepr.fa.gz");
foreach my $chr_fa(@chr_fa){
		my %fasta=&fasta2hash("$ge_dir$chr_fa");
		if(scalar(keys %fasta) !=1){die "ERROR::$chr_fa have more than one fasta file\n".keys(%fasta)."\n";}
		my @chr=keys %fasta;my$chr=$chr[0];
		if(scalar(keys(%{$hc_gff{$chr}}))==0){die "$chr hc_gff have no data?\n";}else{print "doing $chr\n";}
		foreach my $tr_id(keys %{$hc_gff{$chr}}){
				print "$tr_id\n";
				my ($strand,$cds_inf)=($hc_gff{$chr}{$tr_id}{strand},$hc_gff{$chr}{$tr_id}{cds});
				my ($utr5,$utr5_inf,$utr3,$utr3_inf)=("","","","");
				if(defined $hc_gff{$chr}{$tr_id}{'5utr_inf'}){$utr5_inf=$hc_gff{$chr}{$tr_id}{'5utr_inf'};}
				if(defined $hc_gff{$chr}{$tr_id}{'3utr_inf'}){$utr3_inf=$hc_gff{$chr}{$tr_id}{'3utr_inf'};}
				my $cds=$hc_cds{$tr_id};
				if(defined $hc_gff{$chr}{$tr_id}{'5utr'}){
						my @utr=split(/;/,$hc_gff{$chr}{$tr_id}{'5utr'});
						foreach my $utr(@utr){
								my ($start,$end)=split(/-/,$utr);
								if($strand eq "+"){
										$utr5.=substr($fasta{$chr},$start-1,$end-1);
								}else{
										$utr5.= &reverse_trans(substr($fasta{$chr},$start-1,$end-1));
								}}
				}elsif(defined $hc_gff{$chr}{$tr_id}{'3utr'}){
						my @utr=split(/;/,$hc_gff{$chr}{$tr_id}{'3utr'});
						foreach my $utr (@utr){
								my ($start,$end)=split(/-/,$utr);
								if($strand eq "+"){
										$utr3.=substr($fasta{$chr},$start-1,$end-1);
								}else{
										$utr3.= &reverse_trans(substr($fasta{$chr},$start-1,$end-1));
								}}
				}
				if(defined $hc_gff{$chr}{$tr_id}{repr}){
						if($chr=~/A$/){
								print OUTHA ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
								print ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
						}elsif($chr=~/B$/){
								print OUTHB ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
								print ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
						}elsif($chr=~/D$/){
								print OUTHD ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
								print ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
						}else{die "ERROR::what chr?? $chr\n";}
				}else{
						if($chr=~/A$/){
								print OUTHAO ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
								print ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
						}elsif($chr=~/B$/){
								print OUTHBO ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
								print ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
						}elsif($chr=~/D$/){
								print OUTHDO ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
								print ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
						}
				}die"finish";
		}
		foreach my $tr_id(keys %{$lc_gff{$chr}}){
				my ($strand,$cds_inf)=($lc_gff{$chr}{$tr_id}{strand},$lc_gff{$chr}{$tr_id}{cds});
				my ($utr5,$utr5_inf,$utr3,$utr3_inf)=("","","","");
				if(defined $lc_gff{$chr}{$tr_id}{'5utr_inf'}){$utr5_inf=$lc_gff{$chr}{$tr_id}{'5utr_inf'};}
				if(defined $lc_gff{$chr}{$tr_id}{'3utr_inf'}){$utr3_inf=$lc_gff{$chr}{$tr_id}{'3utr_inf'};}
				my $cds=$lc_cds{$tr_id};
				if(defined $lc_gff{$chr}{$tr_id}{'5utr'}){
						my @utr=split(/;/,$lc_gff{$chr}{$tr_id}{'5utr'});
						foreach my $utr(@utr){
								my ($start,$end)=split(/-/,$utr);
								if($strand eq "+"){
										$utr5.=substr($fasta{$chr},$start-1,$end-1);
								}else{
										$utr5.= &reverse_trans(substr($fasta{$chr},$start-1,$end-1));
								}}
				}elsif(defined $lc_gff{$chr}{$tr_id}{'3utr'}){
						my @utr=split(/;/,$lc_gff{$chr}{$tr_id}{'3utr'});
						foreach my $utr (@utr){
								my ($start,$end)=split(/-/,$utr);
								if($strand eq "+"){
										$utr3.=substr($fasta{$chr},$start-1,$end-1);
								}else{
										$utr3.= &reverse_trans(substr($fasta{$chr},$start-1,$end-1));
								}}
				}
				if(defined $lc_gff{$chr}{$tr_id}{repr}){
						if($chr=~/A$/){
								print OUTLA ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
						}elsif($chr=~/B$/){
								print OUTLB ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
						}elsif($chr=~/D$/){
								print OUTLD ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
						}
				}else{
						if($chr=~/A$/){
								print OUTLAO ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
						}elsif($chr=~/B$/){
								print OUTLBO ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
						}elsif($chr=~/D$/){
								print OUTLDO ">$tr_id $strand;$utr5_inf;$cds_inf;$utr3_inf;\n$utr5$cds$utr3\n"; 
						}
				}
		}

}

close OUTHA;
close OUTHB;
close OUTHD;
close OUTHAO;
close OUTHBO;
close OUTHDO;
close OUTLA;
close OUTLB;
close OUTLD;
close OUTLAO;
close OUTLBO;
close OUTLDO;


###############################################################################################
sub reverse_trans($){
		my $seq=$_[0];
		$seq =~ tr/AaTtGgCc/TtAaGgCc/;
		return(reverse($seq));
}


###############################################################################################
sub pick_representative_trid( $ ){
		my $repr_path=$_[0];
		my %repr=();
		open(FA,"unzip -p $repr_path|");
		while(<FA>){
				chomp;
				if($_=~/^>(\S+) gene=\S+$/){
						$repr{$1}="yes";
				}elsif($_=~/[^ATGCNatgcn]/){
						die "ERROR::$repr_path have wrong fasta file? what this line?\n$_\n";
				}
		}
		close FA;
		return(%repr);
}

################################################################################################
sub extract_cdsutr_region_from_gff( $ ){
		my $gff_path=$_[0];
		open(GFF,"unzip -p $gff_path|");
		my %gff=();
		my $l=0;
		while(<GFF>){
				chomp;$l++;
				my @line=split(/\t/,);
				my %info=map{my @inf=split(/=/,$_);($inf[0],$inf[1])}split(/;/,$line[8]);
				if($line[2] eq "exon" || $line[2] eq "gene"){next;
				}elsif($line[2] eq "mRNA"){
						if($info{ID}!~/^$info{Parent}\.\d+$/){print "$info{Parent} transcript_id $info{ID} is not same ids\n";$gene_name_forcal++;}
						$gff{$line[0]}{$info{ID}}{gene}=$info{Parent};
						$gff{$line[0]}{$info{ID}}{strand}=$line[6];
						if(defined $repr{$info{ID}}){$gff{$line[0]}{$info{ID}}{repr}="yes";}
						$gff{$line[0]}{$info{Parent}}{'5utr'}="";
						$gff{$line[0]}{$info{Parent}}{'5utr_inf'}="";
						$gff{$line[0]}{$info{Parent}}{'cds'}="";
						$gff{$line[0]}{$info{Parent}}{'3utr'}="";
						$gff{$line[0]}{$info{Parent}}{'3utr_inf'}="";
				}elsif($line[2] eq "five_prime_UTR"){
						my $length=$line[4] - $line[3] +1;
						$gff{$line[0]}{$info{Parent}}{'5utr'}.="$line[3]-$line[4];";
						$gff{$line[0]}{$info{Parent}}{'5utr_inf'}.="$length,";
				}elsif($line[2] eq "CDS"){
						my $length=$line[4] - $line[3] +1;
						$gff{$line[0]}{$info{Parent}}{'cds'}.="$length,";
				}elsif($line[2] eq "three_prime_UTR"){
						my $length=$line[4] - $line[3] +1;
						$gff{$line[0]}{$info{Parent}}{'3utr'}.="$line[3]-$line[4];";
						$gff{$line[0]}{$info{Parent}}{'3utr_inf'}.="$length,";
				}else{die "ERROR::what gff line?\n$_\n";}
		}
		close GFF;
		return(%gff);
}
#############################################################################################
sub fasta2hash ( $ )
 {
  my ($file,$key,$value,$comment);
  my (%fasta_hash);
  $file=$_[0];
  if ($file =~ /\.zip$/) {
   open (IN,"unzip -p $file |") || die "problem with $file\n";
  } else {
   open (IN,$file) || die "problem with $file\n";
  }
  while (<IN>)
   {
    chomp;
    if (/^>(\S+)/)
     {
	 #my @rec = split(m/\|/,$1);
	# while(scalar(@rec)>0)
	# {
      	#$key = pop(@rec);
	#	last if(length($key)>3);
	# }
      $key = $1;
      $comment = $2;
      if ($key =~ /^.+\|.+\|.+\|(.+)\|$/) {
	$comment = "$key $comment";
      	$key = $1;
      } 
      $fasta_hash{$key}{"comment"}=$comment;
		#if(defined($fasta_hash{$key}{"sequence"}) || $fasta_hash{$key}{"sequence"} ne '')
		#{
		#	print STDERR "sequence $key already exists\n";
			$fasta_hash{$key}{"sequence"}='';
		#}
     } #if (/^>(\w)$/)
    else 
     {
      $key || die "File $file is not a fasta file!\n$key\n$_\n";
      s/\s+//g;
      $fasta_hash{$key}{"sequence"}.=$_;
     } #else 
   } #while (<IN>)
  close IN;
  return (%fasta_hash); 
 } #fasta2hash ( $ )



