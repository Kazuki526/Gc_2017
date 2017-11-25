#!/usr/bin/perl
use warnings;
use strict;

my %reprfa=&fasta2hash($ARGV[0]);
my %allfa =&fasta2hash($ARGV[1]);
open(OUTA,"|gzip -c >$ARGV[2]"."_repr_cdsA.fa.gz");
open(OUTB,"|gzip -c >$ARGV[2]"."_repr_cdsB.fa.gz");
open(OUTD,"|gzip -c >$ARGV[2]"."_repr_cdsD.fa.gz");
open(OUTAO,"|gzip -c >$ARGV[2]"."_notrepr_cdsA.fa.gz");
open(OUTBO,"|gzip -c >$ARGV[2]"."_notrepr_cdsB.fa.gz");
open(OUTDO,"|gzip -c >$ARGV[2]"."_notrepr_cdsD.fa.gz");
foreach my $tr_id(sort(keys %reprfa)){
		if($tr_id =~ /^TraesCS\d([ABD])/){
				my $genome=$1;
				if($genome eq "A"){
						print OUTA ">$tr_id\n$reprfa{$tr_id}{sequence}\n";
				}elsif($genome eq "B"){
						print OUTB ">$tr_id\n$reprfa{$tr_id}{sequence}\n";
				}elsif($genome eq "D"){
						print OUTD ">$tr_id\n$reprfa{$tr_id}{sequence}\n";
				}
		}elsif($tr_id =~ /^TraesCSU/){next;
		}else{die "ERROR::what transcript id? $tr_id\n";}
}
foreach my $tr_id(sort(keys %allfa)){
		if(defined $reprfa{$tr_id}){next;}
		if($tr_id =~ /^TraesCS\d([ABD])/){
				my $genome=$1;
				if($genome eq "A"){
						print OUTAO ">$tr_id\n$allfa{$tr_id}{sequence}\n";
				}elsif($genome eq "B"){
						print OUTBO ">$tr_id\n$allfa{$tr_id}{sequence}\n";
				}elsif($genome eq "D"){
						print OUTDO ">$tr_id\n$allfa{$tr_id}{sequence}\n";
				}
		}elsif($tr_id =~ /^TraesCSU/){next;
		}else{die "ERROR::what transcript id? $tr_id\n";}
}
close OUTA;
close OUTB;
close OUTD;
close OUTAO;
close OUTBO;
close OUTDO;
		





##########################################################################3
sub fasta2hash ( $ )
 {
  my ($file,$key,$value,$comment);
  my (%fasta_hash);
  $file=$_[0];
  if ($file =~ /\.zip$/) {
   open (IN,"unzip -p $file |nkf -wd|") || die "problem with $file\n";
  } else {
   open (IN,$file) || die "problem with $file\n";
  }
  while (<IN>)
   {
    chomp;
    if (/^>(\S+)\s*(\S*)/)
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
