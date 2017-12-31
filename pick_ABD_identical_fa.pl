#/usr/bin/perl
use strict;
use warnings;

######### usage perl pick_ABD_identical_fa.pl ABD_identical___.tsv ##########
#### doing on GC_2017 dir

my %afa = &fasta2hash("singlecopy_CS_cds/mostlysinglecopyA.fa");
my %bfa = &fasta2hash("singlecopy_CS_cds/mostlysinglecopyB.fa");
my %dfa = &fasta2hash("singlecopy_CS_cds/mostlysinglecopyD.fa");

open(IN,"$ARGV[0]");
open(OUT,">ABD_identical.fa");
<IN>;
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		if($line[8] == 0){next;}
		print OUT ">$line[0]\n$afa{$line[0]}{sequence}\n";
		print OUT ">$line[1]\n$bfa{$line[1]}{sequence}\n";
		print OUT ">$line[2]\n$dfa{$line[2]}{sequence}\n";
}
close IN;
close OUT;
exit;

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
