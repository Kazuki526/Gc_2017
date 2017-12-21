#!/usr/bin/perl
use strict;
use warnings;

##useage## perl check_all_aligned_cds_descriminate.pl ALL_ALIGNED.TSV ########

my $all_aligned =$ARGV[0];
-e $all_aligned or die "ERROR::not exist input $all_aligned\n";
my %afa = &fasta2hash_by_file("singlecopy_CS_cds/mostlysinglecopyA.fa");
my %bfa = &fasta2hash_by_file("singlecopy_CS_cds/mostlysinglecopyB.fa");
my %dfa = &fasta2hash_by_file("singlecopy_CS_cds/mostlysinglecopyD.fa");
open(IN,"$all_aligned");
my $header = <IN>;
my %header_posi = ();
my @header = split(/\t/,$header);
for(my $i=0;@header>$i;$i++){
		if($header[$i] eq "aid"){$header_posi{aid}=$i;}
		if($header[$i] eq "bid"){$header_posi{bid}=$i;}
		if($header[$i] eq "did"){$header_posi{did}=$i;}
}
my ($outfa,$outnum,$outnucle)=("temporary.fa",0,0);
open(OUT,">ABD_identical_distinguishable_seq.tsv");
print OUT "aid\tbid\tdid\ta_length\tb_length\td_length\taligned_length\tsame_region\tdistinguishable_length\n";
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		my ($aid,$bid,$did)=($line[$header_posi{aid}],$line[$header_posi{bid}],$line[$header_posi{did}]);
		open(TEMP,">$outfa");
		print TEMP ">$aid\n$afa{$aid}{seq}\n>$bid\n$bfa{$bid}{seq}\n>$did\n$dfa{$did}{seq}\n";
		close TEMP;
		my $multiple_aligned = `mafft $outfa 2>/dev/null`;
		my %multiple_aligned = &fasta2hash_by_carstring($multiple_aligned);
		my ($afa,$bfa,$dfa) = ($multiple_aligned{$aid}{seq},$multiple_aligned{$bid}{seq},$multiple_aligned{$did}{seq});
		my $focal = 0;
		my @focal_start=();
		for(my $i=0;length($afa) - 100 >$i;$i++){
				if((substr($afa,$i,100) ne substr($bfa,$i,100))&&(substr($bfa,$i,100) ne substr($dfa,$i,100))&&(substr($dfa,$i,100) ne substr($afa,$i,100))){
						$focal++;push(@focal_start,$i);
				}
		}
		my($alen,$blen,$dlen,$alignlen)=(length($afa{$aid}{seq}),length($bfa{$bid}{seq}),length($dfa{$did}{seq}),length($afa));
		if($focal==0){
				print OUT "$aid\t$bid\t$did\t$alen\t$blen\t$dlen\t$alignlen\t0\t0\n";
		}else{
				my @focal=();
				my ($start,$end,$focal_length) = (0,0,0);
				foreach my $erst (@focal_start){
						if(($end != 0) && ($end +1 < $erst )){
								$start++;$end++;push(@focal,"$start-$end");$focal_length += ($end - $start +1);
								$start=$erst;$end=$erst+100;
						}elsif($end == 0){
								$start = $erst;$end = $erst +100;
						}else{$end = $erst +100;
						}
				}
				$start++;$end++;push(@focal,"$start-$end");$focal_length += $end - $start +1;
				print OUT "$aid\t$bid\t$did\t$alen\t$blen\t$dlen\t$alignlen\t" . join(";",@focal) . "\t" . $focal_length . "\n";
		}
}
close IN;
close OUT;
unlink $outfa;


###############################################################################
sub fasta2hash_by_carstring ( $ )
 {
  my ($key,$value,$comment);
  my (%fasta_hash);
  my $fasta=$_[0];
  foreach(split(/\n/,$fasta)){
    if (/^>(\S+)\s*(\S*)/){
      $key = $1;
      $comment = $2;
      if ($key =~ /^.+\|.+\|.+\|(.+)\|$/) {
	$comment = "$key $comment";
      	$key = $1;
      } 
      $fasta_hash{$key}{"comment"}=$comment;
      $fasta_hash{$key}{"seq"}='';
     } #if (/^>(\w)$/)
    else 
     {
      $key || die "File $fasta is not a fasta file!\n$key\n$_\n";
      s/\s+//g;
      $fasta_hash{$key}{"seq"}.=$_;
     }
   }
  return (%fasta_hash); 
 } #fasta2hash ( $ )
###############################################################################
sub fasta2hash_by_file ( $ )
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
		#if(defined($fasta_hash{$key}{"seq"}) || $fasta_hash{$key}{"seq"} ne '')
		#{
		#	print STDERR "seq $key already exists\n";
			$fasta_hash{$key}{"seq"}='';
		#}
     } #if (/^>(\w)$/)
    else 
     {
      $key || die "File $file is not a fasta file!\n$key\n$_\n";
      s/\s+//g;
      $fasta_hash{$key}{"seq"}.=$_;
     } #else 
   } #while (<IN>)
  close IN;
  return (%fasta_hash); 
 } #fasta2hash ( $ )
