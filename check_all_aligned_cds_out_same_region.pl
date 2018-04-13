#!/usr/bin/perl
use strict;
use warnings;

##useage## perl check_all_aligned_cds_out_same_region.pl ALL_ALIGNED.TSV ########

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
open(OUT,">ABD_identical_remove_same_seq.tsv");
open(OUTPOSI,">ABD_identical_aligned_position.tsv");
print OUT "aid\tbid\tdid\ta_length\tb_length\td_length\taligned_length\tsame_region\tdistinguishable_length\n";
print OUTPOSI "chr\talined_posi\tposi\tA_genome\tB_genome\tD_genome\n";
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
		my @error_start=();
		my ($aposi,$bposi,$dposi) = (0,0,0);
		for(my $alined=1;length($afa) - 1 >$alined;$alined++){
				my($aposi_out,$bposi_out,$dposi_out);
				if(substr($afa,$alined,1) ne "-"){$aposi++;$aposi_out=$aposi;}else{$aposi_out="NA";}
				if(substr($bfa,$alined,1) ne "-"){$bposi++;$bposi_out=$bposi;}else{$bposi_out="NA";}
				if(substr($dfa,$alined,1) ne "-"){$dposi++;$dposi_out=$dposi;}else{$dposi_out="NA";}
				if((substr($afa,$alined-1,100) eq substr($bfa,$alined-1,100))||
				   (substr($bfa,$alined-1,100) eq substr($dfa,$alined-1,100))||
				   (substr($dfa,$alined-1,100) eq substr($afa,$alined-1,100))){
						$focal++;push(@error_start,$alined);
				}elsif(($alined < 100)&&((substr($afa,0,$alined) eq substr($bfa,0,$alined))||
										 (substr($bfa,0,$alined) eq substr($dfa,0,$alined))||
										 (substr($dfa,0,$alined) eq substr($afa,0,$alined)))){
						$focal++;push(@error_start,$alined);
				}elsif($error_start[$#error_start] + 100  < $alined){
						my ($a_,$b_,$d_) = (substr($afa,$alined-1,1),substr($bfa,$alined-1,1),substr($dfa,$alined-1,1));
						print OUTPOSI "$aid\t$alined\t$aposi_out\t$a_\t$b_\t$d_\n$bid\t$alined\t$bposi_out\t$a_\t$b_\t$d_\n$did\t$alined\t$dposi_out\t$a_\t$b_\t$d_\n";
				}else{$focal++;}
		}
		my($alen,$blen,$dlen,$alignlen)=(length($afa{$aid}{seq}),length($bfa{$bid}{seq}),length($dfa{$did}{seq}),length($afa));
		if($focal==0){
				print OUT "$aid\t$bid\t$did\t$alen\t$blen\t$dlen\t$alignlen\tNA\t$alignlen\n";
		}else{
				my @error=();
				my ($start,$end,$error_length) = (0,0,0);
				foreach my $erst (@error_start){
						if(($end != 0) && ($end +1 < $erst )){
								$start++;$end++;push(@error,"$start-$end");$error_length += ($end - $start +1);
								$start=$erst;$end=$erst+100;
						}elsif($end == 0){
								$start = $erst;$end = $erst +100;
						}else{if($erst +100 <$alignlen){$end = $erst +100;}else{$end = $alignlen;}
						}
				}
				push(@error,"$start-$end");$error_length += $end - $start +1;
				print OUT "$aid\t$bid\t$did\t$alen\t$blen\t$dlen\t$alignlen\t" . join(";",@error) . "\t" . ($alignlen - $error_length) . "\n";
		}
}
close IN;
close OUT;
close OUTPOSI;
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
