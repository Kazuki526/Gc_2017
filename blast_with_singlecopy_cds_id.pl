#!/usr/bin/perl
use strict;
use warnings;

###usage### perl blast_with_singlecopy_cds_id.pl "perfectly" or "mostly" ####
### ^=doing this perl on project dir ###

my $filter = $ARGV[0];
if(($filter ne "perfectly")&&($filter ne "mostly")){die "ERROR::set filter to \"perfectly\" or \"mostly\"\n";}

my $sing_dir = "singlecopy_CS_cds";
-e $sing_dir or die "ERROR::there is not exist $sing_dir\ndoing on wrong directri??\n";

my $a_fasta = &make_fasta("a");
my $b_fasta = &make_fasta("b");
my $d_fasta = &make_fasta("d");

mkdir "$sing_dir/DB";
`makeblastdb -in $a_fasta -out $sing_dir/DB/singlecopyA -dbtype nucl -parse_seqids`;
`makeblastdb -in $b_fasta -out $sing_dir/DB/singlecopyB -dbtype nucl -parse_seqids`;
`makeblastdb -in $d_fasta -out $sing_dir/DB/singlecopyD -dbtype nucl -parse_seqids`;

my $blast_out = "$sing_dir/blast_A2B_$filter.tsv";
print "blast between A & B print out $blast_out\n";
`blastn -db $sing_dir/DB/singlecopyB -query $a_fasta -out $blast_out -outfmt 6 -evalue 1e-6`;
 $blast_out = "$sing_dir/blast_A2D_$filter.tsv";
print "blast between A & D print out $blast_out\n";
`blastn -db $sing_dir/DB/singlecopyD -query $a_fasta -out $blast_out -outfmt 6 -evalue 1e-6`;

###############################################################################
sub make_fasta( $ ){
		my $genome = uc $_[0];
		my $cdsid_file = "$sing_dir/$filter" . "_singlecopy_" . $genome . ".tsv";
		-e $cdsid_file or die "ERROR:: not exist $cdsid_file\n";
		my $HLmix_fasta_file = "CS_cds/DB/HLmix/CS_HLC_repr_cds" . $genome . ".fa";
		-e $HLmix_fasta_file or die "ERROR:: not exist $HLmix_fasta_file\n";
		my %hlc_fa = &fasta2hash($HLmix_fasta_file);
		open(IN,"$cdsid_file") or die "ERROR::problem with $cdsid_file\n";
		<IN>; #remove header line
		my $out_fasta = "$sing_dir/$filter" . "singlecopy" . $genome . ".fa"
		open(OUT,">$out_fasta");
		while(<IN>){
				chomp;
				my @line = split(/\t/,);
				print OUT ">$line[0]\n$hlc_fa{$line[0]}{sequence}\n";
		}
		close IN;
		close OUT;
		return($out_fasta);
}


###############################################################################
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
