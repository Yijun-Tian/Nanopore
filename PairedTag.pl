#!/usr/bin/perl -w
# Extract paired tag from MMML tags, such as 5mC and 5hmC at the same time

use strict;
no warnings 'uninitialized';
my $deletion_string = "D";
#my $samtools="/usr/bin/samtools";
my $samtools="/app/eb/software/SAMtools/1.16.1-GCC-11.3.0/bin/samtools";
sub main{
  my $bamfile = $ARGV[0];
  my $outname = $ARGV[1];
  #my $region = $ARGV[2];
  my $start = time;
  open(OUTFILE, ">$outname.txt") || die("error writing to file $outname.txt\n");
  #open(INFILE, "$samtools view --region-file $region $bamfile | " ) || die("Error running samtools to convert BAM file\n");
  open(INFILE, "$samtools view $bamfile | grep 676ebccf-dcdc-4dd3-bd36-64131e8c4e12 | " ) || die("Error running samtools to convert BAM file\n");
  while(my $line = <INFILE>){
    chomp($line);
    my %alignment = getAlignmentInfo($line);
    #next if (!$alignment{"ml"});
    getMap(\%alignment);
	print $alignment{"seq"},"\n";
	print join(",", getC($alignment{"seq"})),"\n";
	print $alignment{"pseudo_seq"},"\n";
	print $alignment{"ref_match_seq"},"\n";
    my $mm = (split /;/, $alignment{"mm"}) [0] =~ s/MM:Z:C\+h,//r ;
    print $mm,"\n";
    my $ml = $alignment{"ml"} =~ s/ML:B:C,//r ;
    my @mlbc = split(',', $ml);
    my @mlbcA = splice @mlbc, @mlbc/2;
    #next if scalar(@mlbc) ne scalar(@mlbcA);
	#for (0 .. scalar(@mlbc)) {
	#	print OUTFILE $mlbc[$_],"\t",$mlbcA[$_],"\n";
	#}
}
close(OUTFILE);    
}
sub getAlignmentInfo{
	my $line = shift;
	my @tmp = split /\t/, $line;
	my %alignment;
	$alignment{"name"} = $tmp[0];
	$alignment{"flag"} = $tmp[1];
	$alignment{"chr"} = $tmp[2];
	$alignment{"pos"} = $tmp[3];
	$alignment{"mapq"} = $tmp[4];
	$alignment{"cigar"} = $tmp[5];
	#*
	#0
	#0
	$alignment{"seq"} = $tmp[9];
	$alignment{"qual"} = $tmp[10];
	$alignment{"ps"} = (grep { /^PS:i:/ } @tmp) [0];
	$alignment{"hp"} = (grep { /^HP:i:/ } @tmp) [0];
	$alignment{"mm"} = (grep { /^MM:Z:C+/ } @tmp) [0];
	$alignment{"ml"} = (grep { /^ML:B:C,/ } @tmp) [0];
	return %alignment;
}

#The getMap routine will convert the seq fields into a mapped format, with deletion filled with "D", and insertion/softclipping changed into lower case
sub getMap{
	my $alignment = shift;
	my $cigar = $alignment->{"cigar"};
	my $query = $alignment->{"seq"};
	my $qual = $alignment->{"qual"};
	my $pse_seq;
	my $ref_seq;
	my $ref_qual;
	my $i = 0;
	my $j = 0;
	my $nm_i = 0;
	my $nm_d = 0;
	while($i < length($cigar)){
		substr($cigar, $i) =~ /^([0-9]+)/;
		defined($1) || die("Could not parse number at pos $i: '$cigar'");
		my $runlen = $1;
		$i += length($1);
		$i < length($cigar) || die("Bad cigar string : '$cigar'");
		my $op = substr($cigar, $i, 1);
		defined($op) || die("count not parse operation at pos $i: '$cigar'");
		die("Could not understand $op: '$cigar'") if($op !~ m/[MX=DIS]/);
		$op =~ s/[X=]/M/g;
		my ($clip_s, $clip_q);
		if($op eq "M" || $op eq "I" || $op eq "S"){
			$clip_s = substr($query, $j, $runlen);
			$clip_q = substr($qual, $j, $runlen);
			$clip_s = lc($clip_s) if($op eq "I" || $op eq "S"); 
			$nm_i += $runlen if($op eq "I");
			$j += $runlen;
		}else{
			#deletion from reference
			$nm_d += $runlen;
			#length($deletion_string) > $runlen || die("deletion is too long at $runlen: '$cigar'");
			$clip_s = $deletion_string x $runlen;
			$clip_q = $deletion_string x $runlen;
		}
		$i++;
		$pse_seq = $ref_seq . $clip_s if($op =~ m/[MDIS]/);
		$ref_seq = $ref_seq . $clip_s if($op =~ m/[MD]/);
		$ref_qual = $ref_qual . $clip_q if($op =~ m/[MD]/);
	}
	$alignment->{"pseudo_seq"} = $pse_seq;
	$alignment->{"ref_match_seq"} = $ref_seq;
	$alignment->{"ref_match_qual"} = $ref_qual;
	$alignment->{"ref_nm_i"} = $nm_i;
	$alignment->{"ref_nm_d"} = $nm_d;
	$alignment->{"tlen"} = length($ref_seq) - $nm_i if($ref_seq);
}




sub getPos{
	my @ref;
        my $seq = shift;
        my $char = "C";
        my $offset = 0;
        my $result = index($seq, $char, $offset);
        while ($result != -1){
           push(@ref, $result);
           $offset = $result+1;
           $result = index($seq, $char, $offset);
        }
        return @ref;
        }


main();
