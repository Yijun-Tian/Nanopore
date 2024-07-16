#!/usr/bin/perl -w
# This script can substring the MM/ML tags for CpG sites to generate a probability string from a modificaition BAM file by matching the sequence string.
# The probability string can further help to decide the modification model quality from the groundtruth libraries.
# The 2nd argument is the string for searching the background sequences. Here the case is:
# AACAGCACAGCAACTCCCTCGCNNNNNNNNNNNNNGCAGCGACCTGGTGAGCCCG will be AACAGCACAGCAACTCCCTCGC[AGCT]{13}GCAGCGACCTGGTGAGCCCG
# ATCAGATCAGTAGNNNNCGNNNNAGACAGNNNN5hmCGNNNNTAGACTGNNNN5mCGNNNNGCATCGAGTACAGTAC will be ATCAGATCAGTAG[AGCT]{4}CG[AGCT]{4}AGACAG[AGCT]{4}CG[AGCT]{4}TAGACTG[AGCT]{4}CG[AGCT]{4}GCATCGAGTACAGTAC
# remember to remove EcoRI sites
# The 3rd argument is the base modification motif: CG

use strict;
use List::MoreUtils qw(firstidx);
no warnings 'uninitialized';

my $samtools="/app/eb/software/SAMtools/1.16.1-GCC-11.3.0/bin/samtools";
sub main{
	my $bamfile = $ARGV[0];
	my $seq = $ARGV[1];
	# $os is the offset measuring the pattern length
	my $os = $ARGV[2];
	my $outname = $ARGV[3];
	my $motif = $ARGV[4];
	open(OUTFILE, ">$outname.txt") || die("error writing to file $outname.txt\n");
	open(INFILE, "$samtools view $bamfile | " ) || die("Error running samtools to convert BAM file\n");
	while(my $line = <INFILE>){
    		chomp($line);
		my %alignment = getAlignmentInfo($line);
		my $uuid = $alignment{"name"};
		my $SEQ = $alignment{"seq"};
		my $QUAL = $alignment{"qual"};
		my @match = $SEQ =~ /$seq/g;
		next if (!@match);
		my @pos = match_all_positions($seq, $SEQ, $os);
		#print join("|", @pos),"\n";
		next if (!$alignment{"mm"});
		my @mm = split /\;/, $alignment{"mm"};
		my $hp = $mm[0] =~ s/MM:Z:C\+h(\W+)?,//r;
		#print $hp,"\n";
		my $ml = $alignment{"ml"} =~ s/ML:B:C,//r ;
		my @h = split /\,/, $ml;
		my @m = splice @h, @h/2;
		#print join(",",@m),"|",join(",",@h),"\n"; 
	# based on getA(), the routine getCode() return the 0-index of A bases in the sequences
		my @CGref = getCode($SEQ, $motif);
		for (0 .. $#match) {				
				my $phred = substr($QUAL, $pos[$_], $os);
				my $interest = scalar(getCode($match[$_], $motif));
				#print $interest,"\n";
				my $bgnvalue = $pos[$_];
				#print $bgnvalue,"\n";
				#my $bgnidx = (sort {$a <=> $b} grep {$_ > $bgnvalue} @CGref)[0];
				#print $bgnidx,"\n";			
				my $bgnidx = firstidx { $_ > $bgnvalue } @CGref;
				my $endidx = $bgnidx + $interest - 1;
				#print $bgnidx," ",$endidx,"\n";
				my @Hmod = @h[$bgnidx..$endidx];
				my @Mmod = @m[$bgnidx..$endidx];
				
				my $i = 0;
				my $j = 0;
				my @res;
				#print length($match[$_]),"\n";		
				for my $i (0 .. length($match[$_])-1){
					if(substr($match[$_], $i, length($motif)) =~ /$motif/) {
						push @res, join("|", $Hmod[$j],$Mmod[$j]);
						$j ++;
						}else{
						push @res, ("#N/A|#N/A");
						}
						}
				#print $uuid,"\t",$_,"\t",$match[$_], "\t",join(",",@res),"\t",$phred,"\n";
				print OUTFILE $uuid,"\t",$_,"\t",$match[$_], "\t", join(",",@res),"\t",$phred,"\n";
		}
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
	$alignment{"mm"} = (grep { /^MM:Z:/ } @tmp) [0];
	$alignment{"ml"} = (grep { /^ML:B:/ } @tmp) [0];
	return %alignment;
}

sub match_all_positions {
    my ($seq, $SEQ, $os) = @_;
    #my $os = shift;
    my @pos;
    while ($SEQ =~ m/$seq/g) {
        push @pos, pos($SEQ) - $os;
    }
    return @pos
}

sub getCode{
	my @ref;
	my ($seq, $char) = @_;
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
