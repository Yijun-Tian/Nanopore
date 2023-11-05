#!/usr/bin/perl -w
# This script can substring the MM/ML tags to generate a probability string from a modificaition BAM file by matching the sequence string.
# The probability string can further help to decide the modification model quality from the groundtruth libraries.
# The 2nd argument is the string for searching the background sequences. Here the case is:
# AACAGCACAGCAACTCCCTCGCNNNNNNNNNNNNNGCAGCGACCTGGTGAGCCCG will be AACAGCACAGCAACTCCCTCGC[AGCT]{13}GCAGCGACCTGGTGAGCCCG

use strict;
no warnings 'uninitialized';

my $samtools="/usr/local/bin/samtools";
sub main{
	my $bamfile = $ARGV[0];
	my $seq = $ARGV[1];
	my $os = $ARGV[2];
	my $outname = $ARGV[3];
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
		#print join(",",@match)," ", join(",",@pos),"\n";		
		next if (!$alignment{"mm"});
		my $mm = $alignment{"mm"} =~ s/MM:Z:A\+a(\W+)?,//r ;
		#/MM:Z:C\+m(\W+)?,
		#MM:Z:C matches the characters MM:Z:C literally (case sensitive)
		#\+ matches the character + with index 4310 (2B16 or 538) literally (case sensitive)
		#a matches the character m with index 10910 (6D16 or 1558) literally (case sensitive)
		#1st Capturing Group (\W+)?
		#? matches the previous token between zero and one times, as many times as possible, giving back as needed (greedy)
		#\W matches any non-word character (equivalent to [^a-zA-Z0-9_])
		#+ matches the previous token between one and unlimited times, as many times as possible, giving back as needed (greedy)
		#, matches the character , with index 4410 (2C16 or 548) literally (case sensitive)
		#r modifier: regex substitution modifier.
		#print "$mm","\n";
		my $ml = $alignment{"ml"} =~ s/ML:B:C,//r ;
		#print $SEQ,"\n";
		my @Aref = getA($SEQ);
		#print join(",",@Aref), "\n";
		#print $Aref[$0], "\n";
		#print $bgnidx,"\n";
		my @score = split /,/, $ml;
			for (0 .. $#match) {
				
				my $phred = substr($QUAL, $pos[$_], $os);
				my $interest = scalar(getA($match[$_]));
				my $bgnvalue = $pos[$_];
				#print $bgnvalue,"\n";				
				my ($bgnidx) = grep { $Aref[$_] ~~ $bgnvalue } 0 .. $#Aref;
				my $endidx = $bgnidx + $interest - 1;
				#print $bgnidx," ",$endidx,"\n";
				my @mods = @score[$bgnidx..$endidx];
				#print join(",",@mods), "\n";
				my $i = 0;
				my $j = 0;
				my @res;
				#print length($match[$_]),"\n";		
				for my $i (0 .. length($match[$_])-1){
					if(substr($match[$_], $i, 1) =~ /A/) {
						push @res, $mods[$j];
						$j ++;
						}else{
						push @res, "#N/A";
						}
						}
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
	$alignment{"mm"} = (grep { /^MM:Z:A/ } @tmp) [0];
	$alignment{"ml"} = (grep { /^ML:B:C/ } @tmp) [0];
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

sub getA{
	my @ref;
        my $seq = shift;
        my $char = "A";
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
