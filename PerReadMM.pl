#!/usr/bin/perl -w
# ./PerReadMM.pl BAM/acont.sorted.bam 243,243 acont.out.th243
# ./PerReadMM.pl R10P6.Haptag.mod.bam 225,225 out
#./PerReadMM.pl BAM/acont.sorted.bam 204,204 chr22 chr22.acont.out.th204
#./PerReadMM.pl BAM/acont.sorted.bam 204,204 chrM chrM.acont.out.th204

my $samtools="/app/eb/software/SAMtools/1.16.1-GCC-11.3.0/bin/samtools";
sub main{
  my $bamfile = $ARGV[0];
  my @th = split(",", $ARGV[1]);
  my $th5mC = $th[0];
  my $th5mH = $th[1];
  my $chromosome = $ARGV[2];
  my $outname = $ARGV[3];
  open(OUTFILE, ">$outname.txt") || die("error writing to file $outname.txt\n");
  print OUTFILE "readname", "\t", "strand", "\t", "perReadQscore", "\t", "readlength", "\t", "n_C", "\t", "n_call", "\t", "n_5mC", "\t", "n_5hmC", "\n"; 
  open(INFILE, "$samtools view -d MM -F 2304 $bamfile $chromosome | " ) || die("Error running samtools to convert BAM file\n");
  while(my $line = <INFILE>){
    chomp($line);
    my %alignment = getAlignmentInfo($line);
    next if (!$alignment{"ml"});
    my $uuid = $alignment{"name"};
    my @ref = getC($alignment{"seq"});
    my $strand = $alignment{"flag"} & 16 ? "C" : "W";
    my $qs = $alignment{"qs"} =~ s/qs:i://r ;
    my $ml = $alignment{"ml"} =~ s/ML:B:C,//r ;
    my @mlbc = split(',', $ml);
    my @mlbcA = splice @mlbc, @mlbc/2;
    next if scalar(@mlbc) ne scalar(@mlbcA);
    my $m = $h = 0;
	for (0 .. $#mlbc) {
		next if ($mlbc[$_] < $th5mC);
		$h ++;
	}
	for (0 .. $#mlbcA) {
		next if ($mlbcA[$_] < $th5mH);
		$m ++;
	}
	print OUTFILE $uuid, "\t" , $strand, "\t", $qs, "\t", length($alignment{"seq"}), "\t", scalar(@ref), "\t", scalar(@mlbc), "\t", $m, "\t", $h, "\n";
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
	$alignment{"qs"} = (grep { /^qs:i:/ } @tmp) [0];
	return %alignment;
}

sub getC{
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
