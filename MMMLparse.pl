#!/usr/bin/perl -w
# Given MM and ML tags, output the sequence level index for mutation

use strict;
no warnings 'uninitialized';

my $samtools="/home/user/anaconda3/envs/NewSam/bin/samtools";
sub main{
  my $bamfile = $ARGV[0];
  my $outname = $ARGV[1];
  my $prob = $ARGV[2];
  my $start = time;
  open(SAM, ">$outname.sam") || die("error writing to file $outname.sam\n");
  open(HEAD, "$samtools view -H $bamfile |") || die("can't read header\n");
  while(my $line = <HEAD>){
  chomp($line);
  print SAM $line,"\n";
  }
  open(INFILE, "$samtools view -d MM $bamfile | " ) || die("Error running samtools to convert BAM file\n");
  while(my $line = <INFILE>){
    chomp($line);
    my %alignment = getAlignmentInfo($line);
    my $uuid = $alignment{"name"};
    my $SEQ = $alignment{"seq"};
    #print join(",",@ref), "\n";
    next if (!$alignment{"mm"});
    my $mm = $alignment{"mm"} =~ s/MM:Z:C\+m(\W+)?,//r ;
    #/MM:Z:C\+m(\W+)?,
    #MM:Z:C matches the characters MM:Z:C literally (case sensitive)
    #\+ matches the character + with index 4310 (2B16 or 538) literally (case sensitive)
    #m matches the character m with index 10910 (6D16 or 1558) literally (case sensitive)
    #1st Capturing Group (\W+)?
    #? matches the previous token between zero and one times, as many times as possible, giving back as needed (greedy)
    #\W matches any non-word character (equivalent to [^a-zA-Z0-9_])
    #+ matches the previous token between one and unlimited times, as many times as possible, giving back as needed (greedy)
    #, matches the character , with index 4410 (2C16 or 548) literally (case sensitive)
    #Global pattern flags 
    #g modifier: global. All matches (don't return after first match)
    my $ml = $alignment{"ml"} =~ s/ML:B:C,//r ;
    my $strand = $alignment{"flag"} & 16 ? "C":"W";
    #print "$strand","\n";
    #print "$mm", "\t", "$ml", "\n";
    my $seq;
    my $conv;
    my @ref;
    #print "$SEQ","\n";
    if($strand eq "C"){
    $seq = RC($SEQ);
    }
    else{
    $seq = $SEQ;
    }
    @ref = getC($seq);
    #print "$seq","\n";
    my $i = 0;
    my $j = 0;
    my @int;
    my @score = split /,/, $ml;
    my @convcand;
    for (0 .. $#score) {
       if ($score[$_] <= $prob){
           push @convcand, $_;
       }
       }
    while($i < length($mm)){
             substr($mm, $i) =~ /^([0-9]+)/; 
             defined($1) || die("Could not parse cytosine at pos $i: '$mm'");
		my $runlen = $1;
		$i += length($1);
		$i ++;
		$j += $runlen + 1;		
		push @int, $j-1;
    }
    my @mods =  @int[@convcand];
    my @meth = @ref[@mods];
    substr( $seq, $_, 1 ) =~ tr[C][T] for @meth;
    if($strand eq "C"){
    $conv = RC($seq);
    }
    else{
    $conv = $seq;
    }
    #print "$conv","\n";
    if($alignment{"ps"}){
    print SAM $alignment{"name"},"\t",
              $alignment{"flag"},"\t",
              $alignment{"chr"},"\t",
              $alignment{"pos"},"\t",
              $alignment{"mapq"},"\t",
              $alignment{"cigar"},"\t*\t0\t0\t",
              $conv,"\t",
              $alignment{"qual"},"\t",
              $alignment{"mm"},"\t",
              $alignment{"ml"},"\t",
              $alignment{"ps"},"\t",
              $alignment{"hp"},"\n";
     }else{
     print SAM $alignment{"name"},"\t",
              $alignment{"flag"},"\t",
              $alignment{"chr"},"\t",
              $alignment{"pos"},"\t",
              $alignment{"mapq"},"\t",
              $alignment{"cigar"},"\t*\t0\t0\t",
              $conv,"\t",
              $alignment{"qual"},"\t",
              $alignment{"mm"},"\t",
              $alignment{"ml"},"\n";}
}
close(SAM);

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
	$alignment{"mm"} = (grep { /^MM:Z:C/ } @tmp) [0];
	$alignment{"ml"} = (grep { /^ML:B:C/ } @tmp) [0];
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
        
sub RC{
      my $seq = shift;
      my $revcomp = reverse $seq;
      $revcomp =~ tr/ATGC/TACG/;
      return $revcomp;
      }        
main();
