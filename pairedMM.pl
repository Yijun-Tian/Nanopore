#!/usr/bin/perl -w
# Extract paired tag from pairedMM.py output, such as 5mC and 5hmC probability and mapped locations

use strict;
no warnings 'uninitialized';

sub main{
  my $infile = $ARGV[0];
  my $outname = $ARGV[1];
  my $start = time;
  open(OUTFILE, ">$outname.out") || die("error writing to file $outname.out\n");
  open(INFILE, "head -n 3 $infile | " ) || die("Error open $infile\n");
  while(my $line = <INFILE>){
	chomp($line);
	my %alignment = getFields($line);
	next if ($alignment{"tag"} =~ "None");
		my %h = parseTags($alignment{"h"});
		my %m = parseTags($alignment{"m"});
		print $alignment{"chrm"},"|",$alignment{"pos"},"|", $h{"strand"},"|", $h{"coor"},"|",$h{"prob"},"\n";
		#for 
}
close(OUTFILE);    
}

sub getFields{
	my $line = shift;
	my @tmp = split /\t/, $line;
	my %alignment;
	$alignment{"chrm"} = $tmp[0];
	$alignment{"pos"} = $tmp[1];
	$alignment{"tag"} = $tmp[2] =~ s/[{\'\[}]//rg ;
	my @tags = split /\]/, $alignment{"tag"};
	$alignment{"h"} = $tags[0];
	$alignment{"m"} = $tags[1] =~ s/^,\ //r;
	return %alignment;
}

sub parseTags{
	my $tag = shift;
	my @tmp = split /\:/, $tag;
	my %mod;
	$mod{"strand"} = $tmp[0] =~ s/\(C, ([0-9]+), [mh]\)/$1/rg;
	$mod{"coor"} = $tmp[1] =~ s/^ //r =~ s/\(([(0-9]+),\ [0-9)]+\)/$1/rg;
	$mod{"prob"} = $tmp[1] =~ s/^ //r =~ s/\([0-9]+,\ ([0-9]+)\)/$1/rg;
	return %mod;
}

main();
