#!/usr/bin/env perl 

## This program removes barcode sequences from fastq reads and places
## an identifier in the info line.  This now includes a function for
## correcting barcodes that are off by 1.

## zg 4ix11 - This is a new version that will work with 8, 9, or 10 bp barcoes.
## cab 8nov12 -- lines are read in groups of 4, skipping the regular expression to track where we are
## cab 10nov12 -- use Levenshtein distance to find barcode correction (faster)
## cab 16dec12 -- catch instances where raw data are not in a set of four lines
## cab 29nov18 -- changed to use Text::Levenshtein::XS 'distance';
## cab 1oct19 -- took previous GBS 768 barcode script and modified it to run on paired end reads 
## cab 6may20 -- anchor sequences at primer sequence for amplicon, build hash of corrected mids

use warnings;

## to run on teton
##  module load perl-text-levenshtein-xs
## parse_barcodes_pairedend.pl demux_key.csv forward.fq reverse.fq machinename

use Text::Levenshtein::XS 'distance';

unless (scalar @ARGV > 3){
    die "I need four arguments to run: barcodefile fwdfile revfile MACHINENAME ";
}
my $barcodes = shift(@ARGV);
my $forwardfile = shift (@ARGV);
my $reversefile = shift (@ARGV);
my $devicename = shift (@ARGV);

unless($devicename){
	die "Please provide a machine name that appears after \@ on info lines of fastq file";
}

open (FORWARD, $forwardfile) or die "You are a dummy";
open (REVERSE, $reversefile) or die "You are still a dummy";
open (MIDS, $barcodes) or die "Could not open MID file";

my %mids;
my %midsctr;
my %forwardbarcodeshash;
my %reversebarcodeshash;
my @forwardbarcodearray;
my @reversebarcodearray;

<MIDS>; ## get rid of top line in MIDS file
while(<MIDS>){
    chomp;
    @line = split ',', $_;
    $line[0] =~ tr/[a-z]/[A-Z]/; ## catch lower-case barcode input and make it uppercase
    $line[1] =~ tr/[a-z]/[A-Z]/; ## catch lower-case barcode input and make it uppercase

    $forwardbarcodeshash{$line[0]}++;
    $reversebarcodeshash{$line[1]}++;
    $mids{$line[0]}{$line[1]} = $line[2];
    $midsctr{$line[0]}{$line[1]} = 0; ## initialize counters to zero
}
close (MIDS) or die "Could not close MIDS\n";

@forwardbarcodearray = sort keys %forwardbarcodeshash;
@reversebarcodearray = sort keys %reversebarcodeshash;

## simplify file names by dropping leading folder name
$forwardfile =~ s/.*\/([\w.]+fq)$/$1/; 
$forwardfile =~ s/.*\/([\w.]+fa)$/$1/;
$forwardfile =~ s/.*\/([\w.]+fastq)$/$1/;
$reversefile =~ s/.*\/([\w.]+fq)$/$1/;  
$reversefile =~ s/.*\/([\w.]+fa)$/$1/;
$reversefile =~ s/.*\/([\w.]+fastq)$/$1/;

open (FWDSEQ, "> parsed_"."$forwardfile") or 
    die "Could not open FWDSEQ\n";
open (REVSEQ, "> parsed_"."$reversefile") or 
    die "Could not open REVSEQ\n";
open (FWDCRAP, "> truemiderrors_"."$forwardfile") or die "Could not open CRAP\n";
open (REVCRAP, "> truemiderrors_"."$reversefile") or die "Could not open CRAP\n";

open (FWDOTHER, "> phixOther_"."$forwardfile") or die "Could not open FWDOTHER\n";
open (REVOTHER, "> phixOther_"."$reversefile") or die "Could not open REVOTHER\n";

my $getit = 0;
my $seqcnt = 0;
my $fwdgoodmid = 0;
my $revgoodmid = 0;
my $fwdcorrectedmid = 0;
my $revcorrectedmid = 0;
my $fwdnoprimer = 0;
my $revnoprimer = 0;
my $fwdbarcodelength = 0;
my $revbarcodelength = 0;

my $goodmidpairctr = 0;
my $badmidpairctr = 0;
my $otherpairctr = 0;
my $correctedmidctr = 0;
my $adlen;
my $adrem = 0;
my @ad;
my $bclen = 11 ; # 10 barcode
my $qline;

# foreach (keys %mids){
#     print "Forward: $_ "; 
#     print "--(", scalar keys %{$mids{$_}} , ")--> ";
#     foreach (keys %{$mids{$_}}){
# 	print "$_ ";
#     }
#     print "\n";
# }

while (defined($forward = <FORWARD>) && defined($reverse = <REVERSE>)){
    if( !($forward =~ /^\@$devicename/) || 
	!($reverse =~ /^\@$devicename/) ){  ## input file is foobarred here, with data not in sets of 4 lines
	print "parse error -- $seqcnt\n$_\n";
	while(!($forward =~ /^\@$devicename/)){
	    $forward = <FORWARD>;
	    $reverse = <REVERSE>;
	}
	while(!($reverse =~ /^\@$devicename/)){
	    $forward = <FORWARD>;
	    $reverse = <REVERSE>;
	}
	print "Back on track --> $_\n";
    }

    $seqcnt++;
    if(!($seqcnt % 50000)){
	print "$seqcnt\n";
    }
    chomp($forward = <FORWARD>);
    chomp($reverse = <REVERSE>);

    ($fwdgoodmid, $fwdnoprimer, $fwdcorrectedmid, 
     $fwdbarcodelength, $forwardmid, $forward) = lookupmid($forward, 
							   \%forwardbarcodeshash, 
							   \@forwardbarcodearray);
    ($revgoodmid, $revnoprimer, $revcorrectedmid, 
     $revbarcodelength, $reversemid, $reverse) = lookupmid($reverse, 
							   \%reversebarcodeshash,
							   \@reversebarcodearray);

    $correctedmidctr = $correctedmidctr + $fwdcorrectedmid + $revcorrectedmid;
#     print "$fwdgoodmid $forwardmid ($fwdbarcodelength), $revgoodmid $reversemid ($revbarcodelength)\n";
    if(exists $mids{$forwardmid}{$reversemid} && $fwdgoodmid && $revgoodmid){
	$goodmidpairctr++;
	$midsctr{$forwardmid}{$reversemid}++;

	print FWDSEQ "@"."$mids{$forwardmid}{$reversemid}"." -- "." fwd $seqcnt\n";
	print REVSEQ "@"."$mids{$forwardmid}{$reversemid}"." -- "."rev $seqcnt\n";

	print FWDSEQ "$forward\n";
	print REVSEQ "$reverse\n";
    }
    elsif($revnoprimer && $fwdnoprimer){
	print FWDOTHER "$forward\n";
	print REVOTHER "$reverse\n";
	$otherpairctr++;
    }
    else{
	print FWDCRAP "$forward\n";	
	print REVCRAP "$reverse\n";
	$badmidpairctr++;
    }
    chomp($forward = <FORWARD>);  ### should be info line that separates seq and qual, often +
    chomp($reverse = <REVERSE>);  ### should be info line that separates seq and qual, often +

    if(exists $mids{$forwardmid}{$reversemid} && $fwdgoodmid && $revgoodmid){
	print FWDSEQ "$forward\n";
	print REVSEQ "$reverse\n";
    }
    chomp($forward = <FORWARD>);  ### quality line
    chomp($reverse = <REVERSE>);  ### 

    if(exists $mids{$forwardmid}{$reversemid} && $fwdgoodmid && $revgoodmid){
	makeQline(\$forward, $fwdbarcodelength); ## edit qline in place (pass by reference)
	makeQline(\$reverse, $revbarcodelength);	
	print FWDSEQ "$forward\n";
	print REVSEQ "$reverse\n";	
    }
}
    
open(REPORT, "> parsereport_$forwardfile") or die "Failed to open report file";
print REPORT "Good barcode pair count: $goodmidpairctr (corrected $correctedmidctr individual mids)\n";
print REPORT "Pairs that lacked primer sequence: $otherpairctr\n";
print REPORT "Bad barcode pair count: $badmidpairctr\n";

foreach my $fwd (sort keys %mids){
    foreach $rev (sort keys %{$mids{$fwd}}){
 	print REPORT "$fwd,$rev,$midsctr{$fwd}{$rev},$mids{$fwd}{$rev}\n";
    }
}
close(REPORT) or die "Could not close REPORT\n";
close(FWDSEQ) or die "Could not close FWDSEQ\n";
close(REVSEQ) or die "Could not close REVSEQ\n";
close(FWDCRAP) or die "Could not close FWDCRAP\n";
close(REVCRAP) or die "Could not close FWDCRAP\n";
close(FWDOTHER) or die "Could not close FWDOTHER\n";
close(REVOTHER) or die "Could not close REVOTHER\n";

##--------- sub routines -----------------##

sub correctmid{
    my $testseq = $_[0];
    my $arrayref = $_[1];

    my $corrected;
    my $min = 100;
    my $mindex = -1;
    for my $i (0 .. $#{$arrayref}){
	$dist = distance($testseq, ${$arrayref}[$i]);
	if($min > $dist){
	    $min = $dist;
	    $mindex = $i;
	}
    }
    return(${$arrayref}[$mindex], $min);
}

sub lookupmid{
    my $target = $_[0];
    my $hashref = $_[1];
    my $arrayref = $_[2];
    my $goodmid = 0;
    my $noprimerregion=0;
    my $correctedmid = 0;
    my $mid = '';

    $line10 = substr($target, 0, $bclen); ## substr EXPR,OFFSET,LENGTH
    $line9 = substr($target, 0, ($bclen - 1)); # 9 bp barcode
    $line8 = substr($target, 0, ($bclen - 2)); # 8 bp barcode

    if(${$hashref}{$line10}){ ##  barcode is in $target, remove barcode
	$whichL = 11;
	$goodmid = 1;
	$target =~ s/^$line10//;
	$mid = $line10;
    }	
    elsif(${$hashref}{$line9}){ ##  barcode is in $target, remove barcode 
	$whichL = 10;
	$goodmid = 1;
	$target =~ s/^$line9//;
	$mid = $line9;
    }
    elsif(${$hashref}{$line8}){ ##  barcode is in $target, remove barcode 
	$whichL = 9;
	$goodmid = 1;
	$target =~ s/^$line8//;
	$mid = $line8;
    }
    else { ## potential mid error
	## first look for amplicon primers
	if($target =~ /\w{9,11}TG[CT]CAGC[AC]GCCGCGGTAA/ ||
	#if($target =~ /\w{8,10}GTG[CT]CAGC[AC]GCCGCGGTAA/ ||
	   $target =~ /\w{9,11}GACTAC[ACT][ACG]GGGT[AT]TCTAAT/ ||
	   #$target =~ /\w{8,10}GGACTAC[ACT][ACG]GGGT[AT]TCTAAT/ ||
	   $target =~ /\w{9,11}TTGGTCATTTAGAGGAAGTAA/ ||
	  # $target =~ /\w{8,10}CTTGGTCATTTAGAGGAAGTAA/ ||
	   $target =~ /\w{9,11}CTGCGTTCTTCATCGATGC/){
	   #$target =~ /\w{8,10}GCTGCGTTCTTCATCGATGC/){
	    ($line10b, $n10) = correctmid($line10, $arrayref);
	    $minN = $n10;
	    $whichL = 11;
	    if ($minN > 1){
		($line9b, $n9) = correctmid($line9, $arrayref);
		if ($n9 < $minN){
		    $minN = $n9;
		    $whichL = 10;
		}
		if ($minN > 1){
		    ($line8b, $n8) = correctmid($line8, $arrayref);
		    if ($n8 < $minN){
			$minN = $n8;
			$whichL = 9;
		    }
		}
	    }
	    if ($minN > 1){ ## can't correct
		$goodmid = 0;
	    }
	    else { ## has been corrected
		if ($whichL == 11){
		    $mid = $line10b;
		    $target =~ s/^$line10//;
		}
		elsif ($whichL == 10){
		    $mid = $line9b;
		    $target =~ s/^$line9//;
		}
		elsif ($whichL == 9){
		    $mid = $line8b;
		    $target =~ s/$line8//;
		}
		$goodmid = 1;
		$correctedmid = 1;
	    }
	}
	else{
	    $goodmid=0;
	    $noprimerregion = 1;
	}
    }
    return($goodmid, $noprimerregion, $correctedmid, $whichL, $mid, $target);
}

sub makeQline{
    my $qlineref = $_[0];
    my $barcodelength = $_[1];
    
    if ($barcodelength == 11){
	$seqlength = (length ${$qlineref}) - $bclen ;
	${$qlineref} = substr(${$qlineref}, $bclen, $seqlength);
    }
    elsif ($barcodelength == 10){
	$seqlength = (length ${$qlineref}) - $bclen  + 1;
	${$qlineref} = substr(${$qlineref}, ($bclen - 1), $seqlength);
    }
    elsif ($barcodelength == 9){
	$seqlength = (length ${$qlineref}) - $bclen + 2;
	${$qlineref} = substr(${$qlineref}, ($bclen - 2), $seqlength);
    }
}
