#!/usr/bin/perl

## Creates files for each individual/locus from one or more fastq files. Individuals ids are given in a separate file with one ID per line.

use warnings;

$idfile = shift (@ARGV);
$pref = shift(@ARGV);
open (IDS, $idfile) or die;
while (<IDS>){
    chomp;
    $fh = "FQ"."$_";
    open ($fh, "> $_"."$pref".".fastq") or die "Could not write\n";
    $files{$_} = $fh;
}
close (IDS);

foreach $in (@ARGV){
    open (IN, $in) or die;
    while (<IN>){
	chomp;
	if (/^\@/){
	    if(m/ITS/){
		$id = "ITS";
	    }
	    elsif(m/16S/){
		$id = "16S";
            }
            else{
		print "Warning, no match here $_\n";
	    }	
	    ##s/ \-\- /\-/;
	    if (defined $files{$id}){
		$flag = 1;
		print {$files{$id}} "$_\n";
	    }
	    else{
		$flag = 0;
	    }
	    foreach $i (0..2){
		$line = <IN>;
		if ($flag == 1){
		    print {$files{$id}} "$line";
		}	
	    }	
	}
	else{
	    print "Error -- $_\n";
	}
    }
    close (IN);
}
foreach $id (keys %files){
    close ($files{$id});
}

