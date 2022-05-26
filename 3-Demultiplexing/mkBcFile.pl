#!/usr/bin/perl
# combines id files for parse perl script
#

open(IDF, "Barcodes/BarcodesForZach/FullNameWithBase/Barcodes_Forward_FullNameWithBase.fasta") or die "failed to open forward file\n";
open(IDR, "Barcodes/BarcodesForZach/FullNameWithBase/Barcodes_Reverse_FullNameWithBase.fasta") or die "failed to open reverse file\n";
open(OUT, "> idsComb.csv") or die "failed to write\n";

while(<IDF>){
	$idf = $_;
	$idr = <IDR>;
	chomp($idf);
	chomp($idr);
	if($idf == $idr){ ## all is good
		$for = <IDF>;
		$rev = <IDR>;
		chomp($for);
		chomp($rev);
		$idf =~ s/^>//;
		print OUT "$for,$rev,$idf\n";
	}
	else{
		print "WTF: $idf\n$idr\n";
	}
}
close(IDF);
close(IDR);
close(OUT);
