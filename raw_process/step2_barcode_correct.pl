#!/usr/bin/perl
use strict;
use warnings;
use List::MoreUtils qw(uniq);

my($bam,$cell_barode,$output)=@ARGV;

my %barcode;
open BAM,"/public1/Softwares/samtools-1.9/samtools view -@ 8 $bam|" or die $!;
while(<BAM>){
	chomp;
	s/\r//g;
	my @tab = split /\t/;
	my $CR=$tab[18];
	my $CB=$tab[20];
	$CR=~/CR:Z:(\w+)/;
	my $final_CR=$1;
	$CB=~/CB:Z:(\w+)-1/;
	my $final_CB=$1;
	push @{$barcode{$final_CB}},$final_CR;
}

close BAM;

open BARCODE,"$cell_barode" or die $!;
open OUT,">$output" or die $!;
print OUT "CellBarcode\tReadsBarcode\n";
while(<BARCODE>){
	chomp;
	s/\r//g;
	my @tab = split /\,/;
    my @info = split /\-/;
	my $barcode_id = $info[0];
	if(exists $barcode{$barcode_id}){
		my $CRs=join(";",uniq(@{$barcode{$barcode_id}}));
		print OUT "$barcode_id\t$CRs\n";
	}
}
close BARCODE;
close OUT;