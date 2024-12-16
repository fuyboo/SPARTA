#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my($barcode,$R1,$R2,$output);

&GetOptions(
	'barcode=s'   => \$barcode,
	'R1=s'		=>\$R1,
	'R2=s'		=>\$R2,
	'o=s'		=>\$output
);
unless ($barcode || $R1 || $R2) {
	usage();
}

open BAR,$barcode or die $!;
my %cell_index;
my $barcode_length;
while(<BAR>){
	chomp;
    s/\r//g;
    my @tab = split /\t/;
    my @seq_info = split /;/,$tab[1];
    for(my $i=0;$i<@seq_info;$i++){
        $cell_index{$seq_info[$i]}=$tab[0];
    }
    $barcode_length=length($tab[0]);
}
close BAR;

open IN1,"gzip -dc $R1|" or die $!;
open IN2,"gzip -dc $R2|" or die $!;
open OUT,"> $output" or die $!;
while(my $id1=<IN1>){
	my $seq1=<IN1>;
	my $text1=<IN1>;
	my $qual1=<IN1>;
	my $id2=<IN2>;
	my $seq2=<IN2>;
	my $text2=<IN2>;
	my $qual2=<IN2>;
	$id1=~s/\r|\n//g;
	$seq1=~s/\r|\n//g;
	$qual1=~s/\r|\n//g;
	$text1=~s/\r|\n//g;
	$id2=~s/\r|\n//g;
	$seq2=~s/\r|\n//g;
	$text2=~s/\r|\n//g;
	$qual2=~s/\r|\n//g;
	my $cell_barcode=substr($seq1,0,$barcode_length);
    my $umi_reads=substr($seq1,$barcode_length);
	if(exists $cell_index{$cell_barcode}){
        if($seq2=~ /GCAGCTCGGCCC/){
            if(length($seq2)>55){
                $seq2 =~/(\w+)GCAGCTCGGCCC(\w+)/;
                my $adapter=$1;
                print OUT "$cell_index{$cell_barcode}$umi_reads\t$adapter\n";
            }else{
                print OUT "$cell_index{$cell_barcode}$umi_reads\t$seq2\n";
            }
        }else{
            print OUT "$cell_index{$cell_barcode}$umi_reads\t$seq2\n";
        }
	}
}

sub usage {
	die "

Usage:   $0 [Options]
Options:

	--barcode		barcode info,required
	--R1			input R1 reads,required
	--R2			input R2 reads,required
	-o			output required

Example: $0 --barcode barcode.txt --R1 A1.R1.fastq.gz --R2 A1.R2.fastq.gz -o all_barcode_umi_reads.txt

