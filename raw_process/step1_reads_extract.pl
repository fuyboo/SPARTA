#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my($R1,$R2,$out_R1,$out_R2);

&GetOptions(
        'R1=s'          =>\$R1,
        'R2=s'          =>\$R2,
        'o1=s'          =>\$out_R1,
        'o2=s'          =>\$out_R2,
);
unless ($R1 || $R2) {
        usage();
}



open IN1,"gzip -dc $R1|" or die $!;
open IN2,"gzip -dc $R2|" or die $!;
open OUT1,"| gzip > $out_R1" or die $!;
open OUT2,"| gzip > $out_R2" or die $!;
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
		if($seq2 =~ /GCAGCTCGG/){
			$seq2 =~/(\w+)GCAGCTCGG(\w+)/;
			my $adapter=$1;
			my $final_qual2=substr($qual2,0,length($adapter));
			my $final_seq1=substr($seq1,0,26);
			my $final_qual1=substr($qual1,0,26);
			print OUT1 "$id1\n$final_seq1\n$text1\n$final_qual1\n";
			print OUT2 "$id2\n$adapter\n$text2\n$final_qual2\n";
		}
}

close IN1;
close IN2;

sub usage {
	die "

Usage:   $0 [Options]
Options:

	--R1			input R1 reads,required
	--R2			input R2 reads,required
	--o1			output R1 reads,required
	--o2			output R2 reads,required

Example: $0 --R1 A1.R1.fastq.gz --R2 A1.R2.fastq.gz --o1 A1.extra.R1.fastq.gz --o2 A1.extra.R2.fastq.gz


Version: 1.0
Contact:	clge\@lc-bio.com
Last Update:	2023/05/30
";
}