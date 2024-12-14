#!/usr/bin/perl
use List::MoreUtils qw(uniq);

my ($input,$output,$output_matrix_dir)=@ARGV;

open IN,$input or die $!;
my %index;
my %cell_num;
my %cell;
my %aptamer_info;
my %aptamer_id;
my %barcodes;
my $id_num=1;
my $row_num=1;
my $all_umi=0;
while(<IN>){
    chomp;
    $all_umi++;
    my @tab = split /\t/;
    my $barcode=substr($tab[0],0,16);
    if(!(exists $barcodes{$barcode})){
        $barcodes{$barcode}=$row_num;
        $row_num++;
    }
    # my $umi=substr($tab[0],16);
    $index{$barcode}{$tab[1]}++;  ###已做过去重，相同的UMI已过滤，直接通过数量即可
    $cell_num{$barcode}++;
    push @{$cell{$barcode}},$tab[1];
    if(!(exists $aptamer_info{$tab[1]})){
        my $add_info="id.$id_num";
        $aptamer_info{$tab[1]}="$tab[1]\t$add_info";
        $aptamer_id{$tab[1]}="$id_num";
        $id_num++;
    }       
}
close IN;

# print "all_num\t$all_umi\n";
# print "barcode_num\t$row_num\n";
# print "aptamer_num\t$id_num\n";
open OUT,">$output" or die $!;
print OUT "Cell_ID\tadpater_Num\tadpater_info\tadpater_umi\ttotal_umi\n";
foreach my $key(sort { $cell_num{$b} <=> $cell_num{$a} } keys %cell_num){
    my @adapter=uniq(@{$cell{$key}});
    my $adapter_num=@adapter;
    my @adapter_value;
    foreach my $element (@adapter){
        push @adapter_value,$index{$key}{$element};
    }
    my $adapter_info=join("|",@adapter);
    my $adapter_info_umi=join("|",@adapter_value);
    print OUT "$key\t$adapter_num\t$adapter_info\t$adapter_info_umi\t$cell_num{$key}\n";
}

close OUT;


mkdir("$output_matrix_dir")unless(-e "$output_matrix_dir");
my $matrix_file="$output_matrix_dir/matrix.mtx.gz";
my $barcode_file="$output_matrix_dir/barcodes.tsv.gz";
my $aptamer_file="$output_matrix_dir/features.tsv.gz";

my @barcode_order=sort { $barcodes{$a} <=> $barcodes{$b} } keys %barcodes;
my @aptamer_order=sort { $aptamer_id{$a} <=> $aptamer_id{$b} } keys %aptamer_id;
open BAR,"| gzip >$barcode_file" or die $!;
foreach my $cells(@barcode_order){
    print BAR "$cells\n";
}
close BAR;

open APTA,"| gzip >$aptamer_file" or die $!;
foreach my $aptamers(@aptamer_order){
    print APTA "$aptamer_info{$aptamers}\n";
}
close APTA;

open MTX,"| gzip >$matrix_file" or die $!;
print MTX "%%MatrixMarket matrix coordinate integer general\n";
print MTX "%metadata_json: {\"software_version\": \"perl\", \"format_version\": 5}\n";
my $final_cell_num=scalar @barcode_order;
my $final_aptamer_num=scalar @aptamer_order;
print MTX "$final_aptamer_num $final_cell_num $all_umi\n";
foreach my $cell_barcode(@barcode_order){
    foreach my $adapter(@aptamer_order){
        if(exists $index{$cell_barcode}{$adapter}){
            print MTX "$aptamer_id{$adapter} $barcodes{$cell_barcode} $index{$cell_barcode}{$adapter}\n";
        }        
    }
}
close MTX;