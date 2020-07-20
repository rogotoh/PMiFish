#!/usr/bin/perl
use strict;
use warnings;

#Get data names
opendir (DIR, ".\/Results\/4_1_Annotation") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_Representative_seq/){$file{$1}++;}
}
closedir DIR;

print "============================================================\n";
print "            5_1_Fasta_for_Phylogenetic_Analysis             \n";
print "============================================================\n";

#Annotation
mkdir ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis";
my (%merged, %all);
foreach(sort keys %file){
	my $file = $_;
	open (DATA, "<", "./Results/4_1_Annotation\/${file}_Representative_seq.fas") or die("error:$!");
	my $fname;
	while(<DATA>){
		chomp($_);
		if($_ =~ /^>(.+)/){$fname = $1 . "_$file";}
		else{$merged{$_}{$fname}++; $all{$fname} = $_;}
	}
	close(DATA);
}

open (OUT, ">", ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/all_representative_seqs.fas") or die("error:$!");
foreach(sort keys %all){print OUT ">$_\n$all{$_}\n";}
close(OUT);

open (OUT2, ">", ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/merged_seq.fas") or die("error:$!");
open (OUT3, ">", ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/merged_list.txt") or die("error:$!");

my (@narabi1, @narabi2, %haplo);
foreach(keys %merged){
	my $seq = $_;
	my $hash = $merged{$_};
	my %names = %$hash;
	my @keys = keys %names;
	my @narabikae;
	foreach(sort keys %file){
		my $itiji = $_;
		foreach(@keys){
			if($_ =~ /$itiji/){push (@narabikae, $_); last;}
		}
	}
	my $num = @keys;
	if($num > 1){
		my ($temp, $temp2);
		$keys[0] =~ /(.+)_Uniq\d+_\d+_reads/;
		$temp = $1;
		$temp =~ s/_otu\d+//;
		$haplo{$temp}++;
		if($haplo{$temp} > 1){
			push(@narabi1, ">${temp}_h$haplo{$temp}_from_${num}_sites\n$seq\n");
			$temp2 = ">${temp}_h$haplo{$temp}_from_${num}_sites\n";
		}else{
			push(@narabi1, ">${temp}_from_${num}_sites\n$seq\n");
			$temp2 = ">${temp}_from_${num}_sites\n";
		}
		foreach(@narabikae){$temp2 = $temp2 . "\t$_\n";}
		push(@narabi2, $temp2);
	}else{
		push(@narabi1, ">$keys[0]\n$seq\n");
	}
}
foreach(sort @narabi1){print OUT2 "$_";}
foreach(sort @narabi2){print OUT3 "$_";}
my $count = @narabi1;
print "$count Mereged Sequences\n";

close(OUT2);
close(OUT3);
