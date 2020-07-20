#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;

#Setting
my ($phy, $temporary);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Phylogenetic\s*=\s*(\S+)/){$phy = $1;}
}
close(SET);
if($phy =~ /^no$/i){exit;}
unless($phy){$phy = "YES";}

print "============================================================\n";
print "                   6_2_Phylogenetic_trees                   \n";
print "============================================================\n";

#Get data names
opendir (DIR, ".\/Results\/5_3_Fasta_classified_by_Family") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+).fas/){$file{$1}++;}
}
closedir DIR;

my @check;
foreach(sort keys %file){
	my $file = $_;
	open (DATA, "<", ".\/Results\/5_3_Fasta_classified_by_Family\/${file}.fas") or die("error:$!");
	my $count = 0;
	while(<DATA>){
		if($_ =~ /^>/){$count++;}
	}
	if($count > 3){push(@check, $file);}
}

#Phylogenetic_trees�̍쐬
mkdir ".\/Results\/6_2_Phylogenetic_trees";
foreach (@check){
	my $file = $_;
	my $command = "megacc -a \".\/Tools\/muscle_align_nucleotide.mao\" -d \".\/Results\/5_3_Fasta_classified_by_Family\/${file}.fas\" -o $file";
	system $command;
	move (".\/${file}.meg", ".\/Results\/6_2_Phylogenetic_trees");
	unlink ".\/${file}_summary.txt";
	
	$command = "megacc -a \".\/Tools\/infer_NJ_nucleotide.mao\" -d \".\/Results\/6_2_Phylogenetic_trees\/${file}.meg\" -o $file";
	system $command;
	move (".\/${file}.nwk", ".\/Results\/6_2_Phylogenetic_trees");
	unlink ".\/${file}_consensus.nwk";
	unlink ".\/${file}_summary.txt";
	unlink ".\/${file}_partitions.txt";
}
