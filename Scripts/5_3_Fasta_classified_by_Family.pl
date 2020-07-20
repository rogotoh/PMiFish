#!/usr/bin/perl
use strict;
use warnings;

#Setting
my ($db, $depth, $size, $identity, $identity2, $dictionary, $family, $temporary);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^DB\s*=\s*(\S+)/){$db = $1;}
	elsif($_ =~ /^Depth\s*=\s*(\S+)/){$depth = $1;}
	elsif($_ =~ /^Length\s*=\s*(\d+)/){$size = $1;}
	elsif($_ =~ /^UIdentity\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /^LIdentity\s*=\s*(\S+)/){$identity2 = $1;}
	elsif($_ =~ /^Dictionary\s*=\s*(\S+)/){$dictionary = $1;}
	elsif($_ =~ /^Family\s*=\s*(\S+)/){$family = $1;}
	elsif($_ =~ /^Temporary\s*=\s*(\S+)/){$temporary = $1;}
}
close(SET);
unless($db){print "Error!! Check the DB name in Setting.txt.\n"; exit;}
unless($depth){print "Error!! Check the number of Depth in Setting.txt.\n"; exit;}
unless($size){print "Error!! Check the number of Length in Setting.txt.\n"; exit;}
unless($identity){print "Error!! Check the upper cut line of Homology Search in Setting.txt.\n"; exit;}
unless($identity2){print "Error!! Check the lower cut line of Homology Search in Setting.txt.\n"; exit;}
unless($temporary){$temporary = "YES";}
unless(-f ".\/DataBase\/$db"){print "No DataBase file or DataBase name didn't match the name in Setting.txt\n"; exit;}
unless(-f ".\/Dictionary\/$dictionary"){
	if($dictionary =~ /^no$/i){undef($dictionary);}
	else{print "The file of Dictionary for common name isn't or don't match the name in Setting.txt\n"; exit;}
}
unless(-f ".\/Dictionary\/$family"){
	if($family =~ /^no$/i){exit;}
	else{print "The file of Dictionary for Family name isn't or don't match the name in Setting.txt\n"; exit;}
}

unless(-f ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/merged_seq.fas"){exit;}

print "============================================================\n";
print "               5_3_Fasta_classified_by_Family               \n";
print "============================================================\n";

#import_family_list
my %family;
open(FILE, "<", ".\/Dictionary\/$family") or die ("error:$!");
while(<FILE>){
	chomp($_);
	my @temp = split(/\t/, $_);
	unless($temp[1]){next;}
	unless($temp[1] =~ /^cf$/i){$family{$temp[1]} = $temp[0];}
	else{$temp[2] =~ /([^_]+)/; $family{$1} = $temp[0];}
}
close(FILE);

#Fasta_classified_by_family�̍쐬
mkdir ".\/Results\/5_3_Fasta_classified_by_Family";
my (@fasta, $fasta, %class, @merge);
open(FILE, "<", ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/merged_seq.fas") or die ("error:$!");
while(<FILE>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ /^>/){$fasta = $_;}
	else{$fasta = $fasta . "\n$_"; push(@fasta, $fasta); undef($fasta);}
}
close(FILE);

foreach(@fasta){
	if($_ =~ />Nohit/){$class{"Nohit"}{$_}++; push(@merge, $_);}
	elsif($_ =~ />U${identity}_([^_]+)_([^_]+)/){
		my $temp1 = $1;
		my $temp2 = $2;
		if($temp1 =~ /^cf$/i){
			unless($family{$temp2}){$class{"N_A"}{$_}++; $_ =~ s/>/>N_A_/; push(@merge, $_); next;}
			$class{$family{$temp2}}{$_}++;
			$_ =~ s/>/>$family{$temp2}_/;
			push(@merge, $_);
		}else{
			unless($family{$temp1}){$class{"N_A"}{$_}++; $_ =~ s/>/>N_A_/; push(@merge, $_); next;}
			$class{$family{$temp1}}{$_}++;
			$_ =~ s/>/>$family{$temp1}_/;
			push(@merge, $_);
		}
	}elsif($_ =~ />([^_]+)_([^_]+)/){
		my $temp1 = $1;
		my $temp2 = $2;
		if($temp1 =~ /^cf$/i){
			unless($family{$temp2}){$class{"N_A"}{$_}++; $_ =~ s/>/>N_A_/; push(@merge, $_); next;}
			$class{$family{$temp2}}{$_}++;
			$_ =~ s/>/>$family{$temp2}_/;
			push(@merge, $_);
		}else{
			unless($family{$temp1}){$class{"N_A"}{$_}++; $_ =~ s/>/>N_A_/; push(@merge, $_); next;}
			$class{$family{$temp1}}{$_}++;
			$_ =~ s/>/>$family{$temp1}_/;
			push(@merge, $_);
		}
	}
}
open(FILE, ">", ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/merged_seq_with_family_name.fas") or die ("error:$!");
foreach(sort @merge){$_ =~ s/\'//g; $_ =~ s/,//g; print FILE "$_\n";}
close(FILE);

open(LIST, ">", ".\/Results\/5_3_Fasta_classified_by_Family\/list.txt") or die ("error:$!");
foreach(sort keys %class){
	my $class = $_;
	my $temp = $class{$_};
	my $kazu = keys %$temp;
	my $length = length($class);
	$length = 30 - $length;
	my $name = "${class}.fas";
	for(my $v = 0; $v < $length;$v++){$name = $name." ";}
	
	print "$name$kazu sequences\n";
	print LIST "$name$kazu sequences\n";
	
	open(OUT, ">", ".\/Results\/5_3_Fasta_classified_by_Family\/${class}.fas") or die ("error:$!");
	foreach(sort keys %$temp){print OUT "$_\n";}
	close(OUT);
}
close(LIST);

#Representative_seq_with_family_name.fas�̍쐬
undef(@fasta); undef($fasta); undef(%class); undef(@merge);
open(FILE, "<", ".\/Results\/5_2_Summary_table\/Representative_seq.fas") or die ("error:$!");
while(<FILE>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ /^>/){$fasta = $_;}
	else{$fasta = $fasta . "\n$_"; push(@fasta, $fasta); undef($fasta);}
}
close(FILE);

foreach(@fasta){
	if($_ =~ />Nohit/){$class{"Nohit"}{$_}++; push(@merge, $_);}
	elsif($_ =~ />U${identity}_([^_]+)_([^_]+)/){
		my $temp1 = $1;
		my $temp2 = $2;
		if($temp1 =~ /^cf$/i){
			unless($family{$temp2}){$class{"N_A"}{$_}++; $_ =~ s/>/>N_A_/; push(@merge, $_); next;}
			$class{$family{$temp2}}{$_}++;
			$_ =~ s/>/>$family{$temp2}_/;
			push(@merge, $_);
		}else{
			unless($family{$temp1}){$class{"N_A"}{$_}++; $_ =~ s/>/>N_A_/; push(@merge, $_); next;}
			$class{$family{$temp1}}{$_}++;
			$_ =~ s/>/>$family{$temp1}_/;
			push(@merge, $_);
		}
	}elsif($_ =~ />([^_]+)_([^_]+)/){
		my $temp1 = $1;
		my $temp2 = $2;
		if($temp1 =~ /^cf$/i){
			unless($family{$temp2}){$class{"N_A"}{$_}++; $_ =~ s/>/>N_A_/; push(@merge, $_); next;}
			$class{$family{$temp2}}{$_}++;
			$_ =~ s/>/>$family{$temp2}_/;
			push(@merge, $_);
		}else{
			unless($family{$temp1}){$class{"N_A"}{$_}++; $_ =~ s/>/>N_A_/; push(@merge, $_); next;}
			$class{$family{$temp1}}{$_}++;
			$_ =~ s/>/>$family{$temp1}_/;
			push(@merge, $_);
		}
	}
}
open(FILE, ">", ".\/Results\/5_2_Summary_table\/Representative_seq_with_family_name.fas") or die ("error:$!");
foreach(sort @merge){print FILE "$_\n";}
close(FILE);
