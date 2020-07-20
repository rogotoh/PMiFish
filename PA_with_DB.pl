#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;
#version 1.2 (2019/02/13)

#Setting
my ($db, $primer, $diff, $separate, $depth, $size, $identity, $identity2, $dictionary, $family, $temporary, $compress, $trim);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^DB\s*=\s*(\S+)/){$db = $1;}
	elsif($_ =~ /^Primers\s*=\s*(\S+)/){my $temp = $1;if($temp =~ /^no$/i){$primer = "No";}else{$primer = $temp;}}
	elsif($_ =~ /^MaxDiff\s*=\s*(\S+)/){$diff = $1;}
	elsif($_ =~ /^Divide\s*=\s*(\S+)/){$separate = $1;}
	elsif($_ =~ /^Depth\s*=\s*(\S+)/){$depth = $1;}
	elsif($_ =~ /^Length\s*=\s*(\d+)/){$size = $1;}
	elsif($_ =~ /^UIdentity\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /^LIdentity\s*=\s*(\S+)/){$identity2 = $1;}
	elsif($_ =~ /^Dictionary\s*=\s*(\S+)/){$dictionary = $1;}
	elsif($_ =~ /^Family\s*=\s*(\S+)/){$family = $1;}
	elsif($_ =~ /^Temporary\s*=\s*(\S+)/){$temporary = $1;}
	elsif($_ =~ /^Compress\s*=\s*(\S+)/){$compress = $1;}
}
close(SET);
if($diff =~ /^length$/i){$trim = 1;}
unless($diff =~ /\d+/){$diff = 2;}
unless($separate){$separate = 0;}
if($separate =~ /^yes$/i){$separate = 1;}
else{$separate = 0;}
unless($db){print "Error!! Check the DB name in Setting.txt.\n"; exit;}
unless($primer){print "Error!! Check the Primer file name in Setting.txt.\n"; exit;}
unless($depth){print "Error!! Check the number of Depth in Setting.txt.\n"; exit;}
unless($size){print "Error!! Check the number of Length in Setting.txt.\n"; exit;}
unless($identity){print "Error!! Check the upper cut line of Homology Search in Setting.txt.\n"; exit;}
unless($identity2){print "Error!! Check the lower cut line of Homology Search in Setting.txt.\n"; exit;}
unless($temporary){$temporary = "YES";}
unless($compress){$compress = "YES";}
unless(-f ".\/DataBase\/$db"){print "No DataBase file or DataBase name didn't match in Setting.txt\n"; exit;}
unless($primer =~ /^no$/i){
	unless(-f ".\/DataBase\/$primer"){print "No Primers file or file name didn't match in Setting.txt\n"; exit;}
}
unless(-f ".\/Dictionary\/$dictionary"){
	if($dictionary =~ /^no$/i){undef($dictionary);}
	else{print "The file of Dictionary for common name isn't or don't match the name in Setting.txt\n"; exit;}
}
unless(-f ".\/Dictionary\/$family"){
	if($family =~ /^no$/i){undef($family);}
	else{print "The file of Dictionary for Family name isn't or don't match the name in Setting.txt\n"; exit;}
}

unless(-f ".\/Results\/5_2_Summary_table\/Representative_seq_with_family_name.fas"){print "This scripts need \"Representative_seq_with_family_name.fas\" in 5_2_Summary_table folder.\n"; exit;}

print "============================================================\n";
print "             Phylogenetic Analysis with Database            \n";
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

#Representative_seq_with_family_name.fasからデータ収集
mkdir ".\/Results\/5_2_Summary_table\/PA_with_DB";
my ($familyn, $dataname, %data);
open(FILE, "<", ".\/Results\/5_2_Summary_table\/Representative_seq_with_family_name.fas") or die ("error:$!");
while(<FILE>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ />N_A_/){last;}
	if($_ =~ /^>([^_]+_[^_]+)/){$familyn = $1; $dataname = $_;}
	else{$data{$familyn}{$dataname} = $_;}
}
close(FILE);
undef($familyn); undef($dataname);

#Databaseからデータ収集
my (%database);
open(FILE, "<", ".\/Database\/$db") or die ("error:$!");
while(<FILE>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ /^>/){
		$_ =~ s/__/_/g;
		$_ =~ /.+\|([^_]+)_([^_]+)/;
		my $temp1 = $1;
		my $temp2 = $2;
		if($temp1 =~ /^cf$/i){
			unless($family{$temp2}){$familyn = "N_A";}
			else{$familyn = $family{$temp2};}
			$dataname = $_;
		}else{
			unless($family{$temp1}){$familyn = "N_A";}
			else{$familyn = $family{$temp1};}
			$dataname = $_;
		}
	}else{$database{$familyn}{$dataname} = $_;}
}
close(FILE);

#output
foreach(sort keys %data){
	my $count = 0;
	my $familyname = $_;
	open(FILE, ">", ".\/Results\/5_2_Summary_table\/PA_with_DB\/${familyname}.fas") or die ("error:$!");
	my $temp = $data{$_};
	foreach(sort keys %$temp){
		my $filename = $_; $count++;
		my $tempname = $_;
		$tempname =~ s/\'//g;
		$tempname =~ s/,//g;
		$tempname =~ s/&//g;
		print FILE "$tempname\n$data{$familyname}{$filename}\n";
	}
	my $temp2 = $database{$familyname};
	
	#同じ名前かつ同じ配列を持つものを除去
	my (@check, %check);
	foreach(sort keys %$temp2){
		my $filename = $_;
		my $cutname = $_;
		$cutname =~ s/.+\|//;
		unless($check{$database{$familyname}{$filename}}){$check{$database{$familyname}{$filename}}{$cutname}++; push(@check, $filename);}
		else{
			unless($check{$database{$familyname}{$filename}}{$cutname}){$check{$database{$familyname}{$filename}}{$cutname}++; push(@check, $filename);}
		}
	}
	foreach(sort @check){
		my $filename = $_; $count++;
		my $tempname = $_;
		$tempname =~ s/\'//g;
		$tempname =~ s/,//g;
		$tempname =~ s/&//g;
		print FILE "$tempname\n$database{$familyname}{$filename}\n";
	}
	close(FILE);
	
	my $length = length($familyname);
	$length = 30 - $length;
	my $name = "${familyname}.fas";
	for(my $v = 0; $v < $length;$v++){$name = $name." ";}
	print "$name$count sequences\n";
}


#Phylogenetic_trees_with_Database
print "============================================================\n";
print "              Phylogenetic_trees_with_Database              \n";
print "============================================================\n";

#Get data names
opendir (DIR, ".\/Results\/5_2_Summary_table\/PA_with_DB") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+).fas/){$file{$1}++;}
}
closedir DIR;

my @check;
foreach(sort keys %file){
	my $file = $_;
	open (DATA, "<", ".\/Results\/5_2_Summary_table\/PA_with_DB\/${file}.fas") or die("error:$!");
	my $count = 0;
	while(<DATA>){
		if($_ =~ /^>/){$count++;}
	}
	if($count > 3){push(@check, $file);}
}

#Phylogenetic_treesの作成
foreach (@check){
	my $file = $_;
	my $command = "megacc -a \".\/Tools\/muscle_align_nucleotide.mao\" -d \".\/Results\/5_2_Summary_table\/PA_with_DB\/${file}.fas\" -o $file";
	system $command;
	move (".\/${file}.meg", ".\/Results\/5_2_Summary_table\/PA_with_DB");
	unlink ".\/${file}_summary.txt";
	
	$command = "megacc -a \".\/Tools\/infer_NJ_nucleotide.mao\" -d \".\/Results\/5_2_Summary_table\/PA_with_DB\/${file}.meg\" -o $file";
	system $command;
	move (".\/${file}.nwk", ".\/Results\/5_2_Summary_table\/PA_with_DB");
	unlink ".\/${file}_consensus.nwk";
	unlink ".\/${file}_summary.txt";
	unlink ".\/${file}_partitions.txt";
}

