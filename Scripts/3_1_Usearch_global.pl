#!/usr/bin/perl
use strict;
use warnings;

#2018/01/05 add qcovs to outfmt

#Setting
my ($db, $rarefy, $depth, $size, $identity, $identity2, $dictionary, $family, $temporary);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^DB\s*=\s*(\S+)/){$db = $1;}
	elsif($_ =~ /^Rarefaction\s*=\s*(\S+)/){$rarefy = $1;}
	elsif($_ =~ /^Depth\s*=\s*(\S+)/){$depth = $1;}
	elsif($_ =~ /^Length\s*=\s*(\d+)/){$size = $1;}
	elsif($_ =~ /^UIdentity\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /^LIdentity\s*=\s*(\S+)/){$identity2 = $1;}
	elsif($_ =~ /^Dictionary\s*=\s*(\S+)/){$dictionary = $1;}
	elsif($_ =~ /^Family\s*=\s*(\S+)/){$family = $1;}
	elsif($_ =~ /^Temporary\s*=\s*(\S+)/){$temporary = $1;}
}
close(SET);
if($rarefy =~ /^no$/i){undef($rarefy);}
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
	if($family =~ /^no$/i){undef($family);}
	else{print "The file of Dictionary for Family name isn't or don't match the name in Setting.txt\n"; exit;}
}


#Database check
opendir (DIR, ".\/DataBase") or die ("error:$!");
my @database = readdir DIR;
my $dbcheck = 0;
foreach(@database){
	if ($_ =~ /${db}\.udb/){$dbcheck++;}
}
closedir DIR;

opendir (DIR, ".\/Tools") or die ("error:$!");
my @tool = readdir DIR;
my $usearch;
foreach (@tool) {
	if ($_ =~ /(usearch.+)/){$usearch = $1;}
}
closedir DIR;
unless($usearch){print "Error!! Please put an usearch exectable file in Tools directory!\n"; exit;}

unless($dbcheck){
	my $command = ".\/Tools\/$usearch -makeudb_usearch \".\/DataBase\/$db\" -output \".\/DataBase\/${db}\.udb\"";
	system $command;
}

#Get data names
if(-e ".\/Results\/2_4_Rarefaction"){
	opendir (DIR, ".\/Results\/2_4_Rarefaction") or die ("error:$!");
}elsif(-e ".\/Results\/2_3_Separate_chimera"){
	opendir (DIR, ".\/Results\/2_3_Separate_chimera") or die ("error:$!");
}else{
	opendir (DIR, ".\/Results\/2_1_Find_unique") or die ("error:$!");
}
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_zotu_nonchimeras.fa/){$file{$1}++;}
	if ($_ =~ /(.+)_uniques.fa/){$file{$1}++;}
	if ($_ =~ /(.+)_rarefy.fa/){$file{$1}++;}
}
closedir DIR;

unless(-e ".\/Results\/2_3_Separate_chimera"){
	foreach(sort keys %file){
		my ($data, $dseq, @lead);
		my $file = $_;
		open (DATA, "<", ".\/Results\/2_1_Find_unique\/${file}_uniques.fa") or die("error:$!");
		while(<DATA>){
			chomp($_);
			$_ =~ s/\r//g;
			unless($_ =~ /^\n/){
				if($_ =~ /^>(.+\;)/){
					if($data){push(@lead, "$data\n$dseq"); $data = $1; undef($dseq);}
					else{$data = $1;}
				}else{
					if($dseq){$dseq = $dseq . $_;}
					else{$dseq = $_;}
				}
			}
		}
		close(DATA);
		if($data){push(@lead, "$data\n$dseq");}
		
		my $sizecount = 0;
		foreach(@lead){
			if($_ =~ /size=(\d+)/){$sizecount += $1;}
		}
		my $tempdepth = $depth;
		if($depth =~ /(\S+)%/){$tempdepth = int($1/100 * $sizecount);}
		
		open (OUT, ">", ".\/Results\/2_1_Find_unique/${file}_uniques.fa") or die("error:$!");
		foreach(@lead){
			if($_ =~ /size=(\d+)/ and $1 >= $tempdepth){print OUT ">$_\n";}
		}
		close(OUT);
	}
}

print "============================================================\n";
print "                     3_1_Usearch_global                     \n";
print "============================================================\n";

#Usearch_global
mkdir ".\/Results\/3_1_Usearch_global";
$identity2 = $identity2/100;
foreach(sort keys %file){
	print "$_\n";
	my $file = $_;
	if(-f ".\/Results\/2_4_Rarefaction\/${file}_rarefy.fa"){
		my $comand = ".\/Tools\/$usearch -usearch_global \".\/Results\/2_4_Rarefaction\/${file}_rarefy.fa\" -db \".\/DataBase\/${db}.udb\" -id $identity2 -maxaccepts 100 -strand plus -userout \".\/Results\/3_1_Usearch_global\/${file}_usearch_results.txt\" -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+qcov";
		system $comand;
	}elsif(-f ".\/Results\/2_3_Separate_chimera\/${file}_zotu_nonchimeras.fa"){
		my $comand = ".\/Tools\/$usearch -usearch_global \".\/Results\/2_3_Separate_chimera\/${file}_zotu_nonchimeras.fa\" -db \".\/DataBase\/${db}.udb\" -id $identity2 -maxaccepts 100 -strand plus -userout \".\/Results\/3_1_Usearch_global\/${file}_usearch_results.txt\" -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+qcov";
		system $comand;
	}else{
		my $comand = ".\/Tools\/$usearch -usearch_global \".\/Results\/2_1_Find_unique\/${file}_uniques.fa\" -db \".\/DataBase\/${db}.udb\" -id $identity2 -maxaccepts 100 -strand plus -userout \".\/Results\/3_1_Usearch_global\/${file}_usearch_results.txt\" -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+qcov";
		system $comand;
	}
	print "\n";
}
