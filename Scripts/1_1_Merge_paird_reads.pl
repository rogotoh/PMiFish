#!/usr/bin/perl
use strict;
use warnings;

#2018/12/06 use primer.txt
#2018/01/20 output a log file at 1_1 step
#2018/01/05 add the compress option

#Setting
my ($db, $primer, $diff, $depth, $size, $identity, $identity2, $dictionary, $family, $temporary, $compress);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^DB\s*=\s*(\S+)/){$db = $1;}
	elsif($_ =~ /Primers\s*=\s*(\S+)/){my $temp = $1;if($temp =~ /^no$/i){$primer = "No";}else{$primer = $temp;}}
	elsif($_ =~ /MaxDiff\s*=\s*(\d+)/){$diff = $1;}
	elsif($_ =~ /Depth\s*=\s*(\S+)/){$depth = $1;}
	elsif($_ =~ /Length\s*=\s*(\d+)/){$size = $1;}
	elsif($_ =~ /UIdentity\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /LIdentity\s*=\s*(\S+)/){$identity2 = $1;}
	elsif($_ =~ /Dictionary\s*=\s*(\S+)/){$dictionary = $1;}
	elsif($_ =~ /Family\s*=\s*(\S+)/){$family = $1;}
	elsif($_ =~ /Temporary\s*=\s*(\S+)/){$temporary = $1;}
	elsif($_ =~ /Compress\s*=\s*(\S+)/){$compress = $1;}
}
close(SET);
unless($diff =~ /\d+/){$diff = 2;}
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

#Database and Usearch check
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

#Decompressing files
opendir (DIR, ".\/Run") or die ("error:$!");
my @run = readdir DIR;
my $gzip = 0;
foreach (@run) {
	if ($_ =~ /\.gz$/){$gzip++;}
}
closedir DIR;

if($gzip){
	print "Decompressing gz file to fastq file...\n";
	foreach (@run) {
		if ($_ =~ /\.gz$/){
			my $comand = "gzip -d \".\/Run\/$_\"";
			system $comand;
		}
	}
}

#Get data names
opendir (DIR, ".\/Run") or die ("error:$!");
my @read = readdir DIR;
my $read_count = 0;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_R1/){$file{$1}++; $read_count++;}
	if ($_ =~ /(.+)_R2/){$file{$1}++; $read_count++;}
}
closedir DIR;
foreach(keys %file){
	if($file{$_} < 2){print "Error!! Check the Run directory!The analysis needs both R1 and R2 fastq files.\nR1 or R2 fastq files of $_ is not.\n"; exit;}
}
unless(-d ".\/Results"){mkdir ".\/Results";}
print "============================================================\n";
print "                    1_1_Merge_paird_reads                   \n";
print "============================================================\n";
#assemble pair seqs
mkdir ".\/Results\/1_1_Merge_paird_reads";
foreach(sort keys %file){
	my $file = $_;
	my $file1;
	foreach(sort @read){
		if($_ =~ /${file}_R1[\._]/){$file1 = $_;}
	}
	my $comand = ".\/Tools\/$usearch -fastq_mergepairs \".\/Run\/$file1\" -fastqout \".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq\" -log \".\/Results\/1_1_Merge_paird_reads\/${file}_log.txt\"";
	system $comand;
	if($compress =~ /yes/i){
		print "\nCompressing fastq file to gz file...\n";
		foreach(@read){
			if($_ !~ /gz$/ and $_ =~ /$file/){
				$comand = "gzip \".\/Run\/$_\"";
				system $comand;
			}
		}
	}
}
