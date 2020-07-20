#!/usr/bin/perl
use strict;
use warnings;

#2018/12/04 1_2 Strip_primers change to use the Primer_Cleaner.pl

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

#primer”z—ñ‚ÌŽæ“¾
open (DATA2, "<", ".\/DataBase\/$primer") or die("error:$!");
my $countp = 0;
my (%forward, %reverse, $namep);
while(<DATA2>){
	if($_ =~ /\#/){next;}
	chomp($_);
	if($_ =~ /Forward/){$countp = 1;next;}
	if($_ =~ /Reverse/){$countp = 2;next;}
	if($countp == 1){
		if($_ =~ /^>(.+)/){$namep = $1;}
		elsif($_ =~ /^[a-z]/i){$_ =~ tr/[a-z]/[A-Z]/; $forward{$namep} = $_;}
	}
	if($countp == 2){
		if($_ =~ /^>(.+)/){$namep = $1;}
		elsif($_ =~ /^[a-z]/i){
			my $rev = reverse($_);
			$rev =~ tr/ATGCURYMKDHBV/TACGAYRKMHDVB/;
			$reverse{$namep} = $rev;}
	}
}
close(DATA2);
unless(%forward){print "Error!! No primer seqeunce in $primer\n"; exit;}
unless(%reverse){print "Error!! No primer seqeunce in $primer\n"; exit;}

my ($primerF, $primerR);
my $primer_num = keys %forward;
if($primer_num == 1){$separate = 0;}
if($trim and $primer_num == 1){
	foreach(keys %forward){$primerF = length($forward{$_}); $primerR = length($reverse{$_});}
}elsif($trim and $primer_num > 1){
	print "Error! \"MaxDiff = length\" in Setting.txt can use when primer pair is one pair\n"; exit;
}
if($trim){$separate = 0;}

#Usearch check
opendir (DIR, ".\/Tools") or die ("error:$!");
my @tool = readdir DIR;
my $usearch;
foreach (@tool) {
	if ($_ =~ /(usearch.+)/){$usearch = $1;}
}
closedir DIR;
unless($usearch){print "Error!! Please put an usearch exectable file in Tools directory!\n"; exit;}


#Get data names
opendir (DIR, ".\/Results\/1_1_Merge_paird_reads") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_assembled_seq/){$file{$1}++;}
}
closedir DIR;

#Remove_primer
print "============================================================\n";
print "                     1_2_Strip_primers                      \n";
print "============================================================\n";

unless($primer =~ /^no$/i){mkdir ".\/Results\/1_2_Strip_primers";}
foreach(sort keys %file){
	print "\n\n$_\n";
	my $file = $_;
	my $comand;
	unless($primer =~ /^no$/i){
		print "Now processing...\n";
		if($trim){
			$comand = ".\/Tools\/$usearch -fastx_truncate \".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq\" -stripleft $primerF -stripright $primerR -fastqout \".\/Results\/1_2_Strip_primers\/${file}_stripped.fq\"";
		}else{
			$comand = "perl ./Scripts/Primer_Cleaner.pl -f \".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq\" -db \"./DataBase/$primer\" -diff $diff -o \".\/Results\/1_2_Strip_primers\/${file}_stripped.fq\" -s $separate -l";
		}
		system $comand;
		if($temporary =~ /yes/i){unlink ".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq"; unlink ".\/Results\/1_1_Merge_paird_reads\/${file}_log.txt";}
	}else{
		print "The step was skipped (You selected \"Primers = No\" in Setting.txt).\n";
	}
}
if($temporary =~ /yes/i){rmdir ".\/Results\/1_1_Merge_paird_reads";}
