#!/usr/bin/perl
use strict;
use warnings;

#2018/12/04 1_2 Strip_primers change to use the Primer_Cleaner.pl

#Setting
my ($db, $primer, $diff, $depth, $size, $identity, $identity2, $dictionary, $family, $temporary, $compress);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^DB\s*=\s*(\S+)/){$db = $1;}
	elsif($_ =~ /^Primers\s*=\s*(\S+)/){my $temp = $1;if($temp =~ /^no$/i){$primer = "No";}else{$primer = $temp;}}
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
if($primer !~ /^no$/i){
	opendir (DIR, ".\/Results\/1_2_Strip_primers") or die ("error:$!");
}else{
	opendir (DIR, ".\/Results\/1_1_Merge_paird_reads") or die ("error:$!");
}
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($primer !~ /^no$/i and $_ =~ /(.+)_stripped/){$file{$1}++;}
	elsif($_ =~ /(.+)_assembled/){$file{$1}++;}
}
closedir DIR;

print "============================================================\n";
print "                    1_3_Quality_filter                      \n";
print "============================================================\n";

#Quality_filter
mkdir ".\/Results\/1_3_Quality_filter";
foreach(sort keys %file){
	print "\n\n$_\n";
	my $file = $_;
	if($primer !~ /^no$/i){
		my $comand = ".\/Tools\/$usearch -fastq_filter \".\/Results\/1_2_Strip_primers\/${file}_stripped.fq\" -fastq_maxee 1.0 -fastq_minlen $size -fastaout \".\/Results\/1_3_Quality_filter\/${file}_filtered.fa\"";
		system $comand;
		if($temporary =~ /yes/i){unlink ".\/Results\/1_2_Strip_primers\/${file}_stripped.fq";}
	}else{
		my $comand = ".\/Tools\/$usearch -fastq_filter \".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq\" -fastq_maxee 1.0 -fastq_minlen $size -fastaout \".\/Results\/1_3_Quality_filter\/${file}_filtered.fa\" -log \".\/Results\/1_3_Quality_filter\/${file}_log.txt\"";
		system $comand;
		if($temporary =~ /yes/i){unlink ".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq"; unlink ".\/Results\/1_1_Merge_paird_reads\/${file}_log.txt";}
	}
}
if($primer !~ /^no$/i and $temporary =~ /yes/i){rmdir ".\/Results\/1_2_Strip_primers";}
else{rmdir ".\/Results\/1_1_Merge_paird_reads";}
