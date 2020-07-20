#!/usr/bin/perl
use strict;
use warnings;

#Setting
my ($rarefy, $temporary, $timing);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Rarefaction\s*=\s*(\S+)/){$rarefy = $1;}
	elsif($_ =~ /^Timing\s*=\s*(\d+)/){$timing = $1;}
	elsif($_ =~ /^Temporary\s*=\s*(\S+)/){$temporary = $1;}
}
close(SET);
if($rarefy =~ /^no$/i){undef($rarefy);}
unless($temporary){$temporary = "YES";}
unless($timing == 1){undef($rarefy);}

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
if($rarefy){opendir (DIR, ".\/Results\/1_4_Rarefaction") or die ("error:$!");}
else{opendir (DIR, ".\/Results\/1_3_Quality_filter") or die ("error:$!");}
my @read = readdir DIR;
my %file;
foreach (@read) {
	if($rarefy){
		if ($_ =~ /(.+)_rarefy.fa/){$file{$1}++;}
	}else{
		if ($_ =~ /(.+)_filtered.fa/){$file{$1}++;}
	}
}
closedir DIR;

print "============================================================\n";
print "                      2_1_Find_unique                       \n";
print "============================================================\n";

#Find_unique_reads
mkdir ".\/Results\/2_1_Find_unique";
foreach(sort keys %file){
	print "\n\n$_\n";
	my $file = $_;
	if($rarefy){
		my $comand = ".\/Tools\/$usearch -fastx_uniques \".\/Results\/1_4_Rarefaction\/${file}_rarefy.fa\" -sizeout -relabel Uniq -fastaout \".\/Results\/2_1_Find_unique\/${file}_uniques.fa\"";
		system $comand;
	}else{
		my $comand = ".\/Tools\/$usearch -fastx_uniques \".\/Results\/1_3_Quality_filter\/${file}_filtered.fa\" -sizeout -relabel Uniq -fastaout \".\/Results\/2_1_Find_unique\/${file}_uniques.fa\"";
		system $comand;
	}
	if($temporary =~ /yes/i){
		if($rarefy){unlink ".\/Results\/1_4_Rarefaction\/${file}_rarefy.fa";}
		else{unlink ".\/Results\/1_3_Quality_filter\/${file}_filtered.fa"; unlink ".\/Results\/1_3_Quality_filter\/${file}_log.txt";}
	}
}
if($temporary =~ /yes/i){
	if($rarefy){rmdir ".\/Results\/1_4_Rarefaction";}
	else{rmdir ".\/Results\/1_3_Quality_filter";}
}
