#!/usr/bin/perl
use strict;
use warnings;
use List::Util;

#Setting
my ($rarefy, $timing, $temporary);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Rarefaction\s*=\s*(\S+)/){$rarefy = $1;}
	elsif($_ =~ /^Timing\s*=\s*(\d+)/){$timing = $1;}
	elsif($_ =~ /^Temporary\s*=\s*(\S+)/){$temporary = $1;}
}
close(SET);
if($rarefy =~ /^no$/i){exit;}
unless($temporary){$temporary = "YES";}
unless($timing == 1){exit;}

#Get data names
opendir (DIR, ".\/Results\/1_3_Quality_filter") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_filtered.fa/){$file{$1}++;}
}
closedir DIR;

print "============================================================\n";
print "                      1_4_Rarefaction                       \n";
print "============================================================\n";


#Rarefy
mkdir ".\/Results\/1_4_Rarefaction";

my %log;
open (LOG, "<", ".\/Results\/log.txt") or die("error:$!");
my (@logtxt, %raw);
my $min = 0;
while(<LOG>){
	chomp($_);
	$_ =~ s/\r//g;
	if($rarefy =~ /^min$/i){
		my @temp = split(/\t/, $_);
		if($temp[3]){
			$temp[3] =~ s/\s\(.+\)//;
			if($min == 0){$min = $temp[3];}
			if($temp[3] > 0 and $temp[3] < $min){$min = $temp[3];}
		}
	}
	push(@logtxt, $_);
	if($_ =~ /^([^\t]+)\t(\d+)/){$raw{$1} = $2;}
}
close(LOG);
if($rarefy =~ /^min$/i){$rarefy = $min;}

foreach(sort keys %file){
	my $file = $_;
	my ($data, $dseq, @lead);
	open (DATA, "<", ".\/Results\/1_3_Quality_filter\/${file}_filtered.fa") or die("error:$!");
	while(<DATA>){
		chomp($_);
		$_ =~ s/\r//g;
		unless($_ =~ /^\n/){
			if($_ =~ /^>(.+:\d+)/){
				if($data){
					push(@lead, "$data\n$dseq"); $data = $1; undef($dseq);
				}else{$data = $1;}
			}else{
				if($dseq){$dseq = $dseq . $_;}
				else{$dseq = $_;}
			}
		}
	}
	close(DATA);
	push(@lead, "$data\n$dseq");

	open (OUT, ">", ".\/Results\/1_4_Rarefaction\/${file}_rarefy.fa") or die("error:$!");
	my $length = @lead;
	
	srand(1234);
	@lead = List::Util::shuffle @lead;
	@lead = splice (@lead, 0, $rarefy);
	if($length < $rarefy){
		my $per = sprintf("%.2f", $length/$raw{$file}*100);
		$log{$file} = "\t$length ($per%)";
		print "$file\t$length Reads\n\tWarning!! Reads less than $rarefy (Rarefaction)\n";
	}else{
		my $per = sprintf("%.2f", $rarefy/$raw{$file}*100);
		$log{$file} = "\t$rarefy ($per%)";
		print "$file\t$rarefy Reads\n";
	}
	foreach(@lead){print OUT ">$_\n";}
	close(OUT);
	undef($data);undef($dseq);undef(@lead);
	if($temporary =~ /yes/i){unlink ".\/Results\/1_3_Quality_filter\/${file}_filtered.fa"; unlink ".\/Results\/1_3_Quality_filter\/${file}_log.txt";}
}
if($temporary =~ /yes/i){rmdir ".\/Results\/1_3_Quality_filter";}

open (LOG, ">", ".\/Results\/log.txt") or die("error:$!");
foreach(@logtxt){
	my $log = $_;
	foreach(keys %log){
		if($log =~ /$_/){$log = $log . $log{$_};}
	}
	print LOG "$log\n";
}
close(LOG);
