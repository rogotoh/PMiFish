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
unless($timing == 2){exit;}

#Get data names
opendir (DIR, ".\/Results\/2_3_Separate_chimera") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_zotu_nonchimeras.fa/){$file{$1}++;}
}
closedir DIR;

print "============================================================\n";
print "                      2_4_Rarefaction                       \n";
print "============================================================\n";


#Rarefy
mkdir ".\/Results\/2_4_Rarefaction";

my %log;
open (LOG, "<", ".\/Results\/log.txt") or die("error:$!");
my (@logtxt, %raw);
my $min = 0;
while(<LOG>){
	chomp($_);
	$_ =~ s/\r//g;
	if($rarefy =~ /^min$/i){
		my @temp = split(/\t/, $_);
		if($temp[4]){
			$temp[4] =~ s/\s\(.+\)//;
			if($min == 0){$min = $temp[4];}
			if($temp[4] > 0 and $temp[4] < $min){$min = $temp[4];}
		}
	}
	push(@logtxt, $_);
	if($_ =~ /^([^\t]+)\t(\d+)/){$raw{$1} = $2;}
}
close(LOG);
if($rarefy =~ /^min$/i){$rarefy = $min;}

foreach(sort keys %file){
	my $file = $_;
	my ($data, %seq, $lead);
	open (DATA, "<", ".\/Results\/2_3_Separate_chimera\/${file}_zotu_nonchimeras.fa") or die("error:$!");
	while(<DATA>){
		chomp($_);
		$_ =~ s/\r//g;
		unless($_ =~ /^\n/){
			if($_ =~ /^>Uniq(\d+)\;size=(\d+)/){
				$data = $1;
				unless($lead){$lead = "$data\t" x $2;}
				else{$lead = $lead . "$data\t" x $2;}
			}else{
				$seq{$data} = $_;
			}
		}
	}
	close(DATA);
	$lead =~ s/\t$//;
	my @lead = split(/\t/, $lead);

	open (OUT, ">", ".\/Results\/2_4_Rarefaction\/${file}_rarefy.fa") or die("error:$!");
	my $length = @lead;
	
	srand(123);
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
	my %temp;
	foreach(@lead){$temp{$_}++;}
	foreach(sort {$a <=> $b} keys %temp){print OUT ">Uniq$_;size=$temp{$_};amptype=otu;\n$seq{$_}\n";}
	close(OUT);
	undef($data);undef(%seq);undef($lead);
}

open (LOG, ">", ".\/Results\/log.txt") or die("error:$!");
foreach(@logtxt){
	my $log = $_;
	foreach(keys %log){
		if($log =~ /$_/){$log = $log . $log{$_};}
	}
	print LOG "$log\n";
}
close(LOG);
