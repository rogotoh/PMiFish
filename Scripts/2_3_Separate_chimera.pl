#!/usr/bin/perl
use strict;
use warnings;

#Setting
my ($temporary);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Temporary\s*=\s*(\S+)/){$temporary = $1;}
}
close(SET);
unless($temporary){$temporary = "YES";}

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
opendir (DIR, ".\/Results\/2_2_Denoise") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_zotu.fa/){$file{$1}++;}
}
closedir DIR;


print "============================================================\n";
print "                    2_3_Separate_chimera                    \n";
print "============================================================\n";

#Separate_chimera
mkdir ".\/Results\/2_3_Separate_chimera";

my %log;
my (@logtxt, %raw);
if(-f ".\/Results\/log.txt"){
	open (LOG, "<", ".\/Results\/log.txt") or die("error:$!");
	while(<LOG>){
		push(@logtxt, $_);
		if($_ =~ /^([^\t]+)\t(\d+)/){$raw{$1} = $2;}
	}
	close(LOG);
}

foreach(sort keys %file){
	my $file = $_;
	#Separate chimeras
	my ($data, $dseq, @lead);
	open (DATA, "<", ".\/Results\/2_2_Denoise\/${file}_zotu.fa") or die("error:$!");
	open (OUT1, ">", ".\/Results\/2_3_Separate_chimera\/${file}_zotu_chimeras.fa") or die("error:$!");
	open (OUT2, ">", ".\/Results\/2_3_Separate_chimera\/${file}_zotu_nonchimeras.fa") or die("error:$!");
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
	
	my $count1 = 0;
	my $count2 = 0;
	my $read = 0;
	foreach(@lead){
		if($_ =~ /chimera/){print OUT1 ">$_\n";$count1++;}
		else{
			$_ =~ /size=(\d+)/;
			$read = $read + $1;
			print OUT2 ">$_\n";$count2++;
		}
	}
	if(-f ".\/Results\/log.txt"){
		my $per = sprintf("%.2f", $read/$raw{$file}*100);
		$log{$file} = "\t$read ($per%)";
	}
	print "$file\n\t$count2 non-chimeras\n\t$count1 chimeras\n";
	close(OUT1); close(OUT2);
	undef($data);undef($dseq);undef(@lead);
	if($temporary =~ /yes/i){unlink ".\/Results\/2_2_Denoise\/${file}_zotu.fa";}
}

if(-f ".\/Results\/log.txt"){
	open (LOG, ">", ".\/Results\/log.txt") or die("error:$!");
	foreach(@logtxt){
		chomp($_);
		$_ =~ s/\r//g;
		my $log = $_;
		foreach(keys %log){
			if($log =~ /$_\t/){$log = $log . $log{$_};}
		}
		print LOG "$log\n";
	}
	close(LOG);
}

