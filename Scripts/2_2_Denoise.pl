#!/usr/bin/perl
use strict;
use warnings;

#Setting
my ($db, $depth, $size, $identity, $identity2, $dictionary, $family, @setting, $correct);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^DB\s*=\s*(\S+)/){$db = $1;}
	elsif($_ =~ /^Depth\s*=\s*(\S+)/){$depth = $1;}
	elsif($_ =~ /^Length\s*=\s*(\d+)/){$size = $1;}
	elsif($_ =~ /^Correct_error\s*=\s*(\S+)/){$correct = $1;}
	elsif($_ =~ /^UIdentity\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /^LIdentity\s*=\s*(\S+)/){$identity2 = $1;}
	elsif($_ =~ /^Dictionary\s*=\s*(\S+)/){$dictionary = $1;}
	elsif($_ =~ /^Family\s*=\s*(\S+)/){$family = $1;}
	push(@setting, $_);
}
close(SET);
unless($db){print "Error!! Check the DB name in Setting.txt.\n"; exit;}
unless($depth){print "Error!! Check the number of Depth in Setting.txt.\n"; exit;}
unless($size){print "Error!! Check the number of Length in Setting.txt.\n"; exit;}
unless($identity){print "Error!! Check the upper cut line of Homology Search in Setting.txt.\n"; exit;}
unless($identity2){print "Error!! Check the lower cut line of Homology Search in Setting.txt.\n"; exit;}
if($dictionary =~ /^no$/i){undef($dictionary);}
if($family =~ /^no$/i){undef($family);}
if($correct =~ /^no$/i){undef($correct);}

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
opendir (DIR, ".\/Results\/2_1_Find_unique") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_uniques.fa/){$file{$1}++;}
}
closedir DIR;

print "============================================================\n";
print "                        2_2_Denoise                         \n";
print "============================================================\n";

#Denoise
mkdir ".\/Results\/2_2_Denoise";
foreach(sort keys %file){
	my $file = $_;
	print "\n\n$_\n";
	
	open (DATA, "<", ".\/Results\/2_1_Find_unique\/${file}_uniques.fa") or die("error:$!");
	my $count = 0;
	my $sizecount = 0;
	while(<DATA>){if($_ =~ /^>.+size=(\d+)/){$count++; $sizecount += $1;}}
	close(DATA);
	
	if($count < 1){
		open (DATA, ">", ".\/Results\/2_2_Denoise\/${file}_unoise3_result.txt") or die("error:$!");
		close(DATA);
		open (DATA, ">", ".\/Results\/2_2_Denoise\/${file}_zotu.fa") or die("error:$!");
		close(DATA);
		print "$file is no sequences.\n";
		next;
	}
	my $tempdepth = $depth;
	if($depth =~ /(\S+)%/){$tempdepth = int($1/100 * $sizecount);}
	
	my $comand = ".\/Tools\/$usearch -unoise3 \".\/Results\/2_1_Find_unique\/${file}_uniques.fa\" -minsize $tempdepth -ampout \".\/Results\/2_2_Denoise\/${file}_zotu.fa\" -tabbedout \".\/Results\/2_2_Denoise\/${file}_unoise3_result.txt\"";
	system $comand;
	
	#cut low depth seq from unique.fa files
	my ($data, $dseq, @lead);
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
	open (OUT, ">", ".\/Results\/2_1_Find_unique/${file}_uniques.fa") or die("error:$!");
	foreach(@lead){
		if($_ =~ /size=(\d+)/ and $1 >= $tempdepth){print OUT ">$_\n";}
	}
	close(OUT);
	
	#correct_error
	if($correct){
		open (DATA, "<", ".\/Results\/2_2_Denoise\/${file}_unoise3_result.txt") or die("error:$!");
		my %count;
		while(<DATA>){
			if($_ =~ /amp/){
				$_ =~ /(Uniq\d+)/; my $name = $1;
				$_ =~ /size=(\d+)/; my $size = $1;
				$count{$name} = $size;
			}elsif($_ =~ /dqt/){
				$_ =~ /top=(Uniq\d+)/; my $name = $1;
				$_ =~ /size=(\d+)/; my $size = $1;
				$count{$name} += $size;
			}elsif($_ =~ /chfilter/){last;}
		}
		close(DATA);
		
		open (DATA, "<", ".\/Results\/2_2_Denoise\/${file}_zotu.fa") or die("error:$!");
		my @tempdata;
		while(<DATA>){
			if($_ =~ /^>/){
				$_ =~ /(Uniq\d+)/; my $name = $1;
				if($count{$name}){$_ =~ s/size=\d+/size=$count{$name}/;}
				push(@tempdata, $_);
			}else{push(@tempdata, $_);}
		}
		close(DATA);
		
		open (OUT, ">", ".\/Results\/2_2_Denoise\/${file}_zotu.fa") or die("error:$!");
		foreach(@tempdata){print OUT "$_";}
		close(OUT);
	}
}
