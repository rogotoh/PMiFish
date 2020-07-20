#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min);

#Setting
my ($db, $depth, $size, $identity, $identity2, $dictionary, $family, $habitat, $temporary);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^DB\s*=\s*(\S+)/){$db = $1;}
	elsif($_ =~ /^Depth\s*=\s*(\S+)/){$depth = $1;}
	elsif($_ =~ /^Length\s*=\s*(\d+)/){$size = $1;}
	elsif($_ =~ /^UIdentity\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /^LIdentity\s*=\s*(\S+)/){$identity2 = $1;}
	elsif($_ =~ /^Dictionary\s*=\s*(\S+)/){$dictionary = $1;}
	elsif($_ =~ /^Family\s*=\s*(\S+)/){$family = $1;}
	elsif($_ =~ /^Habitat\s*=\s*(\S+)/){$habitat = $1;}
	elsif($_ =~ /^Temporary\s*=\s*(\S+)/){$temporary = $1;}
}
close(SET);
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
if($habitat){
	unless(-f ".\/Dictionary\/$habitat"){
		if($habitat =~ /^no$/i){undef($habitat);}
		else{print "The file of Dictionary for habitat isn't or don't match the name in Setting.txt\n"; exit;}
	}
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
opendir (DIR, ".\/Results\/4_1_Annotation") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_Representative_seq/){$file{$1}++;}
}
closedir DIR;

print "============================================================\n";
print "                     5_2_Summary_table                      \n";
print "============================================================\n";

#Summary_table
mkdir ".\/Results\/5_2_Summary_table";
my (%table, %under, $spname, $title, %phy, $title2, $spname2, %replist);
foreach(sort keys %file){
	my $file = $_;
	open (DATA, "<", ".\/Results\/4_1_Annotation\/${file}_Representative_seq.fas") or die("error:$!");
	while(<DATA>){
		if($_ =~ /^>U${identity}_(.+)_Uniq\d+/){
			$title = $_;
			$spname = $1;
			$spname =~ s/_otu\d+//;
			$spname =~ s/_/ /g;
			next;
		}
		if($_ =~ /^>Nohit/){next;}
		if($_ =~ /^>(.+)_Uniq\d+_(\d+)_reads/){
			$title2 = $_;
			$spname2 = $1;
			my $temp1 = $1;
			my $temp2 = $2;
			$temp1 =~ s/_/ /g;
			$table{$temp1}{$file} = $temp2;
			chomp($title2); chomp($_);
			$title2 =~ s/\r//; $_ =~ s/\r//;
			unless($replist{$spname2}){$replist{$spname2} = $title2;}
			next;
		}
		if($title){
			chomp($title); chomp($_);
			$title =~ s/\r//; $_ =~ s/\r//;
			$under{"U${identity}_$spname"}{"${title}_$file"} = $_;
			undef($title);undef($spname);
		}
		if($title2){
			chomp($_);
			$_ =~ s/\r//;
			unless($phy{$spname2}){$phy{$spname2} = $_;}
			undef($title2);undef($spname2);
		}
	}
	close(DATA);
}

# make global otus;
my %seiri;
foreach (sort {$a cmp $b} keys %under){
	my $spname = $_;
	my $title = $under{$spname};
	my $seqnum = keys %$title;
	if($seqnum > 1){
		my @keys = keys %$title;
		my (%narabi, %tempname, %modosi);
		my $tempcount = 0;
		foreach(@keys){$_ =~ /(\d+)_reads/; $narabi{$1}{$_}++; $tempcount++; $tempname{$_} = "Uniq$tempcount;size=$1;"; $modosi{"Uniq$tempcount;size=$1;"} = $_;}
		undef(@keys);
		open (TEMP, ">", ".\/Results\/5_2_Summary_table\/temp.fas") or die("error:$!");
		foreach(sort {$b <=> $a} keys %narabi){
			my $suji = $narabi{$_};
			foreach(sort keys %$suji){print TEMP ">$tempname{$_}\n$under{$spname}{$_}\n";}
		}
		close(TEMP);
		
		my $tempid = $identity/100;
		my $comand = ".\/Tools\/$usearch -cluster_smallmem \".\/Results\/5_2_Summary_table\/temp.fas\" -id $tempid -uc \".\/Results\/5_2_Summary_table\/temp.uc\" -sortedby size -quiet";
		system $comand;
		
		open (TEMP, "<", ".\/Results\/5_2_Summary_table\/temp.uc") or die("error:$!");
		my $sp = 0;
		my %otus;
		while(<TEMP>){
			chomp($_);
			$_ =~ s/\r//g;
			my @temp = split(/\t/, $_);
			if($temp[0] eq "S"){
				$sp++;
				$seiri{"${spname} gotu$sp"}{$modosi{$temp[8]}}++;
				$otus{$temp[8]} = "${spname} gotu$sp";
				$phy{"${spname} gotu$sp"} = $under{$spname}{$modosi{$temp[8]}};
				$replist{"${spname} gotu$sp"} = $modosi{$temp[8]};
			}elsif($temp[0] eq "H"){
				$seiri{$otus{$temp[9]}}{$modosi{$temp[8]}}++;
				$replist{$otus{$temp[9]}} = $modosi{$temp[8]};
			}
		}
		close(TEMP);
		unlink ".\/Results\/5_2_Summary_table\/temp.fas";
		unlink ".\/Results\/5_2_Summary_table\/temp.uc";
	}else{
		my @keys = keys %$title;
		$seiri{$spname}{$keys[0]}++;
		unless($phy{$spname}){$phy{$spname} = $under{$spname}{$keys[0]};}
		unless($replist{$spname}){$replist{$spname} = $keys[0];}
	}
}

#output cluster list
open(OUT, ">", ".\/Results\/5_2_Summary_table\/cluster_list_of_U${identity}.txt") or die ("error:$!");
foreach(sort {$a cmp $b} keys %seiri){
	my $temp = $seiri{$_};
	my $name = $_;
	print OUT "$_\n";
	foreach(keys %$temp){
		print OUT "\t$_\n";
		$_ =~ /(\d+)_reads_(.+)$/;
		$table{$name}{$2} = $1;
	}
}
close(OUT);

#output Representative_seq.fas
open(OUT, ">", ".\/Results\/5_2_Summary_table\/Representative_seq.fas") or die ("error:$!");
foreach(sort {$a cmp $b} keys %phy){
	my $name = $_;
	$name =~ s/ /_/g;
	print OUT ">$name\n$phy{$_}\n";
}
close(OUT);

#output Representative_list.txt
open(OUT, ">", ".\/Results\/5_2_Summary_table\/Representative_list.txt") or die ("error:$!");
foreach(sort {$a cmp $b} keys %replist){print OUT "$_\n\t$replist{$_}\n";}
close(OUT);

#import_family_list
my %family;
if($family){
	open(FILE, "<", ".\/Dictionary\/$family") or die ("error:$!");
	while(<FILE>){
		chomp($_);
		my @temp = split(/\t/, $_);
		unless($temp[1]){next;}
		unless($temp[1] =~ /^cf$/i){$family{$temp[1]} = $temp[0];}
		else{$temp[2] =~ /([^_]+)/; $family{$1} = $temp[0];}
	}
	close(FILE);
}


#import habitat data
my %habitat;
if($habitat){
	open(FILE, "<", ".\/Dictionary\/$habitat") or die ("error:$!");
	while(<FILE>){
		if($_ =~ /^Genus\t/){next;}
		$_ =~ s/\r//g;
		chomp($_);
		my @temp = split(/\t/, $_);
		my $water;
		if($temp[2] == -1){$water = "F";}
		if($temp[3] == -1){
			unless($water){$water = "B";}
			else{$water = $water . "B";}
		}
		if($temp[4] == -1){
			unless($water){$water = "S";}
			else{$water = $water . "S";}
		}
		unless($water){$water = "NA";}
		$habitat{"$temp[0] $temp[1]"} = "$water\t$temp[5]\t$temp[6]\t$temp[7]";
	}
	close(FILE);
}

#import common name lod id
my (%wamei, %lod, %id, %lod2, %id2);
foreach(sort keys %file){
	my $file = $_;
	open(FILE, "<", ".\/Results\/4_1_Annotation\/${file}_Summary.txt") or die ("error:$!");
	while(<FILE>){
		my @temp = split(/\t/, $_);
		my $temp = @temp;
		if($temp > 3){
			if($dictionary){
				if($temp[0] =~ /^U$identity/){
					my $tempname = $temp[0];
					$tempname =~ s/_/ /g;
					$lod2{$file}{$tempname} = $temp[3];
					$id2{$file}{$tempname} = $temp[4];
				}else{
					$lod{$temp[0]}{$temp[3]}++;
					$id{$temp[0]}{$temp[4]}++;
				}
				$temp[0] =~ s/ otu\d+//;
				$wamei{$temp[0]} = $temp[1];
			}else{
				if($temp[0] =~ /^U$identity/){
					my $tempname = $temp[0];
					$tempname =~ s/_/ /g;
					$lod2{$file}{$tempname} = $temp[2];
					$id2{$file}{$tempname} = $temp[3];
				}else{
					$lod{$temp[0]}{$temp[2]}++;
					$id{$temp[0]}{$temp[3]}++;
				}
			}
		}
	}
	close(FILE);
}

#Output
open (OUT, ">", ".\/Results\/5_2_Summary_table\/Summary_table.tsv") or die("error:$!");
if($family){
	if($dictionary){
		if($habitat){print OUT "Family\tScientific Name\tCommon Name\tAve. Identity\tAve. LOD\tWater area\tHabitat\tDepthS\tDepthD";}
		else{print OUT "Family\tScientific Name\tCommon Name\tAve. Identity\tAve. LOD";}
	}else{
		if($habitat){print OUT "Family\tScientific Name\tAve. Identity\tAve. LOD\tWater area\tHabitat\tDepthS\tDepthD";}
		else{print OUT "Family\tScientific Name\tAve. Identity\tAve. LOD";}
	}
}else{
	if($dictionary){
		if($habitat){print OUT "Scientific Name\tCommon Name\tAve. Identity\tAve. LOD\tWater area\tHabitat\tDepthS\tDepthD";}
		else{print OUT "Scientific Name\tCommon Name\tAve. Identity\tAve. LOD";}
	}else{
		if($habitat){print OUT "Scientific Name\tAve. Identity\tAve. LOD\tWater area\tHabitat\tDepthS\tDepthD";}
		else{print OUT "Scientific Name\tAve. Identity\tAve. LOD";}
	}
}
foreach(sort keys %file){print OUT "\t$_";}
print OUT "\n";
my @output;
foreach(sort {$a cmp $b} keys %table){
	my $spname = $_;
	if($spname =~ /^U$identity/){
		$spname =~ /^U${identity}_([^\s]+)\s([^\s]+)/;
		my $genus = $1;
		my $cf = $2;
		if($genus =~ /^cf$/i){$genus = $cf;}
		
		#family
		my $temp;
		if($family){
			if($family{$genus}){$temp = "$family{$genus}\t$spname";}
			else{$temp = "NA\t$spname";}
		}else{$temp = "$spname";}
		
		my $tempname = $spname;
		$tempname =~ s/\sgotu\d+//;
		if($dictionary){
			if($wamei{$tempname}){$temp = $temp . "\t$wamei{$tempname}";}
			else{$temp = $temp . "\tNA";}
		}
		
		#identity
		if($seiri{$spname}){
			my $total = 0;
			my $score = 0;
			my $title = $seiri{$spname};
			foreach(keys %$title){
				$_ =~ /^>(.+)_Uniq\d+_\d+_reads_(.+)$/;
				my $name = $1; my $filename = $2;
				$name =~ s/_/ /g;
				if($id2{$filename}{$name}){
					$total++;
					$score = $score + $id2{$filename}{$name};
				}
			}
					
			$score = sprintf("%.2f", $score/$total);
			$temp = $temp . "\t$score";
		}else{$temp = $temp . "\tNA";}

		#LOD score
		if($seiri{$spname}){
			my $total = 0;
			my $score = 0;
			my $title = $seiri{$spname};
			foreach(keys %$title){
				$_ =~ /^>(.+)_Uniq\d+_\d+_reads_(.+)$/;
				my $name = $1; my $filename = $2;
				$name =~ s/_/ /g;
				if($lod2{$filename}{$name}){
					$total++;
					if($lod2{$filename}{$name} eq "HIGH"){$score = $score + 3;}
					elsif($lod2{$filename}{$name} eq "MODERATE"){$score = $score + 2;}
					elsif($lod2{$filename}{$name} eq "LOW"){$score = $score + 1;}
				}
			}
			
			$score = int($score/$total + 0.5);
			if($score == 3){$temp = $temp . "\tHIGH";}
			elsif($score == 2){$temp = $temp . "\tMODERATE";}
			elsif($score == 1){$temp = $temp . "\tLOW";}
		}else{$temp = $temp . "\tNA";}
		
		#habitat
		if($habitat){
			my $specific_name = "$genus $cf";
			if($habitat{$specific_name}){$temp = $temp . "\t$habitat{$specific_name}";}
			else{$temp = $temp . "\tNA\tNA\tNA\tNA";}
		}
		
		#reads
		foreach(sort keys %file){
			unless($table{$spname}{$_}){$temp = $temp . "\t0";}
			else{$temp = $temp . "\t$table{$spname}{$_}";}
		}
		
		push(@output, $temp);
	}else{
		$spname =~ /^([^\s]+)\s([^\s]+)/;
		my $genus = $1;
		my $cf = $2;
		if($genus =~ /^cf$/i){$genus = $cf;}
		
		my $temp;
		if($family){
			if($family{$genus}){$temp = "$family{$genus}\t$spname";}
			else{$temp = "NA\t$spname";}
		}else{$temp = "$spname";}
		
		if($dictionary){
			if($wamei{$spname}){$temp = $temp . "\t$wamei{$spname}";}
			else{$temp = $temp . "\tNA";}
		}
		
		if($id{$spname}){
			my $id = $id{$spname};
			my $total = 0;
			my $score = 0;
			foreach(keys %$id){$score = $score + $_ * $id{$spname}{$_}; $total += $id{$spname}{$_};}
			
			$score = sprintf("%.2f", $score/$total);
			$temp = $temp . "\t$score";
		}else{$temp = $temp . "\tNA";}
		
		if($lod{$spname}){
			my $lod = $lod{$spname};
			my $total = 0;
			my $score = 0;
			foreach(keys %$lod){
				if($_ eq "HIGH"){$score = $score + 3 * $lod{$spname}{$_}; $total += $lod{$spname}{$_};}
				elsif($_ eq "MODERATE"){$score = $score + 2 * $lod{$spname}{$_}; $total += $lod{$spname}{$_};}
				elsif($_ eq "LOW"){$score = $score + 1 * $lod{$spname}{$_}; $total += $lod{$spname}{$_};}
			}
			$score = int($score/$total + 0.5);
			if($score == 3){$temp = $temp . "\tHIGH";}
			elsif($score == 2){$temp = $temp . "\tMODERATE";}
			elsif($score == 1){$temp = $temp . "\tLOW";}
		}else{$temp = $temp . "\tNA";}
		
		#habitat
		if($habitat){
			my $specific_name = "$genus $cf";
			if($habitat{$specific_name}){$temp = $temp . "\t$habitat{$specific_name}";}
			else{$temp = $temp . "\tNA\tNA\tNA\tNA";}
		}
		
		foreach(sort keys %file){
			unless($table{$spname}{$_}){$temp = $temp . "\t0";}
			else{$temp = $temp . "\t$table{$spname}{$_}";}
		}
		push(@output, $temp);
	}
}
foreach(sort @output){print OUT "$_\n";}
my $spnumber = @output;
print "$spnumber species were detected\n";
close(OUT);
