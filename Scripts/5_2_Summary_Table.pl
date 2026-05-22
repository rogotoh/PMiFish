#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min);

#Setting
my ($primer, $separate, $identity, $identity2, $dictionary, $family, $habitat, $denoise, $type);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Primers\s*=\s*(\S+)/){my $temp = $1;if($temp =~ /^no$/i){$primer = "No";}else{$primer = $temp;}}
	elsif($_ =~ /^Divide\s*=\s*(\S+)/){$separate = $1;}
	elsif($_ =~ /^UIdentity\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /^LIdentity\s*=\s*(\S+)/){$identity2 = $1;}
	elsif($_ =~ /^Common_name\s*=\s*(\S+)/){$dictionary = $1;}
	elsif($_ =~ /^Family\s*=\s*(\S+)/){$family = $1;}
	elsif($_ =~ /^Habitat\s*=\s*(\S+)/){$habitat = $1;}
	elsif($_ =~ /^Algorithm\s*=\s*(\d+)/){$type = $1;}
	elsif($_ =~ /^Denoise\s*=\s*(\S+)/){$denoise = $1;}
}
close(SET);
unless($separate){$separate = 0;}
if($separate =~ /^yes$/i){$separate = 1;}
else{$separate = 0;}
unless($primer){print "Error: Please check the Primer file name in Setting.txt.\n"; exit;}
unless($identity){print "Error: Please check the upper threshold for homology search in Setting.txt.\n"; exit;}
unless($identity2){print "Error: Please check the lower threshold for homology search in Setting.txt.\n"; exit;}
unless(-f ".\/Dictionary\/$dictionary"){
	if($dictionary =~ /^no$/i){undef($dictionary);}
	else{print "Error: Dictionary file for common names not found or name mismatch in Setting.txt.\n"; exit;}
}
unless($primer =~ /^no$/i){
	unless(-f ".\/DataBase\/$primer"){print "Error: Primer file not found or primer file name mismatch in Setting.txt.\n"; exit;}
}
unless(-f ".\/Dictionary\/$family"){
	if($family =~ /^no$/i){undef($family);}
	else{print "Error: Dictionary file for family names not found or name mismatch in Setting.txt.\n"; exit;}
}
if($habitat){
	unless(-f ".\/Dictionary\/$habitat"){
		if($habitat =~ /^no$/i){undef($habitat);}
		else{print "Error: Dictionary file for habitat names not found or name mismatch in Setting.txt.\n"; exit;}
	}
}
unless($denoise){$denoise = "no";}
unless($type){$type = 1;}


my $as_us;
unless($denoise =~ /^yes$/i){
	$as_us = "Unique";
}else{
	if($type == 1){
		$as_us = "ZOTU";
	}else{
		$as_us = "ASV";
	}
}

#Get data names
opendir (DIR, ".\/Results\/4_1_Annotation") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_Representative_seq/){$file{$1}++;}
}
closedir DIR;

#Get primer sequences
my $countp = 0;
my (%forward, %reverse, $namep, @forward);
my ($primerF, $primerR);
my $primer_num;

unless($primer =~ /^no$/i){
	open (DATA2, "<", ".\/DataBase\/$primer") or die("error:$!");
	while(<DATA2>){
		if($_ =~ /\#/){next;}
		$_ =~ s/\r//g;
		chomp($_);
		if($_ =~ /Forward/){$countp = 1;next;}
		if($_ =~ /Reverse/){$countp = 2;next;}
		if($countp == 1){
			if($_ =~ /^>(.+)/){$forward{$1}++; push(@forward, $1);}
		}
		if($countp == 2){
			if($_ =~ /^>(.+)/){$reverse{$1}++;}
		}
	}
	close(DATA2);
	unless(%forward){
		print "Error: No primer sequence found in $primer. Please check the primer file.\n";
		exit;
	}
	unless(%reverse){
		print "Error: No primer sequence found in $primer. Please check the primer file.\n";
		exit;
	} 
	foreach(keys %forward){
		unless($reverse{$_}){
			print "Error: Paired primers must have the same name. Please check the primer file: $primer.\n";
			exit;
		}
	}
	foreach(keys %reverse){
		unless($forward{$_}){
			print "Error: Paired primers must have the same name. Please check the primer file: $primer.\n";
			exit;
		}
	}
	$primer_num = keys %forward;
	if($primer_num == 1){$separate = 0;}
}

print "============================================================\n";
print "                     5_2_Summary_Table                      \n";
print "============================================================\n";

#Summary_Table
mkdir ".\/Results\/5_2_Summary_Table";
my (%table, %under, $spname, $title, %phy, $title2, $spname2, %replist);
foreach(sort keys %file){
	my $file = $_;
	open (DATA, "<", ".\/Results\/4_1_Annotation\/${file}_Representative_seq.fas") or die("error:$!");
	while(<DATA>){
		my $temp_title;
		if($_ =~  /_${as_us}\d+/){
			$temp_title = "${as_us}";
		}else{
			$temp_title = "OTU";
		}
		if($_ =~ /^>Nohit/){next;}
		if($_ =~ /^>(.+)_$temp_title\d+_(\d+)_reads/){
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
		if($title2){
			chomp($_);
			$_ =~ s/\r//;
			unless($phy{$spname2}){$phy{$spname2} = $_;}
			undef($title2);undef($spname2);
		}
	}
	close(DATA);
}

#output Representative_seq.fas
open(OUT, ">", ".\/Results\/5_2_Summary_Table\/Representative_seq.fas") or die ("error:$!");
foreach(sort {$a cmp $b} keys %phy){
	my $name = $_;
	$name =~ s/ /_/g;
	print OUT ">$name\n$phy{$_}\n";
}
close(OUT);

#output Representative_list.txt
open(OUT, ">", ".\/Results\/5_2_Summary_Table\/Representative_list.txt") or die ("error:$!");
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
my (%wamei, %lod, %id, %lod2, %id2, %lca);
foreach(sort keys %file){
	my $file = $_;
	open(FILE, "<", ".\/Results\/4_1_Annotation\/${file}_Summary.txt") or die ("error:$!");
	while(<FILE>){
		chomp($_);
		$_ =~ s/\r//g;
		my @temp = split(/\t/, $_);
		my $temp = @temp;
		if($temp > 3){
			if($dictionary){
				$id{$temp[0]}{$temp[3]}++;
				$lod{$temp[0]}{$temp[4]}++;
				$lca{$temp[0]}{$temp[5]}++;
				#$temp[0] =~ s/ otu\d+//;
				$wamei{$temp[0]} = $temp[1];
			}else{
				$id{$temp[0]}{$temp[2]}++;
				$lod{$temp[0]}{$temp[3]}++;
				$lca{$temp[0]}{$temp[4]}++;
			}
		}
	}
	close(FILE);
}

#Output
open (OUT, ">", ".\/Results\/5_2_Summary_Table\/Summary_Table.tsv") or die("error:$!");
if($family){
	if($dictionary){
		if($habitat){print OUT "Family\tScientific Name\tCommon Name\tAve. Identity\tAve. LOD\tLCA Rank\tWater area\tHabitat\tDepthS\tDepthD";}
		else{print OUT "Family\tScientific Name\tCommon Name\tAve. Identity\tAve. LOD\tLCA Rank";}
	}else{
		if($habitat){print OUT "Family\tScientific Name\tAve. Identity\tAve. LOD\tLCA Rank\tWater area\tHabitat\tDepthS\tDepthD";}
		else{print OUT "Family\tScientific Name\tAve. Identity\tAve. LOD\tLCA Rank";}
	}
}else{
	if($dictionary){
		if($habitat){print OUT "Scientific Name\tCommon Name\tAve. Identity\tAve. LOD\tLCA Rank\tWater area\tHabitat\tDepthS\tDepthD";}
		else{print OUT "Scientific Name\tCommon Name\tAve. Identity\tAve. LOD\tLCA Rank";}
	}else{
		if($habitat){print OUT "Scientific Name\tAve. Identity\tAve. LOD\tLCA Rank\tWater area\tHabitat\tDepthS\tDepthD";}
		else{print OUT "Scientific Name\tAve. Identity\tAve. LOD\tLCA Rank";}
	}
}
foreach(sort keys %file){print OUT "\t$_";}
print OUT "\n";
my @output;
my $tsv = 0;
foreach(sort {$a cmp $b} keys %table){
	my $spname = $_;
	if($spname =~ /k\s\s/){$tsv = 1;}
	my $temp_title;
	if($_ =~  /_${as_us}\d+/){
		$temp_title = "${as_us}";
	}else{
		$temp_title = "OTU";
	}
	my ($genus, $cf);
	if($spname =~ /U$identity/){
		$spname =~ /^U$identity\s([^\s]+)\s([^\s]+)/;
		$genus = $1;
		$cf = $2;
	}else{
		$spname =~ /^([^\s]+)\s([^\s]+)/;
		$genus = $1;
		$cf = $2;
	}
	unless($genus){$genus = "na";}
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
	
	if($lca{$spname}){
		my $lca = $lca{$spname};
		my $total = 0;
		my @rank = qw(Kingdom >Phylum Phylum >Class Class >Order Order >Family Family >Genus Genus Species);
		foreach(@rank){
			if($lca{$spname}{$_}){
				$temp = $temp . "\t$_";
				last;
			}
		}
	}
	
	
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
foreach(sort @output){print OUT "$_\n";}
my $spnumber = @output;

#divide table
if($separate){
	my $colum = 0;
	foreach(sort @forward){
		my $temp_p = $_;
		open (OUT, ">", ".\/Results\/5_2_Summary_Table\/Summary_Table_${_}.tsv") or die("error:$!");
		open (DATA, "<", ".\/Results\/5_2_Summary_Table\/Summary_Table.tsv") or die("error:$!");
		my $count = 0;
		my %check_c;
		my %fasta;
		while(<DATA>){
			chomp($_);
			$_ =~ s/\r//g;
			my @temp = split(/\t/, $_);
			my $narabi = "$temp[0]";
			if($count == 0){
				for(my $n = 1; $n < @temp; $n++){
					if($temp[$n] =~ /$temp_p/){
						unless($colum){$colum = $n;}
						$check_c{$n}++;
					}
				}
			}
			for(my $n = 1; $n < $colum; $n++){
				$narabi .= "\t$temp[$n]";
			}
			my $sample;
			my $total_read = 0;
			for(my $n = $colum; $n < @temp; $n++){
				if($check_c{$n}){
					if($count == 0){
						unless($sample){$sample = "\t$temp[$n]";}
						else{$sample .= "\t$temp[$n]";}
					}else{
						$total_read += $temp[$n];
						unless($sample){$sample = "\t$temp[$n]";}
						else{$sample .= "\t$temp[$n]";}
					}
				}
			}
			if($count == 0){
				if($sample){print OUT "$narabi$sample\n";}
			}else{
				if($total_read){
					print OUT "$narabi$sample\n";
					my @temp = split(/\t/, $narabi);
					$fasta{$temp[1]}++;
				}
			}
			$count = 1;
		}
		close(OUT);
		close(DATA);
		
		#divide fasta
		open (OUT, ">", ".\/Results\/5_2_Summary_Table\/Representative_seq_${temp_p}.fas") or die("error:$!");
		open (DATA, "<", ".\/Results\/5_2_Summary_Table\/Representative_seq.fas") or die("error:$!");
		my $seq = 0;
		while(<DATA>){
			chomp($_);
			$_ =~ s/\r//g;
			if($_ =~ /^>/){
				$_ =~ s/_from(.+)//;
				$_ =~ s/_/ /g;
				$_ =~ /^>(.+)/;
				if($fasta{$1}){
					$_ =~ s/ /_/g;
					print OUT "$_\n";
					$seq = 1;
				}
			}else{
				if($seq){
					print OUT "$_\n";
					$seq = 0;
				}
			}
		}
		close(OUT);
		close(DATA);
	}
}

#Rank tsv
if($tsv == 1){
	opendir (DIR, ".\/Results\/5_2_Summary_Table") or die ("error:$!");
	my @database = readdir DIR;
	my @tsv;
	foreach(@database){
		chomp($_);
		$_ =~ s/\r//g;
		if ($_ =~ /\.tsv/){
			unless($_ =~ /Rank\.tsv/){push(@tsv, $_);}
		}
	}
	closedir DIR;
	foreach(@tsv){
		my $filename = $_;
		$filename =~ /(.+)\.tsv/;
		my $fname = $1;
		open (DATA, "<", ".\/Results\/5_2_Summary_Table\/$filename") or die("error:$!");
		open (OUT, ">", ".\/Results\/5_2_Summary_Table\/${fname}_Rank.tsv") or die("error:$!");
		my $count = 0;
		my @out;
		while(<DATA>){
			if($count == 0){
				$_ =~ s/Scientific Name/Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies/;
				$count = 1;
				print OUT "$_";
			}else{
				my @temp = split(/\t/, $_);
				$temp[0] =~ s/[kpcofgs]\s\s//g;
				$temp[0] =~ s/\;/\t/g;
				$temp[0] =~ s/(U\d+)\s//;
				my $under = $1;
				if($under){
					my @temp2 = split(/\t/, $temp[0]);
					$temp2[-1] = "$under $temp2[-1]";
					$temp[0] = join("\t", @temp2);
				}
				$_ = join("\t", @temp);
				push(@out, $_);
			}
		}
		foreach(sort {$a cmp $b} @out){
			print OUT "$_";
		}
		close(DATA);
		close(OUT);
	}
}

print "$spnumber species were detected\n";
close(OUT);


