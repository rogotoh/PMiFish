#!/usr/bin/perl
use strict;
use warnings;

#2018/01/05 add coverage filter (not less than 90%)

#Setting
my ($db, $depth, $size, $identity, $rarefy, $identity2, $dictionary, $family, $temporary);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^DB\s*=\s*(\S+)/){$db = $1;}
	elsif($_ =~ /^Rarefaction\s*=\s*(\S+)/){$rarefy = $1;}
	elsif($_ =~ /^Depth\s*=\s*(\S+)/){$depth = $1;}
	elsif($_ =~ /^Length\s*=\s*(\d+)/){$size = $1;}
	elsif($_ =~ /^UIdentity\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /^LIdentity\s*=\s*(\S+)/){$identity2 = $1;}
	elsif($_ =~ /^Dictionary\s*=\s*(\S+)/){$dictionary = $1;}
	elsif($_ =~ /^Family\s*=\s*(\S+)/){$family = $1;}
	elsif($_ =~ /^Temporary\s*=\s*(\S+)/){$temporary = $1;}
}
close(SET);
if($rarefy =~ /^no$/i){undef($rarefy);}
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
opendir (DIR, ".\/Results\/3_1_Usearch_global") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_usearch_results.txt/){$file{$1}++;}
}
closedir DIR;

print "============================================================\n";
print "                      4_1_Annotation                       \n";
print "============================================================\n";

#Annotation
mkdir ".\/Results\/4_1_Annotation";
foreach(sort keys %file){
	print "$_ Annotate...\n";
	my $file = $_;
	my (%jp, %genus);
	if($dictionary){
		open (DATAFILE, "<", ".\/Dictionary\/$dictionary") or die("error:$!");
		while(<DATAFILE>){
			chomp($_);
			$_ =~ s/\r//g;
			my @temp = split(/\t/, $_);
			if($temp[1]){$jp{$temp[0]} = $temp[1];}
			if($temp[1] and $temp[0] and $temp[0] =~ /^([^ ]+) /){
				unless($genus{$1}){
					my $kensaku = $1;
					if($dictionary =~ /Sname_Jname/){$genus{$kensaku} = "${temp[1]}Ç∆ìØÇ∂ëÆÇÃéÌ";}
					else{$genus{$kensaku} = "The species belongs to the same genus as ${temp[1]}";}
				}
			}
		}
		close(DATAFILE);
	}
	
	#import zotu_nonchimeras.fas data 
	my ($data, $dseq, %lead);
	if(-f ".\/Results\/2_4_Rarefaction\/${file}_rarefy.fa"){
		open (DATA, "<", ".\/Results\/2_4_Rarefaction\/${file}_rarefy.fa") or die("error:$!");
	}elsif(-f ".\/Results\/2_3_Separate_chimera\/${file}_zotu_nonchimeras.fa"){
		open (DATA, "<", ".\/Results\/2_3_Separate_chimera\/${file}_zotu_nonchimeras.fa") or die("error:$!");
	}else{
		open (DATA, "<", ".\/Results\/2_1_Find_unique\/${file}_uniques.fa") or die("error:$!");
	}
	while(<DATA>){
		chomp($_);
		$_ =~ s/\r//g;
		unless($_ =~ /^\n/){
			if($_ =~ /^>(.+\;)/){
				if($data){$lead{$data} = $dseq; $data = $1; undef($dseq);}
				else{$data = $1;}
			}else{
				if($dseq){$dseq = $dseq . $_;}
				else{$dseq = $_;}
			}
		}
	}
	if($data){$lead{$data} = $dseq;}
	close(DATA);
	
	#usearch_results.txt
	open (DATA, "<", ".\/Results\/3_1_Usearch_global\/${file}_usearch_results.txt") or die("error:$!");
	my (%twohit, @blastdata, %num);
	while(<DATA>){
		chomp($_);
		$_ =~ s/\r//g;
		$_ =~ s/\t+/\t/g;
		my @temp = split(/\t/, $_);
		$temp[1] =~ s/.+\|//;
		if($temp[-1] < 90){next;} #query cover %
		unless($twohit{$temp[0]}){$twohit{$temp[0]} = $temp[1]; push(@blastdata, $_); $num{$temp[0]}++;}
		else{
			if($num{$temp[0]} == 1 and $twohit{$temp[0]} ne $temp[1]){push(@blastdata, $_); $num{$temp[0]}++;}
		}
	}
	close(DATA);
	
	my (%kekka, %check, %count, %kekka2, %count2, $name, %up, %down, %lowhit);
	foreach(@blastdata){
		chomp($_);
		my @temp = split(/\t/, $_);
		unless($temp[0] =~ /\;$/){$temp[0] = $temp[0] . "\;";}
		$check{$temp[0]}++;
		if($check{$temp[0]} == 1 and $temp[2] >= $identity){
			$temp[1] =~ /\|(.+)\|/;
			my $gbn = $1;
			$temp[1] =~ s/.+\|//;
			$name = $temp[1];
			my $ketu = "$temp[2]\t$temp[3]\t$temp[4]\t$gbn";
			$kekka{$temp[1]}{$temp[0]} = $ketu;
			$temp[0] =~ /size=(\d+)/;
			if($count{$temp[1]}){$count{$temp[1]} += $1;}
			else{$count{$temp[1]} = $1;}
			$up{$temp[0]}++;
		}elsif($check{$temp[0]} == 1 and $temp[2] < $identity){
			$temp[1] =~ /\|(.+)\|/;
			my $gbn = $1;
			$temp[1] =~ s/.+\|//;
			$name = $temp[1];
			my $ketu = "$temp[2]\t$temp[3]\t$temp[4]\t$gbn";
			$kekka2{$temp[1]}{$temp[0]} = $ketu;
			$temp[0] =~ /size=(\d+)/;
			$lowhit{$temp[1]}{$temp[0]} = $lead{$temp[0]};
			$down{$temp[0]}++;
		}elsif($check{$temp[0]} == 2){
			$temp[1] =~ /\|(.+)\|/;
			my $gbn = $1;
			$temp[1] =~ s/.+\|//;
			$temp[1] =~ s/_/ /g;
			my $ketu = "$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$gbn";
			if($up{$temp[0]}){
				$kekka{$name}{$temp[0]} = $kekka{$name}{$temp[0]} . "\t$ketu";}
			elsif($down{$temp[0]}){$kekka2{$name}{$temp[0]} = $kekka2{$name}{$temp[0]} . "\t$ketu";}
		}
	}

	#clusterization under $identity seqs
	my %seiri;
	foreach (sort {$a cmp $b} keys %lowhit){
		my $spname = $_;
		my $title = $lowhit{$spname};
		my $seqnum = keys %$title;
		if($seqnum > 1){
			my @keys = keys %$title;
			my (%narabi, %tempname, %modosi);
			my $tempcount = 0;
			foreach(@keys){$_ =~ /size=(\d+)/; $narabi{$1}{$_}++; $tempcount++; $tempname{$_} = "Uniq$tempcount;size=$1;"; $modosi{"Uniq$tempcount;size=$1;"} = $_;}
			undef(@keys);
			open (TEMP, ">", ".\/Results\/4_1_Annotation\/temp.fas") or die("error:$!");
			foreach(sort {$b <=> $a} keys %narabi){
				my $suji = $narabi{$_};
				foreach(sort keys %$suji){print TEMP ">$tempname{$_}\n$lowhit{$spname}{$_}\n";}
			}
			close(TEMP);
			
			my $tempid = $identity/100;
			my $comand = ".\/Tools\/$usearch -cluster_smallmem \".\/Results\/4_1_Annotation\/temp.fas\" -id $tempid -uc \".\/Results\/4_1_Annotation\/temp.uc\" -sortedby size -quiet";
			system $comand;
			
			open (TEMP, "<", ".\/Results\/4_1_Annotation\/temp.uc") or die("error:$!");
			my $sp = 0;
			my %otus;
			while(<TEMP>){
				chomp($_);
				$_ =~ s/\r//g;
				my @temp = split(/\t/, $_);
				if($temp[0] eq "S"){$sp++; $seiri{"${spname}_otu$sp"}{$modosi{$temp[8]}}++; $otus{$temp[8]} = "${spname}_otu$sp";}
				elsif($temp[0] eq "H"){$seiri{$otus{$temp[9]}}{$modosi{$temp[8]}}++;}
			}
			close(TEMP);
			unlink ".\/Results\/4_1_Annotation\/temp.fas";
			unlink ".\/Results\/4_1_Annotation\/temp.uc";
		}else{
			my @keys = keys %$title;
			$seiri{$spname}{$keys[0]}++;
		}
	}

	#read counts
	foreach(keys %seiri){
		my $spname = $_;
		my $title = $seiri{$spname};
		foreach(keys %$title){
			$_ =~ /size=(\d+)/;
			if($count2{$spname}){$count2{$spname} += $1;}
			else{$count2{$spname} = $1;}
		}
	}
	
	my (%reverse, %reverse2, %nohit);
	foreach(keys %count){$reverse{$count{$_}}{$_}++;}
	foreach(keys %count2){$reverse2{$count2{$_}}{$_}++;}
	my $species = 0;
	my $species2 = 0;
	$species = keys %count;
	$species2 = keys %count2;

	#DB
	my ($ref, %ref);
	open (FAS, "<", ".\/DataBase\/$db") or die("error:$!");
	while(<FAS>){
		chomp($_);
		$_ =~ s/\r//g;
		$_ =~ s/\t//g;
		if($_ =~ /^>/){$_ =~ s/.+\|//; $ref = $_;}
		else{$ref{$_}{$ref}++;}
	}
	close(FAS);

	open (OUT1, ">", ".\/Results\/4_1_Annotation\/${file}_Summary.txt") or die("error:$!");
	open (OUT2, ">", ".\/Results\/4_1_Annotation\/${file}_Detail.txt") or die("error:$!");
	open (OUT3, ">", ".\/Results\/4_1_Annotation\/${file}_Representative_seq.fas") or die("error:$!");
	open (OUT4, ">", ".\/Results\/4_1_Annotation\/${file}_Synonym_list.txt") or die("error:$!");
	open (OUT5, ">", ".\/Results\/4_1_Annotation\/${file}_all_annotated_seq.fas") or die("error:$!");

	my %synonym;
	print OUT1 "===Results of $file===\n\n";
	if(%reverse){
		if($dictionary){print OUT1 "Identity >= ${identity}%\t$species species\nScientific_name\tCommon_name\tReads\tConfidence\tIdentity(%)\tSynonym\n";}
		else{print OUT1 "Identity >= ${identity}%\t$species species\nScientific_name\tReads\tSynonym\n";}
		print OUT2 "Identity >= ${identity}%\t$species species\n";
		foreach (sort {$b <=> $a} keys %reverse){
			my $temp = $_;						#Total number of reads
			my $kame = $reverse{$temp};
			foreach(sort {$a cmp $b} keys %$kame){
				my $temp2 = $_;					#scientific name
				my $cut_cluster = $_;
				$cut_cluster =~ s/_cluster\d+//;
				$cut_cluster =~ s/_/ /g;
				$_ =~ s/_/ /g;
				my $wamei;
				#òañºÇ÷ÇÃïœä∑
				if($dictionary){
					if($jp{$cut_cluster}){$wamei = $jp{$cut_cluster};}
					else{
						$_ =~ /^([^ ]+) /;
						if($genus{$1}){$wamei = $genus{$1};}
						else{
							if($dictionary =~ /Sname_Jname/){$wamei = "äYìñëÆñºÅEòañºÇ»Çµ";}
							else{$wamei = "No applicable name";}
						}
					}
					print OUT1 "$_\t$wamei\t$temp\t";
				}else{print OUT1 "$_\t$temp\t";}
				
				print OUT2 "$_\t$temp reads\n";
				my $kame2 = $kekka{$temp2};
				my %temp;
				foreach(keys %$kame2){
					$_ =~ /Uniq(\d+)\;/;
					$temp{$1} = $_;
				}
				print OUT2 "\tReads\tConfidence\tIdentity(%)\tLOD_score\tAlign_len\tMismatch\tAccession No.\t2nd-sp_name\t2nd_Identity(%)\t2nd_Align_len\t2nd_Mismatch\tAccession No.\tSequence\n";
				my $fasta = 0;
				foreach(sort {$a <=> $b} keys %temp){
					$fasta++;
					my $read = $temp{$_};
					$read =~ /(Uniq\d+)\;size=(\d+)\;/;
					my $id = $1;
					my $cut = $2;
					my $kazu = 0;
					my $kiri = $kekka{$temp2}{$read};
					while($kiri =~ s/\t//i){$kazu++;}
					
					#Make Synonym, Representative, all_annotated_seq
					print OUT5 ">${temp2}_${id}_${cut}_reads\n$lead{$read}\n";
					if($fasta == 1){
						if($ref{$lead{$read}}){
							my $syno = $ref{$lead{$read}};
							my $count =0;
							foreach(sort keys %$syno){
								unless($temp2 eq $_){
									$count++;
									my $temp3;
									if($count == 1){
										$temp3 = $temp2;
										my $cut_cluster = $temp2;
										$cut_cluster =~ s/_cluster\d+//;
										$cut_cluster =~ s/_/ /g;
										$temp3 =~ s/_/ /g;
										my $wamei;
										if($dictionary){
											if($jp{$cut_cluster}){$wamei = $jp{$cut_cluster};}
											else{
												$temp3 =~ /^([^ ]+) /;
												if($genus{$1}){$wamei = $genus{$1};}
												else{
													if($dictionary =~ /Sname_Jname/){$wamei = "äYìñëÆñºÅEòañºÇ»Çµ";}
													else{$wamei = "No applicable name";}
												}
											}
											print OUT4 ">$temp3\t$wamei\n";
											$synonym{$temp3}++;
										}else{
											if($count == 1){print OUT4 ">$temp3\n"; $synonym{$temp3}++;}
										}
									}
									$temp3 = $_;
									my $cut_cluster = $_;
									$cut_cluster =~ s/_cluster\d+//;
									$cut_cluster =~ s/_/ /g;
									$temp3 =~ s/_/ /g;
									my $wamei;
									if($dictionary){
										if($jp{$cut_cluster}){$wamei = $jp{$cut_cluster};}
										else{
											$temp3 =~ /^([^ ]+) /;
											if($genus{$1}){$wamei = $genus{$1};}
											else{
												if($dictionary =~ /Sname_Jname/){$wamei = "äYìñëÆñºÅEòañºÇ»Çµ";}
												else{$wamei = "No applicable name";}
											}
										}
										print OUT4 "\t$temp3\t$wamei\n";
									}else{
										print OUT4 "\t$temp3\n";
									}
								}
							}
						}
						print OUT3 ">${temp2}_${id}_${temp}_reads\n$lead{$read}\n";
					}
					
					#Make Detail file
					if($kazu > 6){
						my @temp = split(/\t/, $kekka{$temp2}{$read});
						my $lod = log((($temp[1])/($temp[2]+1))/(($temp[6])/($temp[7]+1)));
						$lod = sprintf("%.4f", $lod);
						my $conf;
						if($lod >= 0.9){$conf = "HIGH";}
						elsif($lod >= 0.5){$conf = "MODERATE";}
						else{$conf = "LOW";}
						if($fasta == 1){print OUT1 "$conf\t$temp[0]\n";}
						print OUT2 "\t$cut\t$conf\t$temp[0]\t$lod\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\t$lead{$read}\n";}
					else{
						my @temp = split(/\t/, $kekka{$temp2}{$read});
						my $lod = "N.A.";
						if($fasta == 1){print OUT1 "HIGH\t$temp[0]\n";}
						print OUT2 "\t$cut\tHIGH\t$temp[0]\t$lod\t$temp[1]\t$temp[2]\t$temp[3]\tN.A.\tN.A.\tN.A.\tN.A.\tN.A.\t$lead{$read}\n";
					}
					$nohit{$read}++;
				}
		 	}
		}
	}else{
		print OUT1 "\nIdentity >= ${identity}%\t0 species\nNo applicable sequence\n";
		print OUT2 "\nIdentity >= ${identity}%\t0 species\nNo applicable sequence\n";
	}

	#lower identity results
	if(%reverse2){
		if($dictionary){print OUT1 "\n${identity2}% < Identity < ${identity}%\t$species2 species\nScientific_name\tCommon_name\tReads\tConfidence\tIdentity(%)\n";}
		else{print OUT1 "\n${identity2}% < Identity < ${identity}%\t$species2 species\nScientific_name\tReads\tReads\tConfidence\tIdentity(%)\n";}
		print OUT2 "\n${identity2}% < Identity < ${identity}%\t$species2 species\n";
		foreach (sort {$b <=> $a} keys %reverse2){
			my $temp = $_;						#Total number of reads
			my $kame = $reverse2{$temp};
			foreach(sort {$a cmp $b} keys %$kame){
				my $temp2 = $_;					#scientific name
				my $cut_cluster = $_;
				$cut_cluster =~ s/_cluster\d+//;
				$cut_cluster =~ s/_/ /g;
				$_ =~ s/_/ /g;
				my $wamei;
				#òañºÇ÷ÇÃïœä∑
				if($dictionary){
					if($jp{$cut_cluster}){
						if($dictionary =~ /Sname_Jname/){$wamei = "$jp{$cut_cluster}Ç∆ãﬂâèÇ»éÌ";}
						else{$wamei = "The species closely rerated to $jp{$cut_cluster}";}
					}else{
						$_ =~ /^([^ ]+) /;
						if($genus{$1}){
							my $itiji = $1;
							if($dictionary =~ /Sname_Jname/){$wamei = "$genus{$itiji}Ç∆ãﬂâèÇ»éÌ";}
							else{$wamei = "The species closely rerated to $genus{$itiji}";}
						}else{
							if($dictionary =~ /Sname_Jname/){$wamei = "äYìñëÆñºÅEòañºÇ»Çµ";}
							else{$wamei = "No applicable name";}
						}
					}
					print OUT1 "U${identity}_$_\t$wamei\t$temp\t";
				}else{print OUT1 "U${identity}_$_\t$temp\t";}
				
				print OUT2 "U${identity}_$_\t$temp reads\n";
				my $kame2 = $seiri{$temp2};
				my %temp;
				foreach(keys %$kame2){
					$_ =~ /Uniq(\d+)\;/;
					$temp{$1} = $_;
				}
				print OUT2 "\tReads\tConfidence\tIdentity(%)\tLOD_score\tAlign_len\tMismatch\tAccession No.\t2nd-sp_name\t2nd_Identity(%)\t2nd_Align_len\t2nd_Mismatch\tAccession No.\tSequence\n";
				my $fasta = 0;
				foreach(sort {$a <=> $b} keys %temp){
					$fasta++;
					my $read = $temp{$_};
					$read =~ /(Uniq\d+)\;size=(\d+)\;/;
					my $id = $1;
					my $cut = $2;
					my $kazu = 0;
					my $temp3 = $temp2;
					$temp3 =~ s/_otu\d+//;
					my $kiri = $kekka2{$temp3}{$read};
					while($kiri =~ s/\t//i){$kazu++;}
					
					#Make Representative, all_annotated_seq
					print OUT5 ">U${identity}_${temp2}_${id}_${cut}_reads\n$lead{$read}\n";
					if($fasta == 1){print OUT3 ">U${identity}_${temp2}_${id}_${temp}_reads\n$lead{$read}\n";}
					
					#Make Detail file
					if($kazu > 6){
						my @temp = split(/\t/, $kekka2{$temp3}{$read});
						my $lod = log((($temp[1])/($temp[2]+1))/(($temp[6])/($temp[7]+1)));
						$lod = sprintf("%.4f", $lod);
						my $conf;
						if($lod >= 0.9){$conf = "HIGH";}
						elsif($lod >= 0.5){$conf = "MODERATE";}
						else{$conf = "LOW";}
						if($fasta == 1){print OUT1 "$conf\t$temp[0]\n";}
						print OUT2 "\t$cut\t$conf\t$temp[0]\t$lod\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\t$lead{$read}\n";}
					else{
						my @temp = split(/\t/, $kekka2{$temp3}{$read});
						my $lod = "N.A.";
						if($fasta == 1){print OUT1 "HIGH\t$temp[0]\n";}
						print OUT2 "\t$cut\tHIGH\t$temp[0]\t$lod\t$temp[1]\t$temp[2]\t$temp[3]\tN.A.\tN.A.\tN.A.\tN.A.\tN.A.\t$lead{$read}\n";
					}
					$nohit{$read}++;
				}
		 	}
		}
	}else{
		print OUT1 "\n${identity2}% < Identity < ${identity}%\t0 species\nNo applicable sequence\n";
		print OUT2 "\n${identity2}% < Identity < ${identity}%\t0 species\nNo applicable sequence\n";
	}
	

	#Output nohit sequence
	my %temp;
	my @temp = keys %lead;
	foreach(@temp){
		$_ =~ /Uniq(\d+)\;/;
		$temp{$1} = $_;
	}
	my $nohit_n = 0;
	foreach(sort {$a <=> $b} keys %temp){
		unless($nohit{$temp{$_}}){$nohit_n++;}
	}
	
	print OUT1 "\nNohit\t$nohit_n sequences\n";
	print OUT2 "\nNohit\t$nohit_n sequences\n";

	foreach(sort {$a <=> $b} keys %temp){
		my $read = $temp{$_};
		unless($nohit{$read}){
			print OUT1 ">$read\n$lead{$read}\n";
			print OUT2 ">$read\n$lead{$read}\n";
			$read =~ /(Uniq\d+)\;size=(\d+)\;/;
			print OUT3 ">Nohit_${1}_${2}_reads\n$lead{$read}\n";
			print OUT5 ">Nohit_${1}_${2}_reads\n$lead{$read}\n";
		}
	}
	close(OUT1);
	close(OUT2);
	close(OUT3);
	close(OUT4);
	close(OUT5);

	#add synonym data to summary
	open (DATA, "<", ".\/Results\/4_1_Annotation\/${file}_Summary.txt") or die("error:$!");
	my @result;
	while(<DATA>){
		chomp($_);
		$_ =~ s/\r//g;
		my $check;
		if($_ =~ /([^\t]+)\t/){$check = $1;}
		if($check and $synonym{$check}){$_ = $_ . "\tYES";}
		push(@result, $_);
	}
	close(DATA);
	
	open (DATA, ">", ".\/Results\/4_1_Annotation\/${file}_Summary.txt") or die("error:$!");
	foreach(@result){print DATA "$_\n";}
	close(DATA);
	
	#Make Detail html
	open (OUT, ">", ".\/Results\/4_1_Annotation\/${file}_Detail.html") or die("error:$!");
	print OUT <<"EOS";
<html>
<head>
<meta http-equiv="Content-type" content="text/html" charset="Shift_JIS">
<title>${file}_Detail</title>
<style type="text/css">
H1{color: #ffffff; text-align: left; font: 100% Tahoma;}
H2{color: #000000; text-align: left; font: 100% Tahoma;}
H3{text-align: left; font: 100% Tahoma; line-height: 4px;}
</style>
</head>
<body bgcolor="#EFEFFB">
<a name="top"></a>
<font face="Tahoma" size="6">Detailed Result of $file</font><BR>
EOS
	open (DATA, "<", ".\/Results\/4_1_Annotation\/${file}_Detail.txt") or die("error:$!");
	my $iro = 0;
	my $table = 1;
	my $nohit = 0;
	my ($name2, %portal);
	while(<DATA>){
		if($_ =~ /^\n|^\r/){next;}
		chomp($_);
		$_ =~ s/\r//g;
		if($_ =~ /^Identity/){print OUT "<font face=\"Tahoma\" size=\"3\"><B>$_</B></font><BR>\n";next;}
		elsif($_=~ /^${identity2}%/){print OUT "</table>\n<br clear = \"all\"><br>\n<font face=\"Tahoma\" size=\"3\"><B>$_</B></font><BR>\n";$iro = 0; $table = 1; next;}
		elsif($_ =~ /^Nohit/){print OUT "</table>\n<br clear = \"all\"><br>\n<font face=\"Tahoma\" size=\"3\"><B>$_</B></font><BR>\n"; $nohit = 1; next;}
		unless($_ =~ s/^\t//){
			if($nohit){print OUT "<font face=\"Tahoma\" size=\"3\">$_</font><BR>\n";next;}
			$_ =~ /([^\t]+)\t/;
			if($1){$name2 = $1;}
		}else{
			my @temp = split(/\t/, $_);
			if($_ =~ /^Reads/){
				if($table == 1){
					print OUT "<table><table border=\"0\" cellspacing=\"1\" bgcolor=\"\#191970\" border=\"1\" align=\"left\">\n<tr>\n";
					print OUT "<th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>Scientific name&nbsp;&nbsp;</H1></th>\n";
					foreach(@temp){print OUT "<th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>$_&nbsp;&nbsp;</H1></th>\n";}
					print OUT "</tr>\n";
					$table++;
				}else{next;}
			}else{$iro++; $portal{$name2}++; &nakami3($iro, \@temp, $name2, \%portal);}
		}
	}
	close(DATA);
	close(OUT);
	unlink ".\/Results\/4_1_Annotation\/${file}_Detail.txt";
	
	#Make Synonym_list html
	open (OUT, ">", ".\/Results\/4_1_Annotation\/${file}_Synonym_list.html") or die("error:$!");
	print OUT <<"EOS";
<html>
<head>
<meta http-equiv="Content-type" content="text/html" charset="Shift_JIS">
<title>${file}_Synonym_list</title>
<style type="text/css">
H1{color: #ffffff; text-align: left; font: 100% Tahoma;}
H2{color: #000000; text-align: left; font: 100% Tahoma;}
H3{text-align: left; font: 100% Tahoma; line-height: 4px;}
</style>
</head>
<body bgcolor="#EFEFFB">
<a name="top"></a>
<font face="Tahoma" size="6">Synonym list of $file</font><BR>
<font face="Tahoma" size="4">(Species list having the same sequence in DB)</font><BR><br>
EOS
	open (DATA, "<", ".\/Results\/4_1_Annotation\/${file}_Synonym_list.txt") or die("error:$!");
	$table = 0;
	while(<DATA>){
		chomp($_);
		$_ =~ s/\r//g;
		if($_ =~ s/^>//){
			$table++;
			if($table > 1){print OUT "</table>\n<br clear = \"all\"><br>\n";}
			$iro = 0;
			my @temp = split(/\t/, $_);
			if($dictionary){print OUT "<font face=\"Tahoma\" size=\"4\"><B><I>$temp[0]</I>&nbsp;&nbsp;$temp[1]</B></font><BR>\n";}
			else{print OUT "<font face=\"Tahoma\" size=\"4\"><B><I>$temp[0]</I></B></font><BR>\n";}
			print OUT "<table><table border=\"0\" cellspacing=\"1\" bgcolor=\"\#191970\" border=\"1\" align=\"left\">\n<tr>\n";
			if($dictionary){print OUT "<th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>Scientific name </H1></th><th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>Common name </H1></th>\n</tr>\n";}
			else{print OUT "<th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>Scientific name </H1></th>\n</tr>\n";}
		}else{
			$_ =~ s/^\t//;
			my @temp = split(/\t/, $_);
			print OUT "<tr>\n";
			if($dictionary){print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2><I>$temp[0]</I>&nbsp;&nbsp;</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[1]&nbsp;&nbsp;</H2></td>\n";}
			else{print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2><I>$temp[0]</I>&nbsp;&nbsp;</H2></td>\n";}
			print OUT "</tr>\n";
		}
	}
	print OUT "</table>\n<br clear = \"all\"><br>\n";
	close(DATA);
	close(OUT);
	unlink ".\/Results\/4_1_Annotation\/${file}_Synonym_list.txt";
}

sub nakami3{
	my ($iro, $temp, $name, $portal) = @_;
	my @temp = @$temp;
	my %portal = %$portal;
	my $synonym;
	if($iro%2 + 1 == 2){
		print OUT "<tr>\n";
		if($portal{$name} and $portal{$name} == 1){print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2><I>$name</I>&nbsp;&nbsp;</H2></td>\n";}
		else{print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2></H2></td>\n";}
		for(my $n = 0; $n < @temp; $n++){
			if($n == 6 or $n == 11){
				unless($temp[$n] =~ /N.A./){print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><a href=\"https://www.ncbi.nlm.nih.gov/nuccore/$temp[$n]\" target=Åh_blankÅh><H2>$temp[$n]</H2></a></td>\n";}
				else{print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[$n]&nbsp;&nbsp;</H2></td>\n";}
			}elsif($n == 7){print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2><I>$temp[$n]</I>&nbsp;&nbsp;</H2></td>\n";}
			else{print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[$n]&nbsp;&nbsp;</H2></td>\n";}
		}
		print OUT "</tr>\n";
	}else{
		print OUT "<tr>\n";
		if($portal{$name} and $portal{$name} == 1){print OUT "<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2><I>$name</I>&nbsp;&nbsp;</H2></td>\n";}
		else{print OUT "<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2></H2></td>\n";}
		for(my $n = 0; $n < @temp; $n++){
			if($n == 6 or $n == 11){
				unless($temp[$n] =~ /N.A./){print OUT "<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><a href=\"https://www.ncbi.nlm.nih.gov/nuccore/$temp[$n]\" target=Åh_blankÅh><H2>$temp[$n]</H2></a></td>\n";}
				else{print OUT "<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[$n]&nbsp;&nbsp;</H2></td>\n";}
			}elsif($n == 7){print OUT "<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2><I>$temp[$n]</I>&nbsp;&nbsp;</H2></td>\n";}
			else{print OUT "<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[$n]&nbsp;&nbsp;</H2></td>\n";}
		}
		print OUT "</tr>\n";
	}
}
