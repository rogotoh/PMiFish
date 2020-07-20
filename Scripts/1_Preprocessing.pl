#!/usr/bin/perl
use strict;
use warnings;

#2018/12/04 1_2 Strip_primers change to use the Primer_Cleaner.pl
#2018/07/01 1_2_Strip_primers.pl were changed to permit no primer seq
#2018/01/29 change Blast to Usearch_global
#2018/01/20 output a log file at 1_1 step
#2018/01/05 add compress option

#This script execute below scripts in each sample
#1_1 Merge_paird_reads.pl
#1_2_Strip_primers.pl
#1_3_Quality_filter.pl


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

#Get primer sequences
my $countp = 0;
my (%forward, %reverse, $namep);
my ($primerF, $primerR);

unless($primer =~ /^no$/i){
	open (DATA2, "<", ".\/DataBase\/$primer") or die("error:$!");
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
	foreach(keys %forward){
		unless($reverse{$_}){print "Error!! Paird primer need to be the same name!\nCheck $primer\n"; exit;}
	}
	foreach(keys %reverse){
		unless($forward{$_}){print "Error!! Paird primer need to be the same name!\nCheck $primer\n"; exit;}
	}

	my $primer_num = keys %forward;
	if($primer_num == 1){$separate = 0;}
	if($trim and $primer_num == 1){
		foreach(keys %forward){$primerF = length($forward{$_}); $primerR = length($reverse{$_});}
	}elsif($trim and $primer_num > 1){
		print "Error! \"MaxDiff = length\" in Setting.txt can use when primer pair is one pair\n"; exit;
	}
	if($trim){$separate = 0;}
}

#Database and Usearch check
opendir (DIR, ".\/DataBase") or die ("error:$!");
my @database = readdir DIR;
my $dbcheck = 0;
foreach(@database){
	if ($_ =~ /${db}\.udb/){$dbcheck++;}
}
closedir DIR;

opendir (DIR, ".\/Tools") or die ("error:$!");
my @tool = readdir DIR;
my $usearch;
foreach (@tool) {
	if ($_ =~ /(usearch.+)/){$usearch = $1;}
}
closedir DIR;
unless($usearch){print "Error!! Please put an usearch exectable file in Tools directory!\n"; exit;}

unless($dbcheck){
	my $command = ".\/Tools\/$usearch -makeudb_usearch \".\/DataBase\/$db\" -output \".\/DataBase\/${db}\.udb\"";
	system $command;
}


#Decompressing files
opendir (DIR, ".\/Run") or die ("error:$!");
my @run = readdir DIR;
my $gzip = 0;
foreach (@run) {
	if ($_ =~ /\.gz$/){$gzip++;}
}
closedir DIR;

if($gzip){
	print "Decompressing gz file to fastq file...\n";
	foreach (@run) {
		if ($_ =~ /\.gz$/){
			my $comand = "gzip -d \".\/Run\/$_\"";
			system $comand;
		}
	}
}

#Get data names
opendir (DIR, ".\/Run") or die ("error:$!");
my @read = readdir DIR;
my $read_count = 0;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_R1/){$file{$1}++; $read_count++;}
	if ($_ =~ /(.+)_R2/){$file{$1}++; $read_count++;}
}
closedir DIR;
unless(%file){print "None of fastq files in Run file!!\n"; exit;}
foreach(keys %file){
	if($file{$_} < 2){print "Error!! Check the Run directory!The analysis needs both R1 and R2 fastq files.\nR1 or R2 fastq files of $_ is not.\n"; exit;}
}
unless(-d ".\/Results"){mkdir ".\/Results";}


print "============================================================\n";
print "                      1_Preprocessing                       \n";
print "============================================================\n";
#assemble pair seqs
mkdir ".\/Results\/1_1_Merge_paird_reads";
unless($primer =~ /^no$/i){mkdir ".\/Results\/1_2_Strip_primers";}
mkdir ".\/Results\/1_3_Quality_filter";
open (LOG, ">", ".\/Results\/log.txt") or die("error:$!");
my (%log, %raw);

my @youbi = ('Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat');
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
$mon += 1;
if($min < 10){$min = "0$min";}
print LOG "Start: $year/$mon/$mday ($youbi[$wday]) $hour:$min\n";

foreach(sort keys %file){
	my $file = $_;
	print "============================================================\n";
	print "$file    1_1_Merge_paird_reads                   \n";
	print "============================================================\n";
	my $file1;
	foreach(sort @read){
		if($_ =~ /${file}_R1[\._]/){$file1 = $_;}
	}
	my $comand = ".\/Tools\/$usearch -fastq_mergepairs \".\/Run\/$file1\" -fastqout \".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq\" -log \".\/Results\/1_1_Merge_paird_reads\/${file}_log.txt\"";
	system $comand;
	if($compress =~ /yes/i){
		print "\nCompressing fastq file to gz file...\n";
		foreach(@read){
			if($_ !~ /gz$/ and $_ =~ /$file/){
				$comand = "gzip \".\/Run\/$_\"";
				system $comand;
			}
		}
	}
	open (DATA, "<", ".\/Results\/1_1_Merge_paird_reads\/${file}_log.txt") or die("error:$!");
	while(<DATA>){
		chomp($_);
		$_ =~ s/\r//g;
		if($_ =~ /(\d+)\s+Pairs/){$log{$file} = $1; $raw{$file} = $1;}
		if($_ =~ /(\d+)\s+Merged.+\s([^\s]+)%/){$log{$file} = $log{$file} . "\t$1 ($2%)";}
	}
	close(DATA);
	print "\n\n";
	
	print "============================================================\n";
	print "$file    1_2_Strip_primers                      \n";
	print "============================================================\n";
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
	
	if($primer !~ /^no$/i){
		if($separate){
			opendir (DIR, ".\/Results\/1_2_Strip_primers\/") or die ("error:$!");
			my @reads = readdir DIR;
			my %divide;
			foreach (@reads) {
				if ($_ =~ /(${file}.+)_stripped.fq/){$divide{$1}++;}
			}
			closedir DIR;
			foreach(sort keys %divide){
				my $files = $_;
				print "\n\n";
				print "============================================================\n";
				print "$files    1_3_Quality_filter                      \n";
				print "============================================================\n";
				$comand = ".\/Tools\/$usearch -fastq_filter \".\/Results\/1_2_Strip_primers\/${files}_stripped.fq\" -fastq_maxee 1.0 -fastq_minlen $size -fastaout \".\/Results\/1_3_Quality_filter\/${files}_filtered.fa\" -log \".\/Results\/1_3_Quality_filter\/${files}_log.txt\"";
				system $comand;
				my $slog;
				open (DATA, "<", ".\/Results\/1_2_Strip_primers\/${files}_log.txt") or die("error:$!");
				while(<DATA>){
					if($_ =~ />\s+(\d+)\s+reads/){$slog = $1;}
				}
				close(DATA);
				open (DATA, "<", ".\/Results\/1_3_Quality_filter\/${files}_log.txt") or die("error:$!");
				while(<DATA>){
					chomp($_);
					$_ =~ s/\r//g;
					if($_ =~ /(\d+)\s+Filtered/){
						my $sper = sprintf("%.2f", $slog/$raw{$file}*100);
						my $per = sprintf("%.2f", $1/$raw{$file}*100);
						$log{$files} = $log{$file} . "\t$slog ($sper%)\t$1 ($per%)";
					}
				}
				close(DATA);
				if($temporary =~ /yes/i){unlink ".\/Results\/1_2_Strip_primers\/${files}_stripped.fq";}
			}
			undef($log{$file});
		}elsif(-f ".\/Results\/1_2_Strip_primers\/${file}_stripped.fq"){
			print "\n\n";
			print "============================================================\n";
			print "$file    1_3_Quality_filter                      \n";
			print "============================================================\n";
			$comand = ".\/Tools\/$usearch -fastq_filter \".\/Results\/1_2_Strip_primers\/${file}_stripped.fq\" -fastq_maxee 1.0 -fastq_minlen $size -fastaout \".\/Results\/1_3_Quality_filter\/${file}_filtered.fa\" -log \".\/Results\/1_3_Quality_filter\/${file}_log.txt\"";
			system $comand;
			unless($trim){
				my $slog;
				open (DATA, "<", ".\/Results\/1_2_Strip_primers\/${file}_log.txt") or die("error:$!");
				while(<DATA>){
					if($_ =~ />\s+(\d+)\s+reads/){$slog = $1;}
				}
				close(DATA);
				open (DATA, "<", ".\/Results\/1_3_Quality_filter\/${file}_log.txt") or die("error:$!");
				while(<DATA>){
					chomp($_);
					$_ =~ s/\r//g;
					if($_ =~ /(\d+)\s+Filtered/){
						my $sper = sprintf("%.2f", $slog/$raw{$file}*100);
						my $per = sprintf("%.2f", $1/$raw{$file}*100);
						$log{$file} = $log{$file} . "\t$slog ($sper%)\t$1 ($per%)";
					}
				}
				close(DATA);
				if($temporary =~ /yes/i){unlink ".\/Results\/1_2_Strip_primers\/${file}_stripped.fq";}
			}else{
				open (DATA, "<", ".\/Results\/1_3_Quality_filter\/${file}_log.txt") or die("error:$!");
				while(<DATA>){
					chomp($_);
					$_ =~ s/\r//g;
					if($_ =~ /(\d+)\s+Filtered/){my $per = sprintf("%.2f", $1/$raw{$file}*100); $log{$file} = $log{$file} . "\t$1 ($per%)";}
				}
				close(DATA);
				if($temporary =~ /yes/i){unlink ".\/Results\/1_2_Strip_primers\/${file}_stripped.fq";}
			}
		}
	}elsif($primer =~ /^no$/i and -f ".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq"){
		print "\n\n";
		print "============================================================\n";
		print "$file    1_3_Quality_filter                      \n";
		print "============================================================\n";
		$comand = ".\/Tools\/$usearch -fastq_filter \".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq\" -fastq_maxee 1.0 -fastq_minlen $size -fastaout \".\/Results\/1_3_Quality_filter\/${file}_filtered.fa\" -log \".\/Results\/1_3_Quality_filter\/${file}_log.txt\"";
		system $comand;
		open (DATA, "<", ".\/Results\/1_3_Quality_filter\/${file}_log.txt") or die("error:$!");
		while(<DATA>){
			chomp($_);
			$_ =~ s/\r//g;
			if($_ =~ /(\d+)\s+Filtered/){my $per = sprintf("%.2f", $1/$raw{$file}*100); $log{$file} = $log{$file} . "\t$1 ($per%)";}
		}
		close(DATA);
		if($temporary =~ /yes/i){unlink ".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq"; unlink ".\/Results\/1_1_Merge_paird_reads\/${file}_log.txt";}
	}
	print "\n\n";
}
foreach(sort keys %log){
	if($log{$_}){print LOG "$_\t$log{$_}\n";}
}
close(LOG);

if($temporary =~ /yes/i){rmdir ".\/Results\/1_1_Merge_paird_reads";}
