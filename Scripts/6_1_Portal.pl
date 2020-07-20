#!/usr/bin/perl
use strict;
use warnings;

#Setting
#Setting
my ($db, $primer, $diff, $trim, $separate, $depth, $size, $identity, $identity2, $dictionary, $family, $temporary, $compress, $rarefy, $timing);
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
	elsif($_ =~ /^Rarefaction\s*=\s*(\S+)/){$rarefy = $1;}
	elsif($_ =~ /^Timing\s*=\s*(\d+)/){$timing = $1;}
	elsif($_ =~ /^Temporary\s*=\s*(\S+)/){$temporary = $1;}
	elsif($_ =~ /^Compress\s*=\s*(\S+)/){$compress = $1;}
}
close(SET);
if($diff =~ /^length$/i){$trim = 1;}
unless($diff =~ /\d+/){$diff = 2;}
unless($separate){$separate = 1;}
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
unless(-f ".\/Dictionary\/$dictionary"){
	if($dictionary =~ /^no$/i){undef($dictionary);}
	else{print "The file of Dictionary for common name isn't or don't match the name in Setting.txt\n"; exit;}
}


#Get data names
opendir (DIR, ".\/Results\/4_1_Annotation") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_Representative_seq/){$file{$1}++;}
}
closedir DIR;


print "============================================================\n";
print "                        6_1_Portal                          \n";
print "============================================================\n";

my @youbi = ('Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat');
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
$mon += 1;
if($min < 10){$min = "0$min";}

my $start;
if(-f ".\/Results\/log.txt"){
	open (OUT, "<", ".\/Results\/log.txt") or die("error:$!");
	while(<OUT>){
		chomp($_);
		if($_ =~ /Start/){$start = $_; last;}
	}
	close(OUT);
}else{
	$start = "Start: no data";
}

#Make portal
open (OUT, ">", ".\/Results\/Portal.html") or die("error:$!");
print OUT <<"EOS";
<html>
<head>
<meta http-equiv="Content-type" content="text/html" charset="Shift_JIS">
<title>Portal</title>
<style type="text/css">
H1{color: #ffffff; text-align: left; font: 100% Tahoma;}
H2{color: #000000; text-align: left; font: 100% Tahoma;}
H3{text-align: left; font: 100% Tahoma; line-height: 4px;}
</style>
</head>
<body bgcolor="#EFEFFB">
<a name="top"></a>
<font face="Tahoma" size="3">$start</font><BR>
<font face="Tahoma" size="3">End : $year/$mon/$mday ($youbi[$wday]) $hour:$min</font><BR>
<font face="Tahoma" size="4"><a href="./log.html">Log (Reads after each step)</a></font><BR>
<font face="Tahoma" size="6">Summary of Analysis</font><BR>
<font face="Tahoma" size="4"><a href="./5_2_Summary_table/Summary_table.tsv">Summary Table</a></font><BR>
<table><table border="0" cellspacing="1" bgcolor="#191970" border="1" align="left">
<tr>
<th bgcolor="#ff0000" align="left" nowrap><H1>Samples analyzed &nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>Total read #&nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>Identity >= ${identity}% species #&nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>${identity2}% < Identity < ${identity}% species #&nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>Nohit seq. #&nbsp;&nbsp;</H1></th>
</tr>
EOS

my $matome = 0;
foreach(sort keys %file){
	my $file = $_;
	open (DATA, "<", ".\/Results\/4_1_Annotation\/${file}_Representative_seq.fas") or die("error:$!");
	my $count = 0;
	my $tread = 0;
	my $lower = 0;
	my $nohit = 0;
	$matome++;
	while(<DATA>){
		if($_ =~ /^>/){
			$count++;
			if($_ =~ /^>U${identity}_/){$lower++;}
			elsif($_ =~ /^>Nohit/){$nohit++;}
			$_ =~ /(\d+)_reads/;
			$tread  = $tread + $1;
		}
	}
	close(DATA);
	my $upper = $count - $lower - $nohit;
	if($matome%2 + 1 == 2){print OUT "<tr>\n<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><a href=\#$file><H2>$file\&nbsp\;\&nbsp\;</H2></a></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$tread\&nbsp\;\&nbsp\;</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$upper</H2></td></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$lower</H2></td></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$nohit</H2></td>\n</tr>\n";}
	else{print OUT "<tr>\n<td bgcolor=\"\#fceded\" align=\"left\" nowrap><a href=\#$file><H2>$file\&nbsp\;\&nbsp\;</H2></a></td><td bgcolor=\"\#fceded\" align=\"left\" nowrap><H2>$tread\&nbsp\;\&nbsp\;</H2></td><td bgcolor=\"\#fceded\" align=\"left\" nowrap><H2>$upper</H2></td></td><td bgcolor=\"\#fceded\" align=\"left\" nowrap><H2>$lower</H2></td></td><td bgcolor=\"\#fceded\" align=\"left\" nowrap><H2>$nohit</H2></td>\n</tr>\n";}
}
print OUT "</table>\n<br clear = \"all\">\n";
print OUT "<a href=\#top><font face=\"Tahoma\">Back to top<font></a><br>\n";
print OUT "<hr style=\"border:0;border-top:thick dotted white;background-color:red;\">\n";
print OUT "<font face=\"Tahoma\" size=\"6\">Each Result of Sample</font><BR><br>\n";


foreach(sort keys %file){
	my $file = $_;
	print OUT "<a name=\"$_\"><font face=\"Tahoma\" size=\"5\"><B>$_</B></font></a><BR>\n<font face=\"Tahoma\" size=\"4\"><a href=\".\/4_1_Annotation\/${file}_Detail.html\">Detailed Results</a></font>\&nbsp\;\&nbsp\;\&nbsp\;\&nbsp\;";
	print OUT "<font face=\"Tahoma\" size=\"4\"><a href=\".\/4_1_Annotation\/${file}_Representative_seq.fas\">Representative Sequences</a></font><br>\n";
	open (DATA, "<", ".\/Results\/4_1_Annotation\/${file}_Summary.txt") or die("error:$!");
	my $table = 0;
	my $iro = 0;
	while(<DATA>){
		if($_ =~ /^\n/){next;}
		chomp($_);
		if($_ =~ /^Identity/){print OUT "<font face=\"Tahoma\" size=\"3\"><B>$_</B></font><BR>\n";$table = 1;next;}
		if($_ =~ /^${identity2}%/){print OUT "</table>\n<br clear = \"all\"><br>\n<font face=\"Tahoma\" size=\"3\"><B>$_</B></font><BR>\n";$iro = 0;$table = 2;next;}
		if($table == 1){
			if($_ =~ /^No applicable sequence/){print OUT "<font face=\"Tahoma\" size=\"2\"><B>$_</B></font><BR>\n"; $table = 2;next;}
			if($_ =~ /^Scientific/){&table;next;}
			my @temp = split(/\t/, $_);
			$iro++;
			if($dictionary){&nakami1($iro, \@temp, $file);}
			else{&nakami2($iro, \@temp, $file);}
		}
		if($table == 2){
			if($_ =~ /^No applicable sequence/){print OUT "<font face=\"Tahoma\" size=\"2\"><B>$_</B></font><BR>\n"; next;}
			if($_ =~ /^Nohit/){print OUT "</table>\n<br clear = \"all\"><br>\n<font face=\"Tahoma\" size=\"3\"><B>$_</B></font><br><br>\n<a href=\#top><font face=\"Tahoma\">Back to list<font></a><br><br><br>\n\n\n\n"; last;}
			if($_ =~ /^Scientific/){&table;next;}
			my @temp = split(/\t/, $_);
			$iro++;
			if($dictionary){&nakami1($iro, \@temp, $file);}
			else{&nakami2($iro, \@temp, $file);}
		}
	}
	close(DATA);
}
close(OUT);
print "Portal.html was created.\n";

#Make log.html
open (OUT, ">", ".\/Results\/log.html") or die("error:$!");
print OUT <<"EOS";
<html>
<head>
<meta http-equiv="Content-type" content="text/html" charset="Shift_JIS">
<title>Log</title>
<style type="text/css">
H1{color: #ffffff; text-align: left; font: 100% Tahoma;}
H2{color: #000000; text-align: left; font: 100% Tahoma;}
H3{text-align: left; font: 100% Tahoma; line-height: 4px;}
</style>
</head>
<body bgcolor="#EFEFFB">
<a name="top"></a>
<font face="Tahoma" size="6">Reads after each step</font><BR>
<table><table border="0" cellspacing="1" bgcolor="#191970" border="1" align="left">
<tr>
<th bgcolor="#ff0000" align="left" nowrap><H1>Samples analyzed &nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>Raw read #&nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>Merged &nbsp;&nbsp;</H1></th>
EOS

if($trim){$primer = "no";}
unless($primer =~ /^no$/i){
	print OUT "<th bgcolor=\"#ff0000\" align=\"left\" nowrap><H1>Strip primer &nbsp;&nbsp;</H1></th>\n";
}
print OUT "<th bgcolor=\"#ff0000\" align=\"left\" nowrap><H1>Quality filter &nbsp;&nbsp;</H1></th>\n";
unless($rarefy =~ /^no$/i){
	if($timing == 1){
		print OUT "<th bgcolor=\"#ff0000\" align=\"left\" nowrap><H1>Rarefaction &nbsp;&nbsp;</H1></th>\n";
		print OUT "<th bgcolor=\"#ff0000\" align=\"left\" nowrap><H1>Denoise &nbsp;&nbsp;</H1></th>\n</tr>\n";
	}elsif($timing == 2){
		print OUT "<th bgcolor=\"#ff0000\" align=\"left\" nowrap><H1>Denoise &nbsp;&nbsp;</H1></th>\n";
		print OUT "<th bgcolor=\"#ff0000\" align=\"left\" nowrap><H1>Rarefaction &nbsp;&nbsp;</H1></th>\n</tr>\n";
	}
}else{
	print OUT "<th bgcolor=\"#ff0000\" align=\"left\" nowrap><H1>Denoise &nbsp;&nbsp;</H1></th>\n</tr>\n";
}

open (LOG, "<", ".\/Results\/log.txt") or die("error:$!");
my $count = 0;
while(<LOG>){
	chomp($_);
	if($_ =~ /Start/){next;}
	$_ =~ s/\r//g;
	my @log = split(/\t/, $_);
	unless($count){
		print OUT "</tr>\n";
		foreach(@log){print OUT "<td bgcolor=\"#ffffff\" align=\"left\" nowrap><H2>$_&nbsp;&nbsp;</H2></td>"}
		print OUT "<tr>\n";
		$count = 1;
	}else{
		print OUT "</tr>\n";
		foreach(@log){print OUT "<td bgcolor=\"#fceded\" align=\"left\" nowrap><H2>$_&nbsp;&nbsp;</H2></td>"}
		print OUT "<tr>\n";
		$count = 0;
	}
}
close(LOG);
print OUT "</table>\n<br clear = \"all\">\n";
unlink ".\/Results\/log.txt";

#sub
sub table{
	if($dictionary){
		print OUT <<'EOS';
<table><table border="0" cellspacing="1" bgcolor="#191970" border="1" align="left">
<tr>
<th bgcolor="#00008B" align="left" nowrap><H1>Scientific name &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Common name &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Total read #&nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Confidence &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Identity(%) &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Synonym &nbsp;&nbsp;</H1></th>
</tr>
EOS
	}else{
		print OUT <<'EOS';
<table><table border="0" cellspacing="1" bgcolor="#191970" border="1" align="left">
<tr>
<th bgcolor="#00008B" align="left" nowrap><H1>Scientific name &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Total read #&nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Confidence &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Identity(%) &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Synonym &nbsp;&nbsp;</H1></th>
</tr>
EOS
	}
}

sub nakami1{
	my ($iro, $temp, $file) = @_;
	my @temp = @$temp;
	my $synonym;
	if($temp[5]){$synonym = "align=\"left\" nowrap><H2><a href=\".\/4_1_Annotation\/${file}_Synonym_list.html\">$temp[5]</a></H2></td>";}
	else{$synonym = "align=\"left\" nowrap><H2></H2></td>";}
	if($iro%2 + 1 == 2){
		print OUT <<"EOS";
<tr>
<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2><I>$temp[0]\&nbsp\;\&nbsp\;</I></H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[1]\&nbsp\;\&nbsp\;</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[2]</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[3]</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[4]</H2></td><td bgcolor=\"\#ffffff\" $synonym
</tr>
EOS
	}else{
		print OUT <<"EOS";
<tr>
<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2><I>$temp[0]\&nbsp\;\&nbsp\;</I></H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[1]\&nbsp\;\&nbsp\;</H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[2]</H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[3]</H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[4]</H2></td><td bgcolor=\"\#E6E6FA\" $synonym
</tr>
EOS
	}
}

sub nakami2{
	my ($iro, $temp, $file) = @_;
	my @temp = @$temp;
	my $synonym;
	if($temp[4]){$synonym = "align=\"left\" nowrap><H2><a href=\".\/4_1_Annotation\/${file}_Synonym_list.html\">$temp[4]</a></H2></td>";}
	else{$synonym = "align=\"left\" nowrap><H2></H2></td>";}
	if($iro%2 + 1 == 2){
		print OUT <<"EOS";
<tr>
<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2><I>$temp[0]\&nbsp\;\&nbsp\;</I></H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[1]\&nbsp\;\&nbsp\;</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[2]</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[3]</H2></td><td bgcolor=\"\#ffffff\" $synonym
</tr>
EOS
	}else{
		print OUT <<"EOS";
<tr>
<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2><I>$temp[0]\&nbsp\;\&nbsp\;</I></H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[1]\&nbsp\;\&nbsp\;</H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[2]</H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[3]</H2></td><td bgcolor=\"\#E6E6FA\" $synonym
</tr>
EOS
	}
}
