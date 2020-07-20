#!/usr/bin/perl
use strict;
no strict "refs";
use warnings;

#trim primer region
#not account for indels

my ($file, $primer, $diff, $out, $separate, $log);
for (my $v = 0; $v < @ARGV; $v++){
	if($ARGV[$v] =~ /-f$/){$file = $ARGV[$v+1];}
	elsif($ARGV[$v] =~ /-db$/){$primer = $ARGV[$v+1];}
	elsif($ARGV[$v] =~ /-diff$/){$diff = $ARGV[$v+1];}
	elsif($ARGV[$v] =~ /-o$/){$out = $ARGV[$v+1];}
	elsif($ARGV[$v] =~ /-s$/){$separate = $ARGV[$v+1];}
	elsif($ARGV[$v] =~ /-l$/){$log = 1;}
}
unless($file){print "-f option error!\n"; exit;}
unless($primer){print "-db option error!\n"; exit;}
unless($diff =~ /\d+/){print "-diff option error!\n"; exit;}
unless($separate =~ /\d+/){print "-s option error!\n"; exit;}
my $shift = 1;

$file =~ /\/([^\/]+)_assembled/;
my $file2 = $1;

#get reference
open (DATA2, "<", $primer) or die("error:$!");
my $count = 0;
my (%forward, %reverse, $name, @primername);
while(<DATA2>){
	if($_ =~ /\#/){next;}
	$_ =~ s/\r//g;
	chomp($_);
	if($_ =~ /Forward/){$count = 1;next;}
	if($_ =~ /Reverse/){$count = 2;next;}
	if($count == 1){
		if($_ =~ /^>(.+)/){$name = $1;push(@primername, $1);}
		elsif($_ =~ /^[a-z]/i){$_ =~ tr/[a-z]/[A-Z]/; $forward{$name} = $_;}
	}
	if($count == 2){
		if($_ =~ /^>(.+)/){$name = $1;}
		elsif($_ =~ /^[a-z]/i){
			my $rev = reverse($_);
			$rev =~ tr/ATGCURYMKDHBV/TACGAYRKMHDVB/;
			$reverse{$name} = $rev;}
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

#primer set
my (%fcut, %rcut, %fnc, %rnc, %nnum, %nnumr, %flength, %rlength, %plistf, %plistr);
foreach(@primername){
	my $pn = $_;
	my $forward_cut = $forward{$_};
	my $reverse_cut = $reverse{$_};
	$forward_cut =~ s/^N+//;
	$reverse_cut =~ s/N+$//;
	$fnc{$_} = $forward_cut;
	$rnc{$_} = $reverse_cut;
	if($forward_cut =~ /[RYMKDHBVNSW]/){$fcut{$_} = &degenerate($forward_cut);}
	else{$fcut{$_} = $forward_cut;}
	if($reverse_cut =~ /[RYMKDHBVNSW]/){$rcut{$_} = &degenerate($reverse_cut);}
	else{$rcut{$_} = $reverse_cut;}
	
	if($forward{$_} =~ /^(N+)/){$nnum{$_} = length($1);}
	else{$nnum{$_} = 0;}
	$flength{$_} = length($forward{$_});
	
	if($reverse{$_} =~ /(N+)$/){$nnumr{$_} = length($1);}
	else{$nnumr{$_} = 0;}
	$rlength{$_} = length($reverse{$_});

	my @temp = split(//, $fnc{$_});
	for(my $n = 0; $n < @temp; $n++){
		my $moji = $temp[$n];
		unless($moji =~ /[ATGC]/){$moji = &degenerate($moji);}
		$plistf{$pn}{$n} = $moji;
	}
	undef(@temp);
	@temp = split(//, $rnc{$_});
	for(my $n = 0; $n < @temp; $n++){
		my $moji = $temp[$n];
		unless($moji =~ /[ATGC]/){$moji = &degenerate($moji);}
		$plistr{$pn}{$n} = $moji;
	}
}

open (FILE, "<", $file) or die("error:$!");
$count = 0;
my @data;
my $total = 0;
my $after = 0;
my (%kakunin, @mismatch, %pipe, %unlink, %unlink2, %unlink3, %unlink4);

if($separate and $out){
	$out =~ s/_stripped.fq//;
	foreach(@primername){
		my $pipe = $_;
		open ($pipe, ">", "${out}_${_}_stripped.fq") or die("error:$!");
	}
}elsif($out){
	open (OUT, ">", "$out") or die("error:$!");
}else{
	open (OUT, ">", "out.fq") or die("error:$!");
}

while(<FILE>){
	chomp($_);
	$_ =~ s/\r//g;
	unless($count == 3){push(@data, $_); $count++;}
	else{
		$total++;
		push(@data, $_);
		my $next=1;
		foreach(@primername){
			my $pname = $_;
			unless($data[1] =~ /$fcut{$_}/){next;}
			unless($data[1] =~ /$rcut{$_}/){next;}
			if($data[1] =~ /^(\w*)$fcut{$_}(\w+)$rcut{$_}(\w*)$/){
				$data[1] = $2;
				my ($fn, $rn); 
				unless($1){$fn = 0;}
				else{$fn = length($1);}
				unless($3){$rn = 0;}
				else{$rn = length($3);}
				if($fn > $nnum{$_} + $shift or $fn < $nnum{$_} -$shift){next;}
				if($rn > $nnumr{$_} + $shift or $rn < $nnumr{$_} -$shift){next;}
				$after++;
				if($fn != $nnum{$_} or $rn != $nnumr{$_}){$unlink3{$_}++;}
				my $length = length($data[1]);
				my $posi = $fn + length($fcut{$_});
				$data[3] = substr($data[3], $posi, $length);
				$unlink{$pname}++;
				foreach(@data){
					if($separate){
						print $pname "$_\n";
					}else{
						print OUT "$_\n";
					}
				}
				$next = 0;
				$kakunin{$_}++;
				last;
			}
		}
		if($next){
			my @temp = @data;
			push(@mismatch, \@temp);
		}
		undef(@data);
		$count = 0;
	}
}
close(FILE);
my %unlink1 = %unlink;

my $hit = 0;
foreach(@primername){
	if($kakunin{$_}){$hit++;}
}

if($hit and $diff > 0){
	foreach(@mismatch){
		my @data = @$_;
		foreach(@primername){
			my $check = 0;
			my $pname = $_;
			unless($kakunin{$_}){next;}

			my $query = substr($data[1], 0, $flength{$_});
			my $query_cut = substr($data[1], $nnum{$_}, $flength{$_}-$nnum{$_});
			
			if(length($query) < $flength{$_}){next;}
			my @flist = split(//, $query_cut);
			
			my $count1 = 0;
			my $v = 0;
			foreach(@flist){
				unless($_ =~ /$plistf{$pname}{$v}/){$count1++;}
				$v++;
				if($count1 > $diff){last;}
			}

			#Shift check
			if($count1 > $diff){
				for(my $m = -$shift; $m <= $shift; $m++){
					$count1 = 0; $v = 0;
					undef(@flist);
					my $temp_query;
					if($m == 0){next;}
					$temp_query = substr($data[1], $nnum{$_}-$m, $flength{$_}-$nnum{$_});
					@flist = split(//, $temp_query);
					foreach(@flist){
						unless($_ =~ /$plistf{$pname}{$v}/){$count1++;}
						$v++;
						if($count1 > $diff){last;}
					}
					$query = substr($data[1], 0, $flength{$_}+$m);
					$check = 1;
					if($count1 <= $diff){last;}
				}
			}
			
			#Reverse check
			if($count1 <= $diff){
				my $rquery = substr($data[1], -$rlength{$_});
				my $rquery_cut;
				unless($nnumr{$_}){$rquery_cut = substr($data[1], -$rlength{$_});}
				else{$rquery_cut = substr($data[1], -$rlength{$_}, -$nnumr{$_});}
				if(length($rquery) < $rlength{$_}){next;}
				my @rlist = split(//, $rquery_cut);

				my $count2 = 0;
				$v = 0;
				foreach(@rlist){
					unless($_ =~ /$plistr{$pname}{$v}/){$count2++;}
					$v++;
					if($count2 > $diff){last;}
				}
				
				#Shift check
				if($count2 > $diff){
					for(my $m = -$shift; $m <= $shift; $m++){
						$count2 = 0; $v = 0;
						undef(@rlist);
						my $temp_query;
						if($m == 0){next;}
						$temp_query = substr($data[1], -$rlength{$_}+$m, -$nnumr{$_}+$m);
						@rlist = split(//, $temp_query);
						foreach(@rlist){
							unless($_ =~ /$plistr{$pname}{$v}/){$count2++;}
							$v++;
							if($count2 > $diff){last;}
						}
						$rquery = substr($data[1], -$rlength{$_}+$m);
						$check = 1;
						if($count2 <= $diff){last;}
					}
				}
				
				if($count2 <= $diff){
					if($data[1] =~ /(\w*)$query(\w+)$rquery/){
						$after++;
						$data[1] = $2;
						my $length = length($data[1]);
						my $posi = length($1) + length($query);
						$data[3] = substr($data[3], $posi, $length);
						$unlink{$pname}++; 
						if($check == 1){$unlink4{$pname}++;}
						$unlink2{$pname}++;
						foreach(@data){
							if($separate){
								print $pname "$_\n";
							}else{
								print OUT "$_\n";
							}
						}
					}
					last;
				}
			}
		}
		undef(@data);
	}
}
print $total-$after, " reads were discarded ($after/$total)\n";
foreach(@primername){
	if($unlink{$_}){
		print "\tPrimer $_ = Total ",$unlink{$_}, " reads\n";
		if($unlink3{$_}){
			print "\t Diff. = 0 : $unlink1{$_} reads (Shifted $unlink3{$_} reads)\n";
		}else{
			print "\t Diff. = 0 : $unlink1{$_} reads (Shifted 0 reads)\n";
		}
		unless($unlink2{$_}){print "\t Diff. > 0 : 0 reads\n";}
		else{
			unless($unlink4{$_}){print "\t Diff. > 0 : $unlink2{$_} reads (Shifted 0 reads)\n";}
			else{print "\t Diff. > 0 : $unlink2{$_} reads (Shifted $unlink4{$_} reads)\n";}
		}
	}
}

if($separate and $out){
	foreach(@primername){
		my $key = $_;
		close($key);
		unless($unlink{$_}){unlink "${out}_${_}_stripped.fq";}
	}
	my @temp = keys %unlink;
	my $filenum = @temp;
	if($filenum > 1){
		print "\nThe fq file divided into $filenum fq files.\n";
		foreach(keys %unlink){
			print "  ${file2}_${_}_stripped.fq  > $unlink{$_} reads\n";
			if($log){
				open (LOG, ">", "${out}_${_}_log.txt") or die("error:$!");
				print LOG "${file2}_${_}_stripped.fq  > $unlink{$_} reads\n";
				print LOG "\tPrimer $_ = Total ",$unlink{$_}, " reads\n";
				if($unlink3{$_}){
					print LOG "\t Diff. = 0 : $unlink1{$_} reads (Shifted $unlink3{$_} reads)\n";
				}else{
					print LOG "\t Diff. = 0 : $unlink1{$_} reads (Shifted 0 reads)\n";
				}
				unless($unlink2{$_}){print LOG "\t Diff. > 0 : 0 reads\n";}
				else{
					unless($unlink4{$_}){print LOG "\t Diff. > 0 : $unlink2{$_} reads (Shifted 0 reads)\n";}
					else{print LOG "\t Diff. > 0 : $unlink2{$_} reads (Shifted $unlink4{$_} reads)\n";}
				}
				close(LOG);
			}
		}
	}elsif($filenum == 1){
		if($log){
			open (LOG, ">", "${out}_$temp[0]_log.txt") or die("error:$!");
			print LOG "${file2}_$temp[0]_stripped.fq  > $unlink{$temp[0]} reads\n";
			print LOG "\tPrimer $temp[0] = Total ",$unlink{$temp[0]}, " reads\n";
			if($unlink3{$temp[0]}){
				print LOG "\t Diff. = 0 : $unlink1{$temp[0]} reads (Shifted $unlink3{$temp[0]} reads)\n";
			}else{
				print LOG "\t Diff. = 0 : $unlink1{$temp[0]} reads (Shifted 0 reads)\n";
			}
			unless($unlink2{$temp[0]}){print LOG "\t Diff. > 0 : 0 reads\n";}
			else{
				unless($unlink4{$temp[0]}){print LOG "\t Diff. > 0 : $unlink2{$temp[0]} reads (Shifted 0 reads)\n";}
				else{print LOG "\t Diff. > 0 : $unlink2{$temp[0]} reads (Shifted $unlink4{$temp[0]} reads)\n";}
			}
			close(LOG);
		}
	}
}elsif($out){
	if($log){
		$out =~ s/_stripped.fq//;
		open (LOG, ">", "${out}_log.txt") or die("error:$!");
		print LOG "${file2}_stripped.fq  > $after reads\n";
		foreach(@primername){
			if($unlink{$_}){
				print LOG "\tPrimer $_ = Total ",$unlink{$_}, " reads\n";
				if($unlink3{$_}){
					print LOG "\t Diff. = 0 : $unlink1{$_} reads (Shifted $unlink3{$_} reads)\n";
				}else{
					print LOG "\t Diff. = 0 : $unlink1{$_} reads (Shifted 0 reads)\n";
				}
				unless($unlink2{$_}){print LOG "\t Diff. > 0 : 0 reads\n";}
				else{
					unless($unlink4{$_}){print LOG "\t Diff. > 0 : $unlink2{$_} reads (Shifted 0 reads)\n";}
					else{print LOG "\t Diff. > 0 : $unlink2{$_} reads (Shifted $unlink4{$_} reads)\n";}
				}
			}
		}
		close(LOG);
	}
	close(OUT);
}else{
	if($log){
		open (LOG, ">", "log.txt") or die("error:$!");
		print LOG "out.fq  > $after reads\n";
		foreach(@primername){
			if($unlink{$_}){
				print LOG "\tPrimer $_ = Total ",$unlink{$_}, " reads\n";
				if($unlink3{$_}){
					print LOG "\t Diff. = 0 : $unlink1{$_} reads (Shifted $unlink3{$_} reads)\n";
				}else{
					print LOG "\t Diff. = 0 : $unlink1{$_} reads (Shifted 0 reads)\n";
				}
				unless($unlink2{$_}){print LOG "\t Diff. > 0 : 0 reads\n";}
				else{
					unless($unlink4{$_}){print LOG "\t Diff. > 0 : $unlink2{$_} reads (Shifted 0 reads)\n";}
					else{print LOG "\t Diff. > 0 : $unlink2{$_} reads (Shifted $unlink4{$_} reads)\n";}
				}
			}
		}
		close(LOG);
	}
}


#Degenerate
sub degenerate {
	my ($dege) = @_;
	$dege =~ s/B/[CGT]/g; $dege =~ s/D/[AGT]/g; $dege =~ s/H/[ACT]/g; $dege =~ s/K/[GT]/g; 
	$dege =~ s/M/[AC]/g; $dege =~ s/N/[ACGT]/g; $dege =~ s/R/[AG]/g; $dege =~ s/S/[CG]/g;
	$dege =~ s/V/[ACG]/g; $dege =~ s/W/[AT]/g; $dege =~ s/Y/[CT]/g;
	return($dege);
}
