#! /usr/bin/perl -w

use strict;
#my $inputfile = shift;
#open(INPUT,"<$inputfile") or die "Unable to open $inputfile";
#my $outputfile = $inputfile;
#substr($outputfile,-4,0) = '_sbhpair';
#open(OUTPUT,">$outputfile") or die "Unable to write $outputfile";
#while(<INPUT>){
my %previousreadset = ();
while(<>){
	my $currentread = $_;
	if($currentread =~ /([^ \t\n\r]+)\t+([0-9]+)\t([^ \t\n\r@=]+)\t([0-9]+)\t([0-9]+)\t(([0-9]+[MIDNSHP])+|\*)\t([^ \t\n\r@]+)\t([0-9]+)\t(-?[0-9]+)\t([acgtnACGTN.=]+|\*)\t([!-~]+|\*)((\t[A-Z][A-Z0-9]:[AifZH]:[^\t\n\r]+)*)/){
		my($qname,$flag,$rname,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual,$tuples) = ($1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,$13);
		if($tuples =~ /\tRG:Z:(\w+\.?\d?)/){
			my $readgroup = $1;
			if(defined($previousreadset{'readname'}) && $previousreadset{'readname'} eq $qname){
				if(defined($previousreadset{$readgroup})){
					if($previousreadset{$readgroup} =~ /([^ \t\n\r]+)\t+([0-9]+)\t([^ \t\n\r@=]+)\t([0-9]+)\t([0-9]+)\t(([0-9]+[MIDNSHP])+|\*)\t([^ \t\n\r@]+)\t([0-9]+)\t(-?[0-9]+)\t([acgtnACGTN.=]+|\*)\t([!-~]+|\*)((\t[A-Z][A-Z0-9]:[AifZH]:[^\t\n\r]+)*)/){
						my($s2qname,$s2flag,$s2rname,$s2pos,$s2mapq,$s2cigar,$s2mrnm,$s2mpos,$s2isize,$s2seq,$s2qual,$s2tuples) = ($1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,$13);
						if((($flag & 4) == 0) && (($s2flag & 4)==0) && ($tuples =~ /\tX0:i:1\t/) && ($tuples !~ /\tX1:i:[1-9]/) && ($s2tuples =~ /\tX0:i:1\t/) && ($s2tuples !~ /\tX1:i:[1-9]/)){
							print $previousreadset{$readgroup};
							print $currentread;
						}
					}
					undef $previousreadset{$readgroup};
				} else {
					$previousreadset{$readgroup} = $currentread;
				}
			} else {
				# either a new set of querynames or first entry
				%previousreadset = ();
				$previousreadset{'readname'} = $qname;
				$previousreadset{$readgroup} = $currentread;
			}
		} else {
			die "Samfile does not have readgroup information attached.\n";
		}
		#my $secondmate = (<INPUT>);
		#if($secondmate =~ /([^ \t\n\r]+)\t+([0-9]+)\t([^ \t\n\r@=]+)\t([0-9]+)\t([0-9]+)\t(([0-9]+[MIDNSHP])+|\*)\t([^ \t\n\r@]+)\t([0-9]+)\t(-?[0-9]+)\t([acgtnACGTN.=]+|\*)\t([!-~]+|\*)((\t[A-Z][A-Z0-9]:[AifZH]:[^\t\n\r]+)*)/){
		#	my($s2qname,$s2flag,$s2rname,$s2pos,$s2mapq,$s2cigar,$s2mrnm,$s2mpos,$s2isize,$s2seq,$s2qual,$s2tuples) = ($1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,$13);
		#	if($qname eq $s2qname){
		#		if((($flag & 4) == 0) && (($s2flag & 4)==0) && ($tuples =~ /\tX0:i:1\t/) && ($tuples !~ /\tX1:i:[1-9]/) && ($s2tuples =~ /\tX0:i:1\t/) && ($s2tuples !~ /\tX1:i:[1-9]/)){
		#			#print OUTPUT $firstmate;
		#			#print OUTPUT $secondmate;
		#			print $firstmate;
		#			print $secondmate;
		#		}
		#	} else {
		#		die "The paired data is not paired. $qname =/= $s2qname";
		#	}
		#} else {
		#	die "The sam file is not paired or something else is wrong with matching the fields";
		#}
	} else {
		#print OUTPUT $_ ;
		print $_ ;
	}

} 
close INPUT;
#close OUTPUT;
exit 0;
