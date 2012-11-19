package Low_Complexity;
# Low_Complexity.pm
# library for identifying regions of low complexity in DNA sequence
use strict; use warnings;

sub low_comp {
	
	my ($seq, $window, $threshold) = @_; #need die statement in program
	#If window size is larger than remaining sequence, will use largest possible window
	if ($threshold<0 || $threshold>2){
		die "Threshold must have a range between 0 and 2";
	}
	$seq = uc $seq;
	my ($seq_masked) = "";

	for (my $i = 0; $i < length($seq); $i += $window) {
		my $subseq = substr($seq, $i, $window); 
	
		#calculate and store H
		my $r = $subseq =~ tr/R/R/;
		my $y = $subseq =~ tr/Y/Y/;
		my $k = $subseq =~ tr/K/K/;
		my $m = $subseq =~ tr/M/M/;
		my $s = $subseq =~ tr/S/S/;
		my $w = $subseq =~ tr/W/W/;
		my $b = $subseq =~ tr/B/B/;
		my $d = $subseq =~ tr/D/D/;
		my $fh = $subseq =~ tr/H/H/; # FASTA H so as not to confuse w/ entropy
		my $v = $subseq =~ tr/V/V/;
		my $n = $subseq =~ tr/N/N/;
		
		my $a = $subseq =~ tr/A/A/;
		my $c = $subseq =~ tr/C/C/;
		my $g = $subseq =~ tr/G/G/;
		my $t = $subseq =~ tr/T/T/;
		
		my $pseudo_count = 0.0001;
				
		$a += ($r/2 + $m/2 + $w/2 + $n/4 + $d/3 + $fh/3 + $v/3 + $pseudo_count);
		$c += ($y/2 + $m/2 + $s/2 + $b/3 + $fh/3 + $v/3 + $n/4 + $pseudo_count);
		$g += ($r/2 + $k/2 + $s/2 + $b/3 + $d/3 + $v/3 + $n/4 + $pseudo_count);
		$t += ($y/2 + $k/2 + $w/2 + $b/3 + $d/3 + $fh/3 + $n/4 + $pseudo_count);
		
		my $total = ($a + $c + $g + $t);
		$a /= $total;
		$c /= $total;
		$g /= $total;
		$t /= $total;
		
		# print "$a \t $c \t $g \t $t \n";
		
		my $h = - ($a * log($a) + $c * log($c) + $g * log($g) + $t * log($t));
		
		if ($h < $threshold) {
			$subseq =~ tr/A-Z/a-z/;
		}
		
		$seq_masked .= $subseq;
	}
	return $seq_masked;
}

1;