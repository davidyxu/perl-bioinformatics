package Genome;
use strict; use warnings;
my %code = ( 'AAA' => 'K', 'AAC' => 'N', 'AAG' => 'K',
			 'AAT' => 'N', 'ACA' => 'T', 'ACC' => 'T',
			 'ACG' => 'T', 'ACT' => 'T', 'AGA' => 'R',
			 'AGC' => 'S', 'AGG' => 'R', 'AGT' => 'S',
			 'ATA' => 'I', 'ATC' => 'I', 'ATG' => 'M',
			 'ATT' => 'I', 'CAA' => 'Q', 'CAC' => 'H',
			 'CAG' => 'Q', 'CAT' => 'H', 'CCA' => 'P',
			 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',
			 'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R',
			 'CGT' => 'R', 'CTA' => 'L', 'CTC' => 'L',
			 'CTG' => 'L', 'CTT' => 'L', 'GAA' => 'E',
			 'GAC' => 'D', 'GAG' => 'E', 'GAT' => 'D',
			 'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A',
			 'GCT' => 'A', 'GGA' => 'G', 'GGC' => 'G',
			 'GGG' => 'G', 'GGT' => 'G', 'GTA' => 'V',
			 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',
			 'TAA' => '*', 'TAC' => 'Y', 'TAG' => '*',
			 'TAT' => 'Y', 'TCA' => 'S', 'TCC' => 'S',
			 'TCG' => 'S', 'TCT' => 'S', 'TGA' => '*',
			 'TGC' => 'C', 'TGG' => 'W', 'TGT' => 'C',
			 'TTA' => 'L', 'TTC' => 'F', 'TTG' => 'L',
			 'TTT' => 'F' );
			 
my %kyte = ( 'A' => 1.80, 'C' => 2.50, 'D' => -3.5, 'E' => 3.50, 
			 'F' => 2.80, 'G' => -.40, 'H' => -3.2, 'I' => 4.50,
			 'K' => -3.9, 'L' => 3.80, 'M' => 1.90, 'N' => -3.5,
			 'P' => -1.6, 'Q' => -3.5, 'R' => -4.5, 'S' => -.80,
			 'T' => -.70, 'V' => 4.20, 'W' => -.90, 'Y' => -1.3,);
			 
my %hopp = ( 'A' => -.50, 'C' => -1.0, 'D' => 3.00, 'E' => 3.00, 
			 'F' => -2.5, 'G' => 0.00, 'H' => -.50, 'I' => -1.8,
			 'K' => 3.00, 'L' => -1.8, 'M' => -1.3, 'N' => 0.20,
			 'P' => 0.00, 'Q' => 0.20, 'R' => 3.00, 'S' => 0.30,
			 'T' => -.40, 'V' => -1.5, 'W' => -3.4, 'Y' => -2.3,);
			 

sub revcomp {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr[ACGTRYMKWSBDHVacgtrymkwsbdhv]			  [TGCAYRKMWSVHBDtgcayrkmwsvhdb];
	return $seq;
}


sub gccomp {
	die "Usage: <seq>, <size>, <step>\n" unless @_ == 3;
	my ($seq, $size, $step) = @_;
	my $position = 0;
	my @gc;
	my $current = substr($seq, 0, $size);
	$gc[$position] = ($current =~ tr/GC/GC/)/length($current);
	for (my $i = $size; $i<=length($seq); $i+=$step) {
		$position++;
		my $correction=0;
		if ($size>$step) { # sliding window w/ more efficient overlap
			$current = substr($current, $step, length($current)-$step);
		}
		else {
			$current =""; #if no overlap
			$correction = $step-$size; #if step > size, corrects position
		}
		$current .= substr($seq, $i+$correction, $size-length($current));
		$gc[$position] = ($current =~ tr/GC/GC/)/length($current);
	}
	my $printgc = join (', ', @gc);
	print "GC Content: $printgc\n";
	#return @gc;
}


sub hydrophobicity {
	die "Usage: <seq>, <size>, <step> <scale>\n"
		unless (@_==3);
	my ($seq, $size, $step) = @_;
	my (@kyte, @hopp);
	my $position = 0;
	#die "Input sequence is not valid\n" unless 
	my $current = substr($seq, 0, $size);
	for (my $j=0; $j<length($current);$j++){
			my $aa = substr($current, $j, 1);
			if (exists $kyte{$aa}) {
				$kyte[$position] += $kyte{$aa};
				$hopp[$position] += $hopp{$aa};
			}
		}
	
	
	for (my $i=$size; $i<=length($seq); $i+=$step) {
		$position++;
		my $correction=0;
		if ($size>$step) { # sliding window w/ more efficient overlap
			$current = substr($current, $step, length($current)-$step);
		}
		else {
			$current =""; #if no overlap
			$correction = $step-$size; #if step > size, corrects position
		}
		$current .= substr($seq, $i+$correction, $size-length($current));

		for (my $j=0; $j<length($current);$j++){
			my $aa = substr($current, $j, 1);
			if (exists $kyte{$aa}) { #total 
				$kyte[$position] += $kyte{$aa};
				$hopp[$position] += $hopp{$aa};
			}
		}
	}
	my $printkyte = join (', ', @kyte);
	print "Kyte-Doolittle Scale: $printkyte\n";
	my $printhopp = join (', ', @hopp);
	print "Hopp-Woods Scale: $printhopp\n";
}
sub isdna {
	die "Usage: <seq>\n" unless @_==1;
	my ($seq) = @_;
	my $dna = 1;
	if ($seq =~ m/[^ACGTRYMKWSBDHV]/i){
		$dna = 0
	}
	return $dna;
}

sub translate {
	die "Usage: <seq> <frame>\n" unless @_==2;
	my ($seq, $frame) = @_;
	die unless isdna($seq);
	if (abs($frame)>3||$frame==0){die "Inappropriate frame argument.\n";}
	my $prot;
	if ($frame<0){
		$seq=revcomp($seq);
		$frame = abs($frame);
	}
	for (my $i = $frame-1; $i < length($seq); $i +=3) {
		my $codon = uc substr($seq, $i, 3);
		if (exists $code{$codon}) {$prot .= $code{$codon}}
		else					  {$prot .= 'X'}
	}
	return $prot;
}

sub readfasta {
	die "Usage: $0 <fasta file>\n" unless @_ == 1;
	open (my $fh, $_[0]) or die;
	my $header = <$fh>;
	my $seq;
	while(<$fh>) {
		chomp;
		$seq .= $_;
	}
	close $fh;
	$seq = uc $seq;
	return $seq;
}

1;
