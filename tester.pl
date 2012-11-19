#!/usr/bin/perl
use strict; use warnings;
use FAlite;
my ($Fasta_file) = @ARGV;

die "usage:$0 fasta_file" unless @ARGV ==1;
open (my $fasta_file, "<", "$file") or die "Couldn't open $file\n";

my $fasta = new FALite($fasta_file);

while (my $entry = $fasta->nextEntry) {
	my $header = $entry->def;
	my $seq=$entry->seq;
	
	print "$header\n";
	print "$seq\n\n";
}

close($fasta_file);