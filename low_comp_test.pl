#!/usr/bin/perl
use strict; use warnings;
use Low_Complexity;

my $seq = "acgt";
$seq = Low_Complexity::low_comp($seq, 5, 1);
print "$seq\n";