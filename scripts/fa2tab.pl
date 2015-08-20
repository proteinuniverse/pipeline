#! /usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = "Usage: $0 < contigs.tab > contigs.fa\n\n";

my ($help, $min_contig_length);

GetOptions("h|help"  => \$help,
           "m|min=i" => \$min_contig_length);

$help and die $usage;  

while (<STDIN>) {
    my $hdr = $_;
    my $seq = <STDIN>;
    
    if ($hdr = /^>(\S+)\s+len_(\d+)_cov_([0-9.]+)_stdev_([0-9.]+)_GC_([0-9.]+)_seed_(\S+)/ && $seq) {
        print join("\t", $1, $2, $3, $4, $5, $6, $seq) if $2 >= $min_contig_length;
    }
    # megahit defline: >contig_21_1366_length_231_multi_4_in_0_out_0
    elsif ($hdr = /^>(contig_\d+_\d+)_length_(\d+)/ && $seq) {
	my @a = split//, uc($seq);
	my $gc = 0;
	foreach (@a) { $gc++ if /[GC]/ };
	my $gc_pct = int(($gc / int(@a)) * 100);
	print join("\t", $1, $2, "-", "-", $gc_pct, "-", $seq) if $2 >= $min_contig_length;
    }
    # if all else fails, compute length and gc_pct
    elsif ($hdr = /^>(\S+)/ && $seq) {
	my @a = split//, uc($seq);
        my $gc = 0;
        foreach (@a) { $gc++ if /[GC]/ };
        my $gc_pct = int(($gc / int(@a)) * 100);
        print join("\t", $1, length($seq), "-", "-", $gc_pct, "-", $seq) if length($seq)  >= $min_contig_length;
    }
    else {
	die "did not match header";
    }
}
