#!/usr/bin/perl
# Program:     extract_sequences_list
# Programmer:  Sean R. McCorkle
#              Biological, Environmental and Climate Sciences Dept.
#              Brookhaven National Laboratory
# Language:    perl
# 
# Description: Reads a list of identifiers from a file, one per line.
#              Then passes any fasta sequences whose headers start with
#              any of those ids (assuming the id is the first whitespace
#              separated entry on the header).
#              
# Usage:       extract_sequences_list  [-u]  <wanted_list>  [seqs...]
#
# Options:      -u  unique - don't repeatedly pass the same sequence if it
#                   recurs 
#               -v  complementary sense - DON'T pass the sequence if its 
#                   found in the list
#
use strict;

our ( $opt_u, $opt_v );

use Getopt::Std;

getopts( "uv" ) || die "Bad option: only -v and -u recognized\n";

my %want = get_list( (shift) );
my %have;

my $pass = 0;

while ( <> )
   {
    next if ( /^\s*$/ );
    if ( /^>(\S+)+/ )
       { 
        $pass = $want{$1}; 
        $pass = ! $pass if ( $opt_v );
        $pass = 0 if ( $opt_u && $pass && $have{$1} );
        #undef( $want{$1} ) if ( $want{$1} && $opt_u );  # clear after first time
        $have{$1} = 1 if ( $opt_u && $pass );
       }
    print if ( $pass );
   }


sub  get_list
   {
    my $file = shift;
    my %w = ();

    open( LIST, $file ) || die "Can't open $file: $!\n";
    while ( $_ = <LIST> )
       {
        chomp;
        s/^\s*//;
        s/\s*$//;
        $w{$_}++;
       }
    close( LIST );

    return( %w );
   }
