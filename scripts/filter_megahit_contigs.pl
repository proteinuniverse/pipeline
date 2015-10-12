use strict;
use Getopt::Long; 
use Pod::Usage;

my $PRINT = 0;
my $coverage = 0;
my $length = 0;
my ($help, $in, $id);

GetOptions(
        'h'     => \$help,
        'c=f'   => \$coverage,
        'i=s'   => \$in,
	'l=i'   => \$length,
) or pod2usage(0);

pod2usage(-exitstatus => 0,
          -output => \*STDOUT,
          -verbose => 1,
          -noperldoc => 1,
         ) if $help;

pod2usage(-exitstatus => 0,
          -output => \*STDOUT,
          -verbose => 1,
          -noperldoc => 1,
         ) unless ($in);

open IN, $in or die "cannot open $in";

while(<IN>){
  my $cov = 0;
  my $id;
  my $len = 0;

  if (/>/) {
    # look at the id
    $id = $1 if />(\S+)/;
    $id = $1 if />(\S+)_length/;
    die "could not parse id" unless $id;

    # look at coverage
    $cov = $1 if /multi=([\d\.]+)/;
    $cov = $1 if /multi_([\d\.]+)/;
    die "could not parse coverage" unless $cov;

    # look at length
    $len = $1 if /len=(\d+)/;
    $len = $1 if /length_(\d+)/;
    die "could not parse length" unless $len;

    if ($cov >= $coverage && $len >= $length) {
      $PRINT = 1;
    }
    else {
      $PRINT = 0;
    }
  }

 
  if ($PRINT == 1) {
    print;
    print STDERR "$id\t$len\t$cov\n" if />/;
  }
} 


=pod

=head1  NAME

[% kb_method %]

=head1  SYNOPSIS

[% kb_method %] <options>

=head1  DESCRIPTION

The [% kb_method %] command calls the [% kb_method %] method of a [% kb_client %] object.

=head1  OPTIONS

=over

=item   -h, --help

Basic usage documentation

=item   -c

Coverage threshold for inclusion into output.

=item	-l

Minimum length for incusion into the output.

=item   -i

The input megahit assembly file, often named final.contigs.fa

=item   -o

The output file that will contain the filtered contigs

=back

=head1  AUTHORS

[% kb_author %]

=cut

1;
