use strict;

for (my $x=0; $x< 10; $x+=2) {
  download_from_mgrast(2, $x);
}


sub download_from_mgrast {
  my ($limit, $offset) = @_ or die "must provide limit and offset";
  my $cmd = 'curl \'http://api.metagenomics.anl.gov/1/metagenome?limit=';
  $cmd   .= $limit;
  $cmd   .= '&order=name&verbosity=metadata&sequence_type=WGS&offset=';
  $cmd   .= $offset;
  $cmd   .= '\'';
  print "running $cmd\n";

  system($cmd)

}

