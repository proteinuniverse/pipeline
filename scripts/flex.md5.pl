#!/usr/local//bin/perl -w

#  Example/idea  taken from http://perldoc.perl.org/Digest/MD5.html
use Digest::MD5 qw(md5_hex);
use Digest::MD5 qw(md5_base64);

use Math::Fleximal;
  
# Converts a hex representation of a number into
# one that uses more alphanumerics.  (ie base 62)
sub hex2alpha {
  Math::Fleximal->new(
      lc(shift), [0..9, 'a'..'f']
    )->change_flex(
      [0..9,'a'..'z','A'..'Z']
    )->to_str();
}

use strict;

my $file = $ARGV[0] || '';

# Assumptions:
#   1. We have a fasta file where the first line is the first fasta record (no whitespace or comments)
#   2. We have a fasta file where there is a kbase name, no fasta comment, and
#      the sequence is all on on line. For example:

#>kb|g.1031.peg.0
#MIISVFSPKGGVGKTTVALALAESLSKNHRVVALELDFSPGDFVGLLHELDPGKNLLTCKHDILSAVQRPSGKEFDVIIGGYPGEHEHVRREDIKRCIEILKFKYEYIIVDIQPGIVELVIDVLAESDRVLVIAEENFITPIARINAFLDWIQINNLSDLKNFVFVRNKVTNKELVYIDKIKHSLKLVHDIPFYKKLKGYDDKRLQKNIKRLAGVLRNGTVREDKRFWLFRRILGKL
#>kb|g.1031.peg.1
#MKVKYTLSVLVENHPGVLSRVAGLFSRRGFNIDSLAVGVTEDPTISRMTIVVNGDDYIVEQVTKQLNKLIDVIKVKKLNPKEAVERELALIKVNANSQTRSDIIQITEIFRANIVDVSKETLTIEISGDEDKIEALIELLKQYGIREVVRTGLIAIERGNKVISKSKSEEDD

open(MD5, $file) || die("\nERROR: $0, Can't open file, $file: $!\n\n");
my $seq_name  =  '';
my $seq_md5   =  '';
my $seq       =  '';

while(<MD5>)
{
  chomp();
  if (/^>/)
  {
    chomp;
    if($seq_name && $seq) 
    {
      $seq_name   =~ s/^\s*\>\s+|\s+$//;
      $seq =~ s/\*$//;
      $seq_md5  = hex2alpha(md5_hex($seq));
      print("$seq_md5\t$seq_name\t$seq\n");   #  Tab seperated with MD5, Kbase name, Sequence
    }
    $seq_name = $_;
    $seq='';
  }
  else
  {
    $seq .= $_;
  }
}

if($seq_name && $seq)
{
  # get rid of trailing * and leading/trailing whitespace
  $seq =~ s/\*$//;
  $seq_name   =~ s/^\s*\>\s+|\s+$//;
  $seq_md5  = hex2alpha(md5_hex($seq));
  print("$seq_md5\t$seq_name\t$seq\n");   #  Tab seperated with MD5, Kbase name, Sequence
}

close(MD5);
