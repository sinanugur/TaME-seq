package MySam;

use strict;
use warnings;
use Carp;

# sam columns
our $SAM_NAME = 0;
our $SAM_FLAG = 1;
our $SAM_RNAME = 2;
our $SAM_POS = 3;
our $SAM_MAPQ = 4;
our $SAM_CIGAR = 5;
our $SAM_MRNM = 6;
our $SAM_MPOS = 7;
our $SAM_IPOS = 8;
our $SAM_SEQ = 9;
our $SAM_QUAL = 10;
our $SAM_TAG = 11;
our $SAM_ISIZE = 12;
our $SAM_VTYPE = 13;
our $SAM_VALUE = 14;

our $HEADER = '';

#
# perl package to handle sam files
#

sub _get_line(\*) {
  my $fh = shift @_;
  my $line = <$fh>;
  !defined($line) && return ();
  chomp($line);
  return (split(/\t/,$line,-1));
}

sub get_read(\*\@) {
  my $fh = shift @_;
  my $ar = shift @_;
  local $_;
  @{$ar} = _get_line(*{$fh});
  return scalar(@{$ar});
}
sub get_read_pair(\*\@;$) {
  my $fh = shift @_;
  my $ar = shift @_;
  my $verbose = shift @_ || 0;
  local $_;
  my @skip = ();
  @{$ar} = ();
  foreach (0 .. 1) {
    @{$ar->[$_]} = _get_line(*{$fh});
    scalar(@{$ar->[$_]}) || return 0;
  }
  while ($ar->[0][$SAM_NAME] ne $ar->[1][$SAM_NAME]) {
    # die "Read pairing failed: $ar->[0][$SAM_NAME] differs from ",
    #   "$ar->[1][$SAM_NAME]\n";
    $verbose && push (@skip, $ar->[0][$SAM_NAME]);
    @{$ar->[0]} = @{$ar->[1]};
    @{$ar->[1]} = _get_line(*{$fh});
    if (!scalar(@{$ar->[1]})) {
      return 0;
    }
  }
  if (scalar(@skip)) {
    print STDERR "Skipped reads for pair missing: ", join(", ",@skip), "\n";
  }
  return scalar(@{$ar});
}

sub skip_header(\*) {
  my $fh = shift @_;
  local $_;
  my $place = tell($fh);
  while (<$fh>) {
    m/^@/ || last;
    $place = tell($fh);
  }
  seek($fh,$place,0);
  return;
}

# flag and their values
#   0 => "the read is paired in sequencing, no matter whether it is mapped in a pair",
#   1 => "the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment)",
#   2 => "the query sequence itself is unmapped",
#   3 => "the mate is unmapped",
#   4 => "reverse strand", # "strand of the query (0 for forward; 1 for reverse strand)",
#   5 => "strand of the mate",
#   6 => "1st read", # the read is the ﬁrst read in a pair",
#   7 => "2nd read", # "the read is the second read in a pair",
#   8 => "the alignment is not primary (a read having split hits may have multiple primary alignment records)",
#   9 => "the read fails platform/vendor quality checks",
#   10 => "the read is either a PCR duplicate or an optical duplicate",
sub flag_check($@) {
  my $bit = shift @_;
  my @flags = @_;
  local $_;
  my $ret = 0;
  foreach (@flags) {
    if ($_ & 2**$bit) {
      $ret = 1;
      last;
    }
  }
  return $ret;
}

# returns true if flag is asked
sub flag_is_aligned(@) {
  return !flag_check(2, @_);
}
sub flag_is_reverse(@) {
  return flag_check(4, @_);
}
sub flag_is_mate_aligned {
  return !flag_check(3, @_);
}
sub flag_is_first {
  return flag_check(6, @_);
}

sub zero_nonalign (\$\$$) {
  my $pos_r = shift @_;
  my $name_r = shift @_;
  my $flag = shift @_;
  my $ret = 0;
  if (!flag_is_aligned($flag)) {
    $$pos_r = 0;
    $$name_r = 'undef';
    $ret = 1;
  }
  return $ret;
}

sub remove_eq_from_chr (\@) {
  my $row = shift @_;
  $row->[$SAM_RNAME] =
    $row->[$SAM_RNAME] eq '=' ? $row->[$SAM_MRNM] : $row->[$SAM_RNAME];
  $row->[$SAM_MRNM] =
    $row->[$SAM_MRNM] eq '=' ? $row->[$SAM_RNAME] : $row->[$SAM_MRNM];
  return;
}

sub last_pos {
  my @read = @_;
  local $_;
  my ($match,$del) = (0,0);
  while ($read[$SAM_CIGAR] =~ m/(\d+)M/g) {
    $match += $1;
  }
  while ($read[$SAM_CIGAR] =~ m/(\d+)D/g) {
    $del += $1;
  }
  return $read[$SAM_POS] + $match + $del -1;
}

1;
