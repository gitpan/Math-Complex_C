use warnings;
use strict;
use Math::Complex_C qw(:all);

print "1..2\n";

my $c_type = Math::Complex_C::Long::_complex_type();

if($c_type eq 'long double _Complex') {print "ok 1\n"}
else {
  warn "\$c_type: $c_type\n";
  print "not ok 1\n";
}

my $d_type = Math::Complex_C::Long::_double_type();

if($d_type eq "long double") {print "ok 2\n"}
else {
  warn "\$d_type: $d_type\n";
  print "not ok 2\n";
}