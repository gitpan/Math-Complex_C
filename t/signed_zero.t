use warnings;
use strict;

print "1..2\n";

use Math::Complex_C qw(:all);

my $zero = Math::Complex_C->new(0.0, 0.0);
my $rop = Math::Complex_C->new();

mul_c_nv($rop, $zero, -1.1);

my $re = real_c($rop);
my $im = imag_c($rop);

my $correct = 0;

if("$re" eq "-0") {
  $correct = 1;
  print "ok 1\n";
}
else {
  warn "\n  \$re: $re (correct result is '-0')\n";
  print "not ok 1\n";
}

if("$im" eq "-0") {
  print "ok 2\n";
}
else {
  if($correct) {
    warn "\n  \$im: $im\n";
    warn "Correct result is '-0' - this is a known issue with some (old) compilers\n";
    print "not ok 2\n";
  }
  else {
    warn "\n  \$im: $im (correct result is '-0')\n";
    print "not ok 2\n";
  }
}

