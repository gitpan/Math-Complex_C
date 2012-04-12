package Math::Complex_C::Long;
use strict;

require Exporter;
*import = \&Exporter::import;
require DynaLoader;

use overload
    '**'    => \&_overload_pow,
    '*'     => \&_overload_mul,
    '+'     => \&_overload_add,
    '/'     => \&_overload_div,
    '-'     => \&_overload_sub,
    '**='   => \&_overload_pow_eq,
    '*='    => \&_overload_mul_eq,
    '+='    => \&_overload_add_eq,
    '/='    => \&_overload_div_eq,
    '-='    => \&_overload_sub_eq,
    'sqrt'  => \&_overload_sqrt,
    '=='    => \&_overload_equiv,
    '!='    => \&_overload_not_equiv,
    '!'     => \&_overload_not,
    'bool'  => \&_overload_true,
    '='     => \&_overload_copy,
    '""'    => \&_overload_string,
    'abs'   => \&_overload_abs,
    'exp'   => \&_overload_exp,
    'log'   => \&_overload_log,
    'sin'   => \&_overload_sin,
    'cos'   => \&_overload_cos,
    'atan2' => \&_overload_atan2;

$Math::Complex_C::Long::VERSION = '0.05';

DynaLoader::bootstrap Math::Complex_C::Long $Math::Complex_C::Long::VERSION;

@Math::Complex_C::Long::EXPORT = ();
@Math::Complex_C::Long::EXPORT_OK = ();

sub dl_load_flags {0} # Prevent DynaLoader from complaining and croaking

sub _overload_string {
     return "(" . real_cl($_[0]) . " " . imag_cl($_[0]) . ")";
}

sub new {

    my $ret = create_cl();

    # This function caters for 2 possibilities:
    # 1) that 'new' has been called OOP style - in which 
    #    case there will be a maximum of 3 args
    # 2) that 'new' has been called as a function - in
    #    which case there will be a maximum of 2 args.
    # If there are no args, then we just want to return a
    # Math::Complex_C::Long object

    if(!@_) {return $ret}
   
    if(@_ > 3) {die "Too many arguments supplied to new()"}

    # If 'new' has been called OOP style, the first arg is the string
    # "Math::Complex_C::Long" which we don't need - so let's remove it.

    if($_[0] eq "Math::Complex_C::Long") {
      shift;
      if(!@_) {return $ret}
    }

    if(@_ > 2) {die "Bad argument list supplied to new()"} 

    if(@_ == 2) {assign_cl($ret, $_[0], $_[1])}
    else {assign_cl($ret, $_[0], 0)}

    return $ret;    
}


1;