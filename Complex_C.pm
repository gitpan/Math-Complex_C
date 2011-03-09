## This file generated by InlineX::C2XS (version 0.17) using Inline::C (version 0.48)
package Math::Complex_C;

require Exporter;
*import = \&Exporter::import;
require DynaLoader;

use Math::Complex_C::Long; # 'long double _Complex' support.

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
    'not'   => \&_overload_not,
    'bool'  => \&_overload_true,
    '='     => \&_overload_copy,
    '""'    => \&_overload_string,
    'abs'   => \&_overload_abs,
    'exp'   => \&_overload_exp,
    'log'   => \&_overload_log,
    'sin'   => \&_overload_sin,
    'cos'   => \&_overload_cos,
    'atan2' => \&_overload_atan2;

$Math::Complex_C::VERSION = '0.03';

DynaLoader::bootstrap Math::Complex_C $Math::Complex_C::VERSION;

@Math::Complex_C::EXPORT = ();
@Math::Complex_C::EXPORT_OK = qw(
    create_c assign_c mul_c mul_c_nv mul_c_iv mul_c_uv div_c div_c_nv div_c_iv div_c_uv add_c 
    add_c_nv add_c_iv add_c_uv sub_c sub_c_nv sub_c_iv sub_c_uv real_c 
    imag_c arg_c abs_c conj_c acos_c asin_c atan_c cos_c sin_c tan_c acosh_c asinh_c atanh_c 
    cosh_c sinh_c tanh_c exp_c log_c sqrt_c proj_c pow_c get_nan get_inf is_nan is_inf

    create_cl assign_cl mul_cl mul_c_nvl mul_c_ivl mul_c_uvl div_cl div_c_nvl div_c_ivl div_c_uvl add_cl 
    add_c_nvl add_c_ivl add_c_uvl sub_cl sub_c_nvl sub_c_ivl sub_c_uvl real_cl 
    imag_cl arg_cl abs_cl conj_cl acos_cl asin_cl atan_cl cos_cl sin_cl tan_cl acosh_cl asinh_cl atanh_cl
    cosh_cl sinh_cl tanh_cl exp_cl log_cl sqrt_cl proj_cl pow_cl get_nanl get_infl is_nanl is_infl
    );

%Math::Complex_C::EXPORT_TAGS = (all => [qw(
    create_c assign_c mul_c mul_c_nv mul_c_iv mul_c_uv div_c div_c_nv div_c_iv div_c_uv add_c 
    add_c_nv add_c_iv add_c_uv sub_c sub_c_nv sub_c_iv sub_c_uv real_c 
    imag_c arg_c abs_c conj_c acos_c asin_c atan_c cos_c sin_c tan_c acosh_c asinh_c atanh_c 
    cosh_c sinh_c tanh_c exp_c log_c sqrt_c proj_c pow_c get_nan get_inf is_nan is_inf

    create_cl assign_cl mul_cl mul_c_nvl mul_c_ivl mul_c_uvl div_cl div_c_nvl div_c_ivl div_c_uvl add_cl 
    add_c_nvl add_c_ivl add_c_uvl sub_cl sub_c_nvl sub_c_ivl sub_c_uvl real_cl 
    imag_cl arg_cl abs_cl conj_cl acos_cl asin_cl atan_cl cos_cl sin_cl tan_cl acosh_cl asinh_cl atanh_cl
    cosh_cl sinh_cl tanh_cl exp_cl log_cl sqrt_cl proj_cl pow_cl get_nanl get_infl is_nanl is_infl
    )]);

*create_cl = \&Math::Complex_C::Long::create_cl;
*assign_cl = \&Math::Complex_C::Long::assign_cl;
*mul_cl    = \&Math::Complex_C::Long::mul_cl;
*mul_c_nvl = \&Math::Complex_C::Long::mul_c_nvl;
*mul_c_ivl = \&Math::Complex_C::Long::mul_c_ivl;
*mul_c_uvl = \&Math::Complex_C::Long::mul_c_uvl ;
*div_cl    = \&Math::Complex_C::Long::div_cl;
*div_c_nvl = \&Math::Complex_C::Long::div_c_nvl;
*div_c_ivl = \&Math::Complex_C::Long::div_c_ivl;
*div_c_uvl = \&Math::Complex_C::Long::div_c_uvl;
*add_cl    = \&Math::Complex_C::Long::add_cl;
*add_c_nvl = \&Math::Complex_C::Long::add_c_nvl;
*add_c_ivl = \&Math::Complex_C::Long::add_c_ivl;
*add_c_uvl = \&Math::Complex_C::Long::add_c_uvl ;
*sub_cl    = \&Math::Complex_C::Long::sub_cl;
*sub_c_nvl = \&Math::Complex_C::Long::sub_c_nvl;
*sub_c_ivl = \&Math::Complex_C::Long::sub_c_ivl;
*sub_c_uvl = \&Math::Complex_C::Long::sub_c_uvl;
*real_cl   = \&Math::Complex_C::Long::real_cl;
*imag_cl   = \&Math::Complex_C::Long::imag_cl;
*arg_cl    = \&Math::Complex_C::Long::arg_cl;
*abs_cl    = \&Math::Complex_C::Long::abs_cl;
*conj_cl   = \&Math::Complex_C::Long::conj_cl;
*acos_cl   = \&Math::Complex_C::Long::acos_cl;
*asin_cl   = \&Math::Complex_C::Long::asin_cl;
*atan_cl   = \&Math::Complex_C::Long::atan_cl;
*cos_cl    = \&Math::Complex_C::Long::cos_cl;
*sin_cl    = \&Math::Complex_C::Long::sin_cl;
*tan_cl    = \&Math::Complex_C::Long::tan_cl;
*acosh_cl  = \&Math::Complex_C::Long::acosh_cl;
*asinh_cl  = \&Math::Complex_C::Long::asinh_cl;
*atanh_cl  = \&Math::Complex_C::Long::atanh_cl;
*cosh_cl   = \&Math::Complex_C::Long::cosh_cl;
*sinh_cl   = \&Math::Complex_C::Long::sinh_cl;
*tanh_cl   = \&Math::Complex_C::Long::tanh_cl;
*exp_cl    = \&Math::Complex_C::Long::exp_cl;
*log_cl    = \&Math::Complex_C::Long::log_cl;
*sqrt_cl   = \&Math::Complex_C::Long::sqrt_cl;
*proj_cl   = \&Math::Complex_C::Long::proj_cl;
*pow_cl    = \&Math::Complex_C::Long::pow_cl;
*get_nanl  = \&Math::Complex_C::Long::get_nanl;
*get_infl  = \&Math::Complex_C::Long::get_infl;
*is_nanl   = \&Math::Complex_C::Long::is_nanl;
*is_infl   = \&Math::Complex_C::Long::is_infl;

sub dl_load_flags {0} # Prevent DynaLoader from complaining and croaking

sub _overload_string {
     return "(" . real_c($_[0]) . " " . imag_c($_[0]) . ")";
}

sub new {

    my $ret = create_c();

    # This function caters for 2 possibilities:
    # 1) that 'new' has been called OOP style - in which 
    #    case there will be a maximum of 3 args
    # 2) that 'new' has been called as a function - in
    #    which case there will be a maximum of 2 args.
    # If there are no args, then we just want to return a
    # Math::Complex_C object

    if(!@_) {return $ret}
   
    if(@_ > 3) {die "Too many arguments supplied to new()"}

    # If 'new' has been called OOP style, the first arg is the string
    # "Math::Complex_C" which we don't need - so let's remove it.

    if($_[0] eq "Math::Complex_C") {
      shift;
      if(!@_) {return $ret}
    }

    if(@_ > 2) {die "Bad argument list supplied to new()"} 

    if(@_ == 2) {assign_c($ret, $_[0], $_[1])}
    else {assign_c($ret, $_[0], 0)}

    return $ret;    
}

1;

__END__

=head1 NAME

Math::Complex_C - perl interface to C's complex.h functions.

=head1 DEPENDENCIES

   In order to compile this module, a C compiler that provides
   complex.h is needed.

=head1 DESCRIPTION
   This module wraps the 'double _Complex' type (as a Math::Complex_C
   object) and the 'long double _Complex' type (as a Math::Complex_C::Long
   object).


   use warnings;
   use strict;
   use Math::Complex_C qw(:all);
   my $c =    Math::Complex_C->new(12.5, 1125); # 'double _Complex' type
   my $root = Math::Complex_C->new();
   my $cl = Math::Complex_C::Long->new(5.9,1.1); # 'long double _Complex'
   my $rootl = Math::Complex_C::Long->new();
   sqrt_c($root, $c);
   sqrt_cl($rootl, $cl);
   print "Square root of $c is $root\n";
   print "Square root of $c1 is $rootl\n";

   On many perls, the values printed out by the above code will be
   identical - see the README for some elaboration.

   Note that Math::Complex_C and Math::Complex_C::Long objects use
   different functions (sqrt_c vs sqrt_cl, in the above example).

   See also the Math::Complex_C test suite for some (simplistic) examples
   of usage.

=head1 FUNCTIONS

   $rop = Math::Complex_C->new([$re, [$im]]);
   $rop = Math::Complex_C::new([$re, [$im]]);
   $rop = new Math::Complex_C([$re, [$im]]);
   $ropl = Math::Complex_C::Long->new([$re, [$im]]);
   $ropl = Math::Complex_C::Long::new([$re, [$im]]);
   $ropl = new Math::Complex_C::Long([$re, [$im]]);
    $rop/$ropl is a Math::Complex_C/Math::Complex_C::Long object;
    $re and $im are the real and imaginary values (respectively) that $rop
    holds. They can be an integer (signed or unsigned) or a floating point
    value.Integer values (IV/UV) will be converted to floating point (NV)
    before being assigned.
    Note that the two arguments ($re $im) are optional - ie can be omitted.
    If no arguments are supplied to new, then $rop will be assigned NaN for
    both the real and imaginary parts.
    If only one argument is supplied, then $rop will be assigned that value
    for the real part, and 0 for the imaginary part.
    
   $rop = create_c();
   $ropl = create_cl();
    $rop/$ropl is a Math::Complex_C/Math::Complex_C::Long object (created
    with both real and imaginary values set to NaN).

   assign_c($rop, $re, $im);
   assign_cl($ropl, $re, $im); 
    The real part of $rop/$ropl is set to the value of $re, the imaginary
    part is set to the value of $im.
    $re and $im can be an integer (signed or unsigned) or a floating point
    value. Integer values (IV/UV) will be converted to floating point (NV)
    before being assigned.

   mul_c($rop, $op1, $op2);
   mul_c_iv($rop, $op1, $si);
   mul_c_uv($rop, $op1, $ui);
   mul_c_nv($rop, $op1, $nv);
   mul_cl($ropl, $op1, $op2);
   mul_c_ivl($ropl, $op1, $si);
   mul_c_uvl($ropl, $op1, $ui);
   mul_c_nvl($ropl, $op1, $nv);
    Multiply $op1 by the 3rd arg, and store the result in $rop/$ropl.
    The "3rd arg" is (respectively, from the top) a Math::Complex_C object,
    a signed integer value (IV), an unsigned integer value (UV), and a
    floating point value (NV).

   add_c($rop, $op1, $op2);
   add_c_iv($rop, $op1, $si);
   add_c_uv($rop, $op1, $ui);
   add_c_nv($rop, $op1, $nv);
   add_cl($ropl, $op1, $op2);
   add_c_ivl($ropl, $op1, $si);
   add_c_uvl($ropl, $op1, $ui);
   add_c_nvl($ropl, $op1, $nv);
    As for mul_c(), etc., but performs addition.

   div_c($rop, $op1, $op2);
   div_c_iv($rop, $op1, $si);
   div_c_uv($rop, $op1, $ui);
   div_c_nv($rop, $op1, $nv);
   div_cl($ropl, $op1, $op2);
   div_c_ivl($ropl, $op1, $si);
   div_c_uvl($ropl, $op1, $ui);
   div_c_nvl($ropl, $op1, $nv);
    As for mul_c(), etc., but performs division.

   sub_c($rop, $op1, $op2);
   sub_c_iv($rop, $op1, $si);
   sub_c_uv($rop, $op1, $ui);
   sub_c_nv($rop, $op1, $nv);
   sub_cl($ropl, $op1, $op2);
   sub_c_ivl($ropl, $op1, $si);
   sub_c_uvl($ropl, $op1, $ui);
   sub_c_nvl($ropl, $op1, $nv);
    As for mul_c(), etc., but performs subtraction.

   $nv = real_c($op);
   $nv = real_cl($op);
    Returns the (floating point) value of the real component of $op.
    Wraps C's 'creal/creall' function.

   $nv = imag_c($op);
   $nv = imag_cl($op);
    Returns the (floating point) value of the imaginary component of $op.
    Wraps C's 'cimag/cimagl' function.

   $nv = arg_c($op);
   $nv = arg_cl($op);
    Returns the (floating point) value of the argument of $op.
    Wraps C's 'carg/cargl' function.

   $nv = abs_c($op);
   $nv = abs_cl($op);
    Returns the (floating point) absolute value of $op.
    Wraps C's 'cabs/cabsl' function.

   conj_c($rop, $op);
   conj_cl($ropl, $op);
    Sets $rop/$ropl to the conjugate of $op.
    Wraps C's 'conj/conjl' function.

   acos_c($rop, $op);
   acos_cl($ropl, $op);
    Sets $rop/$ropl to acos($op).
    Wraps C's 'cacos/cacosl' function.

   asin_c($rop, $op);
   asin_cl($ropl, $op);
    Sets $rop/$ropl to asin($op).
    Wraps C's 'casin/casinl' function.

   atan_c($rop, $op);
   atan_cl($ropl, $op);
    Sets $rop/$ropl to atan($op).
    Wraps C's 'catan/catanl' function.

   cos_c($rop, $op);
   cos_cl($ropl, $op);
    Sets $rop/$ropl to cos($op).
    Wraps C's 'ccos/ccosl' function.

   sin_c($rop, $op);
   sin_cl($ropl, $op);
    Sets $rop/$ropl to sin($op).
    Wraps C's 'csin/csinl' function.

   tan_c($rop, $op);
   tan_cl($ropl, $op);
    Sets $rop/$ropl to tan($op).
    Wraps C's 'ctan/ctanl' function.

   acosh_c($rop, $op);
   acosh_cl($ropl, $op);
    Sets $rop/$ropl to acosh($op).
    Wraps C's 'cacosh/cacoshl' function.

   asinh_c($rop, $op);
   asinh_cl($ropl, $op);
    Sets $rop/$ropl to asinh($op).
    Wraps C's 'casinh/casinhl' function.

   atanh_c($rop, $op);
   atanh_cl($ropl, $op);
    Sets $rop/$ropl to atanh($op).
    Wraps C's 'catanh/catanhl' function.

   cosh_c($rop, $op);
   cosh_cl($ropl, $op);
    Sets $rop/$ropl to cosh($op).
    Wraps C's 'ccosh/ccoshl' function.

   sinh_c($rop, $op);
   sinh_cl($rop, $op);
    Sets $rop/$ropl to sinh($op).
    Wraps C's 'csinh/csinhl' function.

   tanh_c($rop, $op);
   tanh_cl($ropl, $op);
    Sets $rop/$ropl to tanh($op).
    Wraps C's 'ctanh/ctanhl' function.

   exp_c($rop, $op);
   exp_cl($rop1, $op);
    Sets $rop/$ropl to e ** $op.
    Wraps C's 'cexp/cexpl' function.

   log_c($rop, $op);
   log_cl($ropl, $op);
    Sets $rop/$ropl to log($op).
    Wraps C's 'clog/clogl' function.

   pow_c($rop, $op1, $op2);
   pow_cl($ropl, $op1, $op2);
    Sets $rop/$ropl to $op1 ** $op2.
    Wraps C's 'cpow/cpowl' function.

   sqrt_c($rop, $op);
   sqrt_cl($ropl, $op);
    Sets $rop/$ropl to sqrt($op).
    Wraps C's 'csqrt/csqrtl' function.

   proj_c($rop, $op);
   proj_cl($ropl, $op);
    Sets $rop/$ropl to a projection of $op onto the Riemann sphere.
    Wraps C's 'cproj/cprojl' function.

   $rop = get_nan();
   $ropl = get_nanl();
    Sets $rop/$ropl to NaN.

   $rop = get_inf();
   $ropl = get_infl();
    Sets $rop/$ropl to Inf.

   $bool = is_nan($op);
   $bool = is_nanl($op);
    Returns true if $op is a NaN - else returns false

   $bool = is_inf($op);
   $bool = is_infl($op);
    Returns true if $op is -Inf or +Inf - else returns false

=head1 OPERATOR OVERLOADING

   Both Math::Complex_C and Math::Complex_C::Long overload the
   following operators:
    *, +, /, -, **,
    *=, +=, /=, -=, **=,
    not, !, bool,
    ==, !=,
    =, "",
    abs, exp, log, cos, sin, atan2, sqrt

    Cross-class overloading is not impllemented - that is, you
    cannot use both a Math::Complex_C::Long object and a 
    Math::Complex_C object in the same overloaded operation.

    Note: For the purposes of the overloaded 'not', '!' and 'bool'
    operators, a "false" Math::Complex_C object is one with real 
    and imaginary parts that are both "false" - where "false"
    currently means either 0 (including -0) or NaN.
    (A "true" Math::Complex_C object is, of course, simply one
    that is not "false".)

    If the value passed to the overloaded '**', '**=' or
    'equivalence' operators is an IV/UV (integer), it will be
    converted to an NV (floating point value) before the calculation
    is performed.

=head1 LICENSE

   This program is free software; you may redistribute it and/or 
   modify it under the same terms as Perl itself.
   Copyright 2011 Sisyphus.

=head1 AUTHOR

   Sisyphus <sisyphus at(@) cpan dot (.) org>

=cut
