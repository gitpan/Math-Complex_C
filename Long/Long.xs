#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include <complex.h>

#define MATH_COMPLEX long double _Complex
#define MATH_COMPLEX_DOUBLE long double
#define creal creall
#define cimag cimagl
#define carg cargl
#define cabs cabsl
#define conj conjl
#define cacos cacosl
#define casin casinl
#define catan catanl
#define ccos ccosl
#define csin csinl
#define ctan ctanl
#define cacosh cacoshl
#define casinh casinhl
#define catanh catanhl
#define ccosh ccoshl
#define csinh csinhl
#define ctanh ctanhl
#define cexp cexpl
#define clog clogl
#define cpow cpowl
#define csqrt csqrtl
#define cproj cprojl

#ifdef OLDPERL
#define SvUOK SvIsUV
#endif

int _is_nan(long double x) {
    if(x == x) return 0;
    return 1;
}

int _is_inf(long double x) {
    if(x == 0) return 0;
    if(x / x == x / x) return 0;
    return 1;
}

long double _get_nan() {
    long double _Complex c;
    long double nan, inf;
    __real__ c = 1.0;
    __imag__ c = 0.0;

    inf = creal(c) / cimag(c);
    nan = inf / inf;
    return nan;       
}

long double _get_inf() {
    long double _Complex c;
    long double inf;
    __real__ c = 1.0;
    __imag__ c = 0.0;

    inf = creal(c) / cimag(c);
    return inf;      
}

SV * create_cl() {

     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in create() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C::Long");

     __real__ *pc = _get_nan();
     __imag__ *pc = _get_nan();

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);
     return obj_ref;

}

void assign_cl(SV * rop, SV * d1, SV * d2) {
     long double _d1, _d2;

     _d1 = SvNV(d1);
     _d2 = SvNV(d2);

     __real__ *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = _d1;
     __imag__ *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = _d2;
}

void mul_cl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) *
                                                   *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op2))));
}

void mul_c_nvl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) * SvNV(op2);
}

void mul_c_ivl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) * SvIV(op2);
}

void mul_c_uvl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) * SvUV(op2);
}

void div_cl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) /
                                                   *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op2))));
}

void div_c_nvl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) / SvNV(op2);
}

void div_c_ivl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) / SvIV(op2);
}

void div_c_uvl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) / SvUV(op2);
}

void add_cl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) +
                                                   *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op2))));
}

void add_c_nvl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) + SvNV(op2);
}

void add_c_ivl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) + SvIV(op2);
}

void add_c_uvl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) + SvUV(op2);
}

void sub_cl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) -
                                                   *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op2))));
}

void sub_c_nvl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) - SvNV(op2);
}

void sub_c_ivl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) - SvIV(op2);
}

void sub_c_uvl(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) - SvUV(op2);
}

void DESTROY(SV *  rop) {
     Safefree(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))));
}

SV * real_cl(SV * rop) {
     return newSVnv(creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))));
}

SV * imag_cl(SV * rop) {
     return newSVnv(cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))));
}

SV * arg_cl(SV * rop) {
     return newSVnv(carg(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))));
}

SV * abs_cl(SV * rop) {
     return newSVnv(cabs(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))));
}

void conj_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = conj(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void acos_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = cacos(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void asin_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = casin(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void atan_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = catan(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void cos_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = ccos(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void sin_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = csin(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void tan_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = ctan(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void acosh_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = cacosh(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void asinh_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = casinh(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void atanh_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = catanh(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void cosh_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = ccosh(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void sinh_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = csinh(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void tanh_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = ctanh(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void exp_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))= cexp(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void log_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = clog(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void sqrt_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = csqrt(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void proj_cl(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = cproj(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void pow_cl(SV * rop, SV * op, SV * exp) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = cpow(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))),
                                                        *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(exp)))));
}

SV * _overload_true(SV * rop, SV * second, SV * third) {
     if (_is_nan(creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))))) &&
         _is_nan(cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))))) return newSVuv(0);
     if(creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))) ||
        cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))))) return newSVuv(1);
     return newSVuv(0);
}

SV * _overload_not(SV * rop, SV * second, SV * third) {
     if (_is_nan(creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))))) &&
         _is_nan(cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))))) return newSVuv(1);
     if(creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))) ||
        cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))))) return newSVuv(0);
     return newSVuv(1);
}

SV * _overload_equiv(SV * a, SV * b, SV * third) {
     if(SvUOK(b) || SvIOK(b) || SvNOK(b)) {
       if(SvNV(b) == creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))))) &&
          0       == cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))))) return newSVuv(1);
       return newSVuv(0);
     }
     if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C::Long")) {
         if(creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))))) == creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))))) &&
            cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))))) == cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))))))
              return newSVuv(1);
         return newSVuv(0);
       }
     }
     croak("Invalid argument supplied to Math::Complex_C::Long::_overload_equiv function");
}

SV * _overload_not_equiv(SV * a, SV * b, SV * third) {
     if(SvUOK(b) || SvIOK(b) || SvNOK(b)) {
       if(SvNV(b) == creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))))) &&
          0       == cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))))) return newSVuv(0);
       return newSVuv(1);
     }

     if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C::Long")) {
         if(creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))))) == creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))))) &&
            cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))))) == cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))))))
              return newSVuv(0);
         return newSVuv(1);
       }
     }
     croak("Invalid argument supplied to Math::Complex_C::Long::_overload_not_equiv function");
}


SV * _overload_pow(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc, t;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_pow() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C::Long");

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);

     if(SvNOK(b) || SvIOK(b) || SvUOK(b) ) {
       __real__ t = SvNV(b);
       __imag__ t = 0.0;
       *pc = cpow(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))), t);
       return obj_ref;
     }
     else if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C::Long")) {
         *pc = cpow(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))), *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b)))));
         return obj_ref;
       }
     }
     else croak("Invalid argument supplied to Math::Complex_C::Long::_overload_pow function");
}

SV * _overload_mul(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_mul() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C::Long");

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);

     if(SvUOK(b)) {
       *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) * SvUV(b);
       return obj_ref;
     }

     if(SvIOK(b)) {
       *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) * SvIV(b);
       return obj_ref;
     }

     if(SvNOK(b)) {
       *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) * SvNV(b);
       return obj_ref;
     }

     if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C::Long")) {
         *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) * *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return obj_ref;
       }
     }

     croak("Invalid argument supplied to Math::Complex_C::Long::_overload_mul function");
}

SV * _overload_add(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_add() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C::Long");

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);

     if(SvUOK(b)) {
       *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) + SvUV(b);
       return obj_ref;
     }

     if(SvIOK(b)) {
       *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) + SvIV(b);
       return obj_ref;
     }

     if(SvNOK(b)) {
       *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) + SvNV(b);
       return obj_ref;
     }

     if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C::Long")) {
         *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) + *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return obj_ref;
       }
     }

     croak("Invalid argument supplied to Math::Complex_C::Long::_overload_add function");
}

SV * _overload_div(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_div() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C::Long");

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);

     if(SvUOK(b)) {
       if(SvUV(third)) *pc = SvUV(b) / *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))));
       else *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) / SvUV(b);
       return obj_ref;
     }

     if(SvIOK(b)) {
       if(SvUV(third)) *pc = SvIV(b) / *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))));
       else *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) / SvIV(b);
       return obj_ref;
     }

     if(SvNOK(b)) {
       if(SvUV(third)) *pc = SvNV(b) / *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))));
       else *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) / SvNV(b);
       return obj_ref;
     }

     if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C::Long")) {
         *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) / *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return obj_ref;
       }
     }

     croak("Invalid argument supplied to Math::Complex_C::Long::_overload_div function");
}

SV * _overload_sub(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_sub() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C::Long");

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);

     if(SvUOK(b)) {
       if(SvUV(third)) *pc = SvUV(b) - *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))));
       else *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) - SvUV(b);
       return obj_ref;
     }

     if(SvIOK(b)) {
       if(SvUV(third)) *pc = SvIV(b) - *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))));
       else *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) - SvIV(b);
       return obj_ref;
     }

     if(SvNOK(b)) {
       if(SvUV(third)) *pc = SvNV(b) - *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))));
       else *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) - SvNV(b);
       return obj_ref;
     }

     if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C::Long")) {
         *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) - *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return obj_ref;
       }
     }

     croak("Invalid argument supplied to Math::Complex_C::Long::_overload_sub function");
}

SV * _overload_sqrt(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_sqrt() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C::Long");

     *pc = csqrt(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))));

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);
     return obj_ref;
}

SV * _overload_pow_eq(SV * a, SV * b, SV * third) {
     MATH_COMPLEX t;
     SvREFCNT_inc(a);     

     if(SvNOK(b) || SvIOK(b) || SvUOK(b)) {
       __real__ t = SvNV(b);
       __imag__ t = 0.0;
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) = cpow(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))), t);
       return a;
     }
     else if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C::Long")) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) = cpow(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))),
                                                        *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b)))));
         return a;
       }
     }
     else {
       SvREFCNT_dec(a);
       croak("Invalid argument supplied to Math::Complex_C::Long::_overload_pow_eq function");
     }
}

SV * _overload_mul_eq(SV * a, SV * b, SV * third) {
     SvREFCNT_inc(a);     

     if(SvUOK(b)) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) *= SvUV(b);
       return a;
     }

     if(SvIOK(b)) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) *= SvIV(b);
       return a;
     }

     if(SvNOK(b)) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) *= SvNV(b);
       return a;
     }

     if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C::Long")) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) *= *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return a;
       }
     }

     SvREFCNT_dec(a);
     croak("Invalid argument supplied to Math::Complex_C::Long::_overload_mul_eq function");
}

SV * _overload_add_eq(SV * a, SV * b, SV * third) {
     SvREFCNT_inc(a);     

     if(SvUOK(b)) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) += SvUV(b);
       return a;
     }

     if(SvIOK(b)) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) += SvIV(b);
       return a;
     }

     if(SvNOK(b)) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) += SvNV(b);
       return a;
     }

     if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C::Long")) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) += *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return a;
       }
     }

     SvREFCNT_dec(a);
     croak("Invalid argument supplied to Math::Complex_C::Long::_overload_add_eq function");
}

SV * _overload_div_eq(SV * a, SV * b, SV * third) {
     SvREFCNT_inc(a);     

     if(SvUOK(b)) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) /= SvUV(b);
       return a;
     }

     if(SvIOK(b)) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) /= SvIV(b);
       return a;
     }

     if(SvNOK(b)) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) /= SvNV(b);
       return a;
     }

     if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C::Long")) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) /= *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return a;
       }
     }

     SvREFCNT_dec(a);
     croak("Invalid argument supplied to Math::Complex_C::Long::_overload_div_eq function");
}

SV * _overload_sub_eq(SV * a, SV * b, SV * third) {
     SvREFCNT_inc(a);     

     if(SvUOK(b)) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) -= SvUV(b);
       return a;
     }

     if(SvIOK(b)) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) -= SvIV(b);
       return a;
     }

     if(SvNOK(b)) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) -= SvNV(b);
       return a;
     }

     if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C::Long")) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) -= *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return a;
       }
     }

     SvREFCNT_dec(a);
     croak("Invalid argument supplied to Math::Complex_C::Long::_overload_sub_eq function");
}

SV * _overload_copy(SV * a, SV * second, SV * third) {

     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_copy() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C::Long");

     __real__ *pc = __real__ *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))));
     __imag__ *pc = __imag__ *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))));

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);
     return obj_ref;

}

SV * _overload_abs(SV * rop, SV * second, SV * third) {
     return newSVnv(cabs(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))));
}

SV * _overload_exp(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_exp() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C::Long");

     *pc = cexp(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))));

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);
     return obj_ref;
}

SV * _overload_log(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_log() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C::Long");

     *pc = clog(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))));

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);
     return obj_ref;
}

SV * _overload_sin(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_sin() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C::Long");

     *pc = csin(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))));

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);
     return obj_ref;
}

SV * _overload_cos(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_cos() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C::Long");

     *pc = ccos(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))));

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);
     return obj_ref;
}

SV * _overload_atan2(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_atan2() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C::Long");

     *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) / 
           *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b)))) ;

     *pc = catan(*pc);

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);
     return obj_ref;
}

SV * get_nanl() {
     return newSVnv(_get_nan());
}

SV * get_infl() {
     return newSVnv(_get_inf());
}

SV * is_nanl(SV * a) {
     if(SvNV(a) == SvNV(a)) return newSVuv(0);
     return newSVuv(1);
}

SV * is_infl(SV * a) {
     if(SvNV(a) == 0) return newSVuv(0);
     if(SvNV(a) / SvNV(a) == SvNV(a) / SvNV(a)) return newSVuv(0);
     return newSVuv(1);
}

SV * _complex_type() {
    return newSVpv("long double _Complex", 0);
}

SV * _double_type() {
    return newSVpv("long double", 0);
}

SV * _get_nv(SV * x) {
     return newSVnv(SvNV(x));
}

SV * _which_package(SV * b) {
     if(sv_isobject(b)) return newSVpv(HvNAME(SvSTASH(SvRV(b))), 0);
     return newSVpv("Not an object", 0);
}

MODULE = Math::Complex_C::Long	PACKAGE = Math::Complex_C::Long	

PROTOTYPES: DISABLE


SV *
create_cl ()

void
assign_cl (rop, d1, d2)
	SV *	rop
	SV *	d1
	SV *	d2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	assign_cl(rop, d1, d2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
mul_cl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	mul_cl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
mul_c_nvl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	mul_c_nvl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
mul_c_ivl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	mul_c_ivl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
mul_c_uvl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	mul_c_uvl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
div_cl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	div_cl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
div_c_nvl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	div_c_nvl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
div_c_ivl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	div_c_ivl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
div_c_uvl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	div_c_uvl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
add_cl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	add_cl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
add_c_nvl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	add_c_nvl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
add_c_ivl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	add_c_ivl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
add_c_uvl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	add_c_uvl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sub_cl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sub_cl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sub_c_nvl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sub_c_nvl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sub_c_ivl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sub_c_ivl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sub_c_uvl (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sub_c_uvl(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
DESTROY (rop)
	SV *	rop
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	DESTROY(rop);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

SV *
real_cl (rop)
	SV *	rop

SV *
imag_cl (rop)
	SV *	rop

SV *
arg_cl (rop)
	SV *	rop

SV *
abs_cl (rop)
	SV *	rop

void
conj_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	conj_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
acos_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	acos_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
asin_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	asin_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
atan_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	atan_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
cos_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	cos_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sin_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sin_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
tan_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	tan_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
acosh_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	acosh_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
asinh_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	asinh_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
atanh_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	atanh_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
cosh_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	cosh_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sinh_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sinh_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
tanh_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	tanh_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
exp_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	exp_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
log_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	log_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sqrt_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sqrt_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
proj_cl (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	proj_cl(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
pow_cl (rop, op, exp)
	SV *	rop
	SV *	op
	SV *	exp
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	pow_cl(rop, op, exp);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

SV *
_overload_true (rop, second, third)
	SV *	rop
	SV *	second
	SV *	third

SV *
_overload_not (rop, second, third)
	SV *	rop
	SV *	second
	SV *	third

SV *
_overload_equiv (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_not_equiv (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_pow (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_mul (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_add (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_div (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_sub (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_sqrt (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_pow_eq (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_mul_eq (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_add_eq (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_div_eq (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_sub_eq (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_copy (a, second, third)
	SV *	a
	SV *	second
	SV *	third

SV *
_overload_abs (rop, second, third)
	SV *	rop
	SV *	second
	SV *	third

SV *
_overload_exp (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_log (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_sin (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_cos (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
_overload_atan2 (a, b, third)
	SV *	a
	SV *	b
	SV *	third

SV *
get_nanl ()

SV *
get_infl ()

SV *
is_nanl (a)
	SV *	a

SV *
is_infl (a)
	SV *	a

SV *
_complex_type ()

SV *
_double_type ()

SV *
_get_nv (x)
	SV *	x

SV *
_which_package (b)
	SV *	b

