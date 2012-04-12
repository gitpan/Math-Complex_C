#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include <complex.h>

#define MATH_COMPLEX double _Complex
#define MATH_COMPLEX_DOUBLE double

#ifdef OLDPERL
#define SvUOK SvIsUV
#endif

#ifndef Newx
#  define Newx(v,n,t) New(0,v,n,t)
#endif

#ifndef Newxz
#  define Newxz(v,n,t) Newz(0,v,n,t)
#endif

int _is_nan(double x) {
    if(x == x) return 0;
    return 1;
}

int _is_inf(double x) {
    if(x == 0) return 0;
    if(x / x == x / x) return 0;
    return 1;
}

double _get_nan(void) {
    double _Complex c;
    double nan, inf;
    __real__ c = 1.0;
    __imag__ c = 0.0;

    inf = creal(c) / cimag(c);
    nan = inf / inf;
    return nan;       
}

double _get_inf(void) {
    double _Complex c;
    double inf;
    __real__ c = 1.0;
    __imag__ c = 0.0;

    inf = creal(c) / cimag(c);
    return inf;      
}

SV * create_c(void) {

     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in create() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C");

     __real__ *pc = _get_nan();
     __imag__ *pc = _get_nan();

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);
     return obj_ref;

}

void assign_c(SV * rop, SV * d1, SV * d2) {
     double _d1, _d2;

     _d1 = SvNV(d1);
     _d2 = SvNV(d2);

     __real__ *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = _d1;
     __imag__ *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = _d2;
}

void mul_c(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) *
                                                   *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op2))));
}

void mul_c_nv(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) * SvNV(op2);
}

void mul_c_iv(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) * SvIV(op2);
}

void mul_c_uv(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) * SvUV(op2);
}

void div_c(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) /
                                                   *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op2))));
}

void div_c_nv(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) / SvNV(op2);
}

void div_c_iv(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) / SvIV(op2);
}

void div_c_uv(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) / SvUV(op2);
}

void add_c(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) +
                                                   *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op2))));
}

void add_c_nv(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) + SvNV(op2);
}

void add_c_iv(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) + SvIV(op2);
}

void add_c_uv(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) + SvUV(op2);
}

void sub_c(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) -
                                                   *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op2))));
}

void sub_c_nv(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) - SvNV(op2);
}

void sub_c_iv(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) - SvIV(op2);
}

void sub_c_uv(SV * rop, SV * op1, SV * op2) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op1)))) - SvUV(op2);
}

void DESTROY(SV *  rop) {
     Safefree(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))));
}

SV * real_c(SV * rop) {
     return newSVnv(creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))));
}

SV * imag_c(SV * rop) {
     return newSVnv(cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))));
}

SV * arg_c(SV * rop) {
     return newSVnv(carg(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))));
}

SV * abs_c(SV * rop) {
     return newSVnv(cabs(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))));
}

void conj_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = conj(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void acos_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = cacos(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void asin_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = casin(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void atan_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = catan(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void cos_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = ccos(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void sin_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = csin(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void tan_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = ctan(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void acosh_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = cacosh(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void asinh_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = casinh(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void atanh_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = catanh(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void cosh_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = ccosh(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void sinh_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = csinh(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void tanh_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = ctanh(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void exp_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop))))= cexp(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void log_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = clog(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void sqrt_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = csqrt(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void proj_c(SV * rop, SV * op) {
     *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(rop)))) = cproj(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(op)))));
}

void pow_c(SV * rop, SV * op, SV * exp) {
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
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C")) {
         if(creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))))) == creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))))) &&
            cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))))) == cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))))))
              return newSVuv(1);
         return newSVuv(0);
       }
     }
     croak("Invalid argument supplied to Math::Complex_C::_overload_equiv function");
}

SV * _overload_not_equiv(SV * a, SV * b, SV * third) {
     if(SvUOK(b) || SvIOK(b) || SvNOK(b)) {
       if(SvNV(b) == creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))))) &&
          0       == cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))))) return newSVuv(0);
       return newSVuv(1);
     }
     if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C")) {
         if(creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))))) == creal(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))))) &&
            cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a))))) == cimag(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))))))
              return newSVuv(0);
         return newSVuv(1);
       }
     }
     croak("Invalid argument supplied to Math::Complex_C::_overload_not_equiv function");
}


SV * _overload_pow(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc, t;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_pow() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C");

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);

     if(SvNOK(b) || SvIOK(b) || SvUOK(b) ) {
       __real__ t = SvNV(b);
       __imag__ t = 0.0;
       *pc = cpow(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))), t);
       return obj_ref;
     }
     else if(sv_isobject(b)) {
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C")) {
         *pc = cpow(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))), *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b)))));
         return obj_ref;
       }
     }
     else croak("Invalid argument supplied to Math::Complex_C::_overload_pow function");
}

SV * _overload_mul(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_mul() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C");

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
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C")) {
         *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) * *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return obj_ref;
       }
     }

     croak("Invalid argument supplied to Math::Complex_C::_overload_mul function");
}

SV * _overload_add(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_add() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C");

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
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C")) {
         *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) + *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return obj_ref;
       }
     }

     croak("Invalid argument supplied to Math::Complex_C::_overload_add function");
}

SV * _overload_div(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_div() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C");

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
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C")) {
         *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) / *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return obj_ref;
       }
     }

     croak("Invalid argument supplied to Math::Complex_C::_overload_div function");
}

SV * _overload_sub(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_sub() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C");

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
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C")) {
         *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) - *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return obj_ref;
       }
     }

     croak("Invalid argument supplied to Math::Complex_C::_overload_sub function");
}

SV * _overload_sqrt(SV * a, SV * b, SV * third) {
     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_sqrt() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C");

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
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C")) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) = cpow(*(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))),
                                                        *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b)))));
         return a;
       }
     }
     else {
       SvREFCNT_dec(a);
       croak("Invalid argument supplied to Math::Complex_C::_overload_pow_eq function");
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
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C")) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) *= *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return a;
       }
     }

     SvREFCNT_dec(a);
     croak("Invalid argument supplied to Math::Complex_C::_overload_mul_eq function");
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
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C")) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) += *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return a;
       }
     }

     SvREFCNT_dec(a);
     croak("Invalid argument supplied to Math::Complex_C::_overload_add_eq function");
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
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C")) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) /= *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return a;
       }
     }

     SvREFCNT_dec(a);
     croak("Invalid argument supplied to Math::Complex_C::_overload_div_eq function");
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
       if(strEQ(HvNAME(SvSTASH(SvRV(b))), "Math::Complex_C")) {
       *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) -= *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b))));
         return a;
       }
     }

     SvREFCNT_dec(a);
     croak("Invalid argument supplied to Math::Complex_C::_overload_sub_eq function");
}

SV * _overload_copy(SV * a, SV * second, SV * third) {

     MATH_COMPLEX *pc;
     SV * obj_ref, * obj;

     New(42, pc, 1, MATH_COMPLEX);
     if(pc == NULL) croak("Failed to allocate memory in _overload_copy() function");

     obj_ref = newSV(0);
     obj = newSVrv(obj_ref, "Math::Complex_C");

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
     obj = newSVrv(obj_ref, "Math::Complex_C");

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
     obj = newSVrv(obj_ref, "Math::Complex_C");

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
     obj = newSVrv(obj_ref, "Math::Complex_C");

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
     obj = newSVrv(obj_ref, "Math::Complex_C");

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
     obj = newSVrv(obj_ref, "Math::Complex_C");

     *pc = *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(a)))) / 
           *(INT2PTR(MATH_COMPLEX *, SvIV(SvRV(b)))) ;

     *pc = catan(*pc);

     sv_setiv(obj, INT2PTR(IV,pc));
     SvREADONLY_on(obj);
     return obj_ref;
}

SV * get_nan(void) {
     return newSVnv(_get_nan());
}

SV * get_inf(void) {
     return newSVnv(_get_inf());
}

SV * is_nan(SV * a) {
     if(SvNV(a) == SvNV(a)) return newSVuv(0);
     return newSVuv(1);
}

SV * is_inf(SV * a) {
     if(SvNV(a) == 0) return newSVuv(0);
     if(SvNV(a) / SvNV(a) == SvNV(a) / SvNV(a)) return newSVuv(0);
     return newSVuv(1);
}

SV * _complex_type(void) {
    return newSVpv("double _Complex", 0);
}

SV * _double_type(void) {
    return newSVpv("double", 0);
}

SV * _get_nv(SV * x) {
     return newSVnv(SvNV(x));
}

SV * _which_package(SV * b) {
     if(sv_isobject(b)) return newSVpv(HvNAME(SvSTASH(SvRV(b))), 0);
     return newSVpv("Not an object", 0);
}

SV * _ivsize(void) {
     return newSViv(sizeof(IV));
}

SV * _nvsize(void) {
     return newSViv(sizeof(NV));
}

SV * _doublesize(void) {
     return newSViv(sizeof(double));
}

SV * _longdoublesize(void) {
     return newSViv(sizeof(long double));
}

SV * _double_Complexsize(void) {
     return newSViv(sizeof(double _Complex));
}

SV * _longdouble_Complexsize(void) {
     return newSViv(sizeof(long double _Complex));
}

SV * is_neg_zero(SV * x) {
     char * buffer;

     if(!SvNOK(x)) {
      warn("Argument passed to is_neg_zero function is not an NV - therefore can't be '-0'");
      return newSVuv(0);
     }

     if(SvNV(x) != 0) return newSVuv(0);

     Newz(42, buffer, 2, char);
     if(buffer == NULL) croak("Failed to allocate memory in is_neg_zero");

     sprintf(buffer, "%.0f", (double)SvNV(x));

     if(strEQ(buffer, "-0")) {
       Safefree(buffer);
       return newSVuv(1);
     }   

     Safefree(buffer);
     return newSVuv(0);
}

/* Attempt to return -0 */
SV * _get_neg_zero(void) {
     double d = 0.0;
     d *= -1.1;

     return newSVnv(d);
}    

MODULE = Math::Complex_C	PACKAGE = Math::Complex_C	

PROTOTYPES: DISABLE


int
_is_nan (x)
	double	x

int
_is_inf (x)
	double	x

double
_get_nan ()
		

double
_get_inf ()
		

SV *
create_c ()
		

void
assign_c (rop, d1, d2)
	SV *	rop
	SV *	d1
	SV *	d2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	assign_c(rop, d1, d2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
mul_c (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	mul_c(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
mul_c_nv (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	mul_c_nv(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
mul_c_iv (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	mul_c_iv(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
mul_c_uv (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	mul_c_uv(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
div_c (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	div_c(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
div_c_nv (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	div_c_nv(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
div_c_iv (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	div_c_iv(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
div_c_uv (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	div_c_uv(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
add_c (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	add_c(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
add_c_nv (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	add_c_nv(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
add_c_iv (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	add_c_iv(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
add_c_uv (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	add_c_uv(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sub_c (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sub_c(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sub_c_nv (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sub_c_nv(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sub_c_iv (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sub_c_iv(rop, op1, op2);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sub_c_uv (rop, op1, op2)
	SV *	rop
	SV *	op1
	SV *	op2
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sub_c_uv(rop, op1, op2);
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
real_c (rop)
	SV *	rop

SV *
imag_c (rop)
	SV *	rop

SV *
arg_c (rop)
	SV *	rop

SV *
abs_c (rop)
	SV *	rop

void
conj_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	conj_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
acos_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	acos_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
asin_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	asin_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
atan_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	atan_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
cos_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	cos_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sin_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sin_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
tan_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	tan_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
acosh_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	acosh_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
asinh_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	asinh_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
atanh_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	atanh_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
cosh_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	cosh_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sinh_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sinh_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
tanh_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	tanh_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
exp_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	exp_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
log_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	log_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
sqrt_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	sqrt_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
proj_c (rop, op)
	SV *	rop
	SV *	op
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	proj_c(rop, op);
	if (PL_markstack_ptr != temp) {
          /* truly void, because dXSARGS not invoked */
	  PL_markstack_ptr = temp;
	  XSRETURN_EMPTY; /* return empty stack */
        }
        /* must have used dXSARGS; list context implied */
	return; /* assume stack size is correct */

void
pow_c (rop, op, exp)
	SV *	rop
	SV *	op
	SV *	exp
	PREINIT:
	I32* temp;
	PPCODE:
	temp = PL_markstack_ptr++;
	pow_c(rop, op, exp);
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
get_nan ()
		

SV *
get_inf ()
		

SV *
is_nan (a)
	SV *	a

SV *
is_inf (a)
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

SV *
_ivsize ()
		

SV *
_nvsize ()
		

SV *
_doublesize ()
		

SV *
_longdoublesize ()
		

SV *
_double_Complexsize ()
		

SV *
_longdouble_Complexsize ()
		

SV *
is_neg_zero (x)
	SV *	x

SV *
_get_neg_zero ()
		
