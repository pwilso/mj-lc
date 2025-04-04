// use C++ mode in Emacs: -*- C++ -*-

#ifndef PMIFloat_HEADER
#define PMIFloat_HEADER

static const char PMIFloat_SccsId[] =
  "$Id: //tcad/gemini/D-2010.03/src/dessis/src/PMIFloat.h#1 $, "
  "Copyright Synopsys, Inc. 2010-2010";

// definitions for the physical model interface in Sentaurus Device

#include <iostream>
#include <limits>

template <class des_t_float> class pmi_float_value;
class dd_real;
class qd_real;

#ifdef isnan
#undef isnan
#endif

#ifdef isinf
#undef isinf
#endif



/// precision of floating point arithmetic
enum pmi_e_precision {
  pmi_c_np,  ///< normal precision (double): -ExtendedPrecision
  pmi_c_xp,  ///< extended precision (long double): ExtendedPrecision
  pmi_c_dd,  ///< double-double: ExtendedPrecision(128)
  pmi_c_qd   ///< quad-double: ExtendedPrecision(256)
};



/// floating point class supporting extended precision and automatic differentiation
class pmi_float {

private:

  static pmi_e_precision precision;

  union {
    pmi_float_value<double>* value_np;
    pmi_float_value<long double>* value_xp;
    pmi_float_value<dd_real>* value_dd;
    pmi_float_value<qd_real>* value_qd;
  };

public:
  // constructors
  pmi_float ();
  pmi_float (const int a);
  pmi_float (const long int a);
  pmi_float (const size_t a);
  pmi_float (const double a);
  pmi_float (const long double a);
  pmi_float (const dd_real a);
  pmi_float (const qd_real a);

  /// constructor with gradient of size 'size_gradient', gradient is 0
  template <class des_t_float> pmi_float (const des_t_float a, size_t size_gradient);

  /// constructor with gradient of size 'size_gradient', component 'index' of gradient is 1
  template <class des_t_float> pmi_float (const des_t_float a, size_t size_gradient, size_t index);

  /// copy constructor
  pmi_float (const pmi_float& a);

  /// destructor
  ~pmi_float ();

  /// precision of floating point arithmetic
  static pmi_e_precision get_precision () { return precision; }

  // only for internal use in Sentaurus Device
  static void set_precision (pmi_e_precision prec) { precision = prec; }

  // constants
  static pmi_float min ();  ///< smallest positive number
  static pmi_float max ();  ///< largest positive number
  static pmi_float epsilon ();  ///< precision
  static pmi_float inf ();  ///< infinity
  static pmi_float nan ();  ///< not a number (quiet NaN)
  static pmi_float signaling_nan ();  ///< not a number (signaling NaN)

  // access to values
  template <class des_t_float> void set_value (const des_t_float a);
  template <class des_t_float> void set_value (const des_t_float a, size_t size_gradient);
  template <class des_t_float> void set_value (const des_t_float a, size_t size_gradient, size_t index);
  template <class des_t_float> void set_gradient (size_t i, const des_t_float a);

  template <class des_t_float> des_t_float get_value () const;
  size_t size_gradient () const;
  template <class des_t_float> des_t_float get_gradient (size_t i) const;

  // assignment
  pmi_float& operator = (const pmi_float& a);

  // unary operators
  const pmi_float& operator + () const;  ///< prefix unary plus
  pmi_float operator - () const;  ///< prefix unary minus

  pmi_float& operator += (const pmi_float& a);
  pmi_float& operator -= (const pmi_float& a);
  pmi_float& operator *= (const pmi_float& a);
  pmi_float& operator /= (const pmi_float& a);

  pmi_float& add (const pmi_float& a);
  pmi_float& sub (const pmi_float& a);
  pmi_float& mul (const pmi_float& a);
  pmi_float& div (const pmi_float& a);

  pmi_float& abs ();
  pmi_float& acos ();
  pmi_float& acosh ();
  pmi_float& asin ();
  pmi_float& asinh ();
  pmi_float& atan ();
  pmi_float& atanh ();
  pmi_float& cos ();
  pmi_float& cosh ();
  pmi_float& erf ();
  pmi_float& erfc ();
  pmi_float& exp ();
  pmi_float& expm1 ();
  pmi_float& ldexp (const int exp);
  pmi_float& log ();
  pmi_float& log1p ();
  pmi_float& log10 ();
  pmi_float& pow (const pmi_float& a);
  pmi_float& pow_int (const int n);
  pmi_float& sin ();
  pmi_float& sinh ();
  pmi_float& sqrt ();
  pmi_float& tan ();
  pmi_float& tanh ();

  // binary operators
  friend pmi_float operator + (const pmi_float& a, const pmi_float& b);
  friend pmi_float operator - (const pmi_float& a, const pmi_float& b);
  friend pmi_float operator * (const pmi_float& a, const pmi_float& b);
  friend pmi_float operator / (const pmi_float& a, const pmi_float& b);

  // comparisons
  friend bool operator == (const pmi_float& a, const pmi_float& b);
  friend bool operator != (const pmi_float& a, const pmi_float& b);
  friend bool operator < (const pmi_float& a, const pmi_float& b);
  friend bool operator <= (const pmi_float& a, const pmi_float& b);
  friend bool operator > (const pmi_float& a, const pmi_float& b);
  friend bool operator >= (const pmi_float& a, const pmi_float& b);

  // functions
  friend pmi_float abs (const pmi_float& a);
  friend pmi_float acos (const pmi_float& a);
  friend pmi_float acosh (const pmi_float& a);
  friend pmi_float asin (const pmi_float& a);
  friend pmi_float asinh (const pmi_float& a);
  friend pmi_float atan (const pmi_float& a);
  friend pmi_float atanh (const pmi_float& a);
  friend pmi_float atan2 (const pmi_float& y, const pmi_float& x);
  friend pmi_float cos (const pmi_float& a);
  friend pmi_float cosh (const pmi_float& a);
  friend pmi_float erf (const pmi_float& a);
  friend pmi_float erfc (const pmi_float& a);
  friend pmi_float exp (const pmi_float& a);
  friend pmi_float expm1 (const pmi_float& a);
  friend pmi_float hypot (const pmi_float& x, const pmi_float& y);
  friend bool isinf (const pmi_float& a);
  friend bool isnan (const pmi_float& a);
  friend pmi_float ldexp (const pmi_float& a, const int exp);
  friend pmi_float log (const pmi_float& a);
  friend pmi_float log1p (const pmi_float& a);
  friend pmi_float log10 (const pmi_float& a);
  friend pmi_float pow (const pmi_float& x, const pmi_float& y);
  friend pmi_float pow_int (const pmi_float& x, const int n);
  friend pmi_float sin (const pmi_float& a);
  friend pmi_float sinh (const pmi_float& a);
  friend pmi_float sqrt (const pmi_float& a);
  friend pmi_float tan (const pmi_float& a);
  friend pmi_float tanh (const pmi_float& a);

  // output
  friend std::ostream& operator << (std::ostream& s, const pmi_float& f);

};



namespace std {

  template <> class numeric_limits<pmi_float> {

  public:
    static const bool is_specialized = true;

    static const bool is_signed = true;
    static const bool is_integer = false;
    static const bool is_exact = false;

    static pmi_float min () { return pmi_float::min (); }
    static pmi_float max () { return pmi_float::max (); }
    static pmi_float epsilon () { return pmi_float::epsilon (); }

    static pmi_float infinity () { return pmi_float::inf (); }
    static pmi_float quiet_NaN () { return pmi_float::nan (); }
    static pmi_float signaling_NaN () { return pmi_float::signaling_nan (); }

  };

}

#endif
