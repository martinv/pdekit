#ifndef PDEKIT_Common_PDEKit_hpp
#define PDEKIT_Common_PDEKit_hpp

#include <complex>

/// Definition of PDEKit namespace.
namespace pdekit
{

// ----------------------------------------------------------------------------

/// Definition of the basic types for possible portability conflicts

/// float
typedef float Float;
/// double
typedef double Double;
/// long double
typedef long double LDouble;
/// int
typedef int Int;
/// unsigned int
typedef unsigned int Uint;
/// long int
typedef long int Lint;
/// long long int
typedef long long int LLint;
/// long unsigned int
typedef long unsigned int LUint;
/// long long unsigned int
typedef long long unsigned int LLUint;
/// short unsigned int
typedef short unsigned int SUint;
/// short signed int
typedef short int SInt;
/// char
typedef char Char;
/// unsigned char
typedef unsigned char Uchar;

/// Enumeration of the dimensions
enum SCDim
{
  DIM_0D,
  DIM_1D,
  DIM_2D,
  DIM_3D
};
/// Enumeration of the coordinates indexes
enum CoordXYZ
{
  X,
  Y,
  Z
};
/// Enumeration of the coordinate indexes
enum CoordX0X1X2
{
  X0,
  X1,
  X2
};
/// Enumeration of the reference coordinates indexes
enum CoordRefXiEtaZeta
{
  KSI,
  ETA,
  ZTA
};
/// Enumeration of the reference coordinates indexes
enum CoordRefXi0Xi1Xi2
{
  XI0,
  XI1,
  XI2
};

// ----------------------------------------------------------------------------

/// Definition of the default precision
#ifdef PDEKIT_PRECISION_LONG_DOUBLE
typedef LDouble Real;
#else
#ifdef PDEKIT_PRECISION_DOUBLE
typedef Double Real;
#else
#ifdef SC_PRECISION_SINGLE
typedef Float Real;
#endif
#endif
#endif
// if nothing defined, use double
#if !defined PDEKIT_PRECISION_DOUBLE && !defined PDEKIT_PRECISION_SINGLE &&                        \
    !defined PDEKIT_PRECISION_LONG_DOUBLE
typedef Double Real;
#endif

typedef std::complex<Real> Complex;

// ----------------------------------------------------------------------------

const Real PI  = 3.14159265358979323846;
const Real EXP = 2.71828182845904523536;

#if 0
/* Some useful constants.  */
#if defined __USE_BSD || defined __USE_XOPEN
#define M_E 2.7182818284590452354         /* e */
#define M_LOG2E 1.4426950408889634074     /* log_2 e */
#define M_LOG10E 0.43429448190325182765   /* log_10 e */
#define M_LN2 0.69314718055994530942      /* log_e 2 */
#define M_LN10 2.30258509299404568402     /* log_e 10 */
#define M_PI 3.14159265358979323846       /* pi */
#define M_PI_2 1.57079632679489661923     /* pi/2 */
#define M_PI_4 0.78539816339744830962     /* pi/4 */
#define M_1_PI 0.31830988618379067154     /* 1/pi */
#define M_2_PI 0.63661977236758134308     /* 2/pi */
#define M_2_SQRTPI 1.12837916709551257390 /* 2/sqrt(pi) */
#define M_SQRT2 1.41421356237309504880    /* sqrt(2) */
#define M_SQRT1_2 0.70710678118654752440  /* 1/sqrt(2) */
#endif

/* The above constants are not adequate for computation using `long double's.
   Therefore we provide as an extension constants with similar names as a
   GNU extension.  Provide enough digits for the 128-bit IEEE quad.  */
#ifdef __USE_GNU
#define M_El 2.718281828459045235360287471352662498L        /* e */
#define M_LOG2El 1.442695040888963407359924681001892137L    /* log_2 e */
#define M_LOG10El 0.434294481903251827651128918916605082L   /* log_10 e */
#define M_LN2l 0.693147180559945309417232121458176568L      /* log_e 2 */
#define M_LN10l 2.302585092994045684017991454684364208L     /* log_e 10 */
#define M_PIl 3.141592653589793238462643383279502884L       /* pi */
#define M_PI_2l 1.570796326794896619231321691639751442L     /* pi/2 */
#define M_PI_4l 0.785398163397448309615660845819875721L     /* pi/4 */
#define M_1_PIl 0.318309886183790671537767526745028724L     /* 1/pi */
#define M_2_PIl 0.636619772367581343075535053490057448L     /* 2/pi */
#define M_2_SQRTPIl 1.128379167095512573896158903121545172L /* 2/sqrt(pi) */
#define M_SQRT2l 1.414213562373095048801688724209698079L    /* sqrt(2) */
#define M_SQRT1_2l 0.707106781186547524400844362104849039L  /* 1/sqrt(2) */
#endif

#endif

// ----------------------------------------------------------------------------

} // namespace pdekit

// ----------------------------------------------------------------------------

#endif // PDEKIT_Common_PDEKit_hpp
