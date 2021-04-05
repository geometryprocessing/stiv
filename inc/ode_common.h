/*************************************************************************
 *                                                                       *
 * Open Dynamics Engine, Copyright (C) 2001,2002 Russell L. Smith.       *
 * All rights reserved.  Email: russ@q12.org   Web: www.q12.org          *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of EITHER:                                  *
 *   (1) The GNU Lesser General Public License as published by the Free  *
 *       Software Foundation; either version 2.1 of the License, or (at  *
 *       your option) any later version. The text of the GNU Lesser      *
 *       General Public License is included with this library in the     *
 *       file LICENSE.TXT.                                               *
 *   (2) The BSD-style license that is included with this library in     *
 *       the file LICENSE-BSD.TXT.                                       *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the files    *
 * LICENSE.TXT and LICENSE-BSD.TXT for more details.                     *
 *                                                                       *
 *************************************************************************/

#ifndef _ODE_COMMON_H_
#define _ODE_COMMON_H_


/* Define to one of `_getb67', `GETB67', `getb67' for Cray-2 and Cray-YMP
   systems. This function is required for `alloca.c' support on those systems.
   */
/* #undef CRAY_STACKSEG_END */

/* Define to 1 if using `alloca.c'. */
/* #undef C_ALLOCA */

/* Define to 1 if you have `alloca', as a function or macro. */
#define HAVE_ALLOCA 1

/* Define to 1 if you have <alloca.h> and it should be used (not on Ultrix).
   */
#define HAVE_ALLOCA_H 1

/* Use the Apple OpenGL framework. */
#define HAVE_APPLE_OPENGL_FRAMEWORK 1

/* Define to 1 if you have the `atan2f' function. */
#define HAVE_ATAN2F 1

/* Define to 1 if you have the `copysign' function. */
#define HAVE_COPYSIGN 1

/* Define to 1 if you have the `copysignf' function. */
#define HAVE_COPYSIGNF 1

/* Define to 1 if you have the `cosf' function. */
#define HAVE_COSF 1

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
/* #undef HAVE_DOPRNT */

/* Define to 1 if you have the `fabsf' function. */
#define HAVE_FABSF 1

/* Define to 1 if you have the <float.h> header file. */
#define HAVE_FLOAT_H 1

/* Define to 1 if you have the `floor' function. */
#define HAVE_FLOOR 1

/* Define to 1 if you have the `fmodf' function. */
#define HAVE_FMODF 1

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <GL/glext.h> header file. */
/* #undef HAVE_GL_GLEXT_H */

/* Define to 1 if you have the <GL/glu.h> header file. */
/* #undef HAVE_GL_GLU_H */

/* Define to 1 if you have the <GL/gl.h> header file. */
/* #undef HAVE_GL_GL_H */

/* Define to 1 if you have the <ieeefp.h> header file. */
/* #undef HAVE_IEEEFP_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `isnan' function. */
#define HAVE_ISNAN 1

/* Define to 1 if you have the `isnanf' function. */
/* #undef HAVE_ISNANF */

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#define HAVE_MALLOC 1

/* Define to 1 if you have the <malloc.h> header file. */
/* #undef HAVE_MALLOC_H */

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the `memmove' function. */
#define HAVE_MEMMOVE 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#define HAVE_MEMSET 1

/* Define to 1 if libc includes obstacks. */
/* #undef HAVE_OBSTACK */

/* Define to 1 if your system has a GNU libc compatible `realloc' function,
   and to 0 otherwise. */
#define HAVE_REALLOC 1

/* Define to 1 if you have the `select' function. */
#define HAVE_SELECT 1

/* Define to 1 if you have the `sinf' function. */
#define HAVE_SINF 1

/* Define to 1 if you have the `snprintf' function. */
#define HAVE_SNPRINTF 1

/* Define to 1 if you have the `sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 if you have the `sqrtf' function. */
#define HAVE_SQRTF 1

/* Use SSE Optimizations */
/* #undef HAVE_SSE */

/* Define to 1 if you have the <stdarg.h> header file. */
#define HAVE_STDARG_H 1

/* Define to 1 if stdbool.h conforms to C99. */
#define HAVE_STDBOOL_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/select.h> header file. */
#define HAVE_SYS_SELECT_H 1

/* Define to 1 if you have the <sys/socket.h> header file. */
#define HAVE_SYS_SOCKET_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the <values.h> header file. */
/* #undef HAVE_VALUES_H */

/* Define to 1 if you have the `vprintf' function. */
#define HAVE_VPRINTF 1

/* Define to 1 if you have the `vsnprintf' function. */
#define HAVE_VSNPRINTF 1

/* Define to 1 if the system has the type `_Bool'. */
#define HAVE__BOOL 1

/* Define to 1 if you have the `_isnan' function. */
/* #undef HAVE__ISNAN */

/* Define to 1 if you have the `_isnanf' function. */
/* #undef HAVE__ISNANF */

/* Define to 1 if you have the `__isnan' function. */
#define HAVE___ISNAN 1

/* Define to 1 if you have the `__isnanf' function. */
#define HAVE___ISNANF 1

/* Name of package */
#define PACKAGE "ODE"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "ode@ode.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "ODE"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "ODE 0.9.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "ode"

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.9.0"

/* is this a pentium on a gcc-based platform? */
#define PENTIUM 1

/* Define to the type of arg 1 for `select'. */
#define SELECT_TYPE_ARG1 int

/* Define to the type of args 2, 3 and 4 for `select'. */
#define SELECT_TYPE_ARG234 (fd_set *)

/* Define to the type of arg 5 for `select'. */
#define SELECT_TYPE_ARG5 (struct timeval *)

/* The size of `char', as computed by sizeof. */
#define SIZEOF_CHAR 1

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long int', as computed by sizeof. */
#define SIZEOF_LONG_INT 4

/* The size of `short', as computed by sizeof. */
#define SIZEOF_SHORT 2

/* The size of `void*', as computed by sizeof. */
#define SIZEOF_VOIDP 4

/* The extension for shared libraries. */
#define SO_EXT ".dylib"

/* If using the C implementation of alloca, define if you know the
   direction of stack growth for your system; otherwise it will be
   automatically deduced at runtime.
	STACK_DIRECTION > 0 => grows toward higher addresses
	STACK_DIRECTION < 0 => grows toward lower addresses
	STACK_DIRECTION = 0 => direction of growth unknown */
/* #undef STACK_DIRECTION */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "0.9.0"

/* Define to 1 if your processor stores words with the most significant byte
   first (like Motorola and SPARC, unlike Intel and VAX). */
/* #undef WORDS_BIGENDIAN */

/* is this a X86_64 system on a gcc-based platform? */
/* #undef X86_64_SYSTEM */

/* Define to 1 if the X Window System is missing or not being used. */
/* #undef X_DISPLAY_MISSING */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Use double precision */
#define dDOUBLE 

/* dEpsilon Constant */
#define dEpsilon DBL_EPSILON

/* Use gyroscopic terms */
#define dGYROSCOPIC 

/* dInfinity Constant */
#define dInfinity DBL_MAX

/* Disable debug output */
#define dNODEBUG 

/* Use single precision */
/* #undef dSINGLE */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */

/* Define to rpl_realloc if the replacement function should be used. */
/* #undef realloc */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to empty if the keyword `volatile' does not work. Warning: valid
   code using `volatile' can become incorrect without. Disable with care. */
/* #undef volatile */



#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif
#if defined(HAVE_IEEEFP_H) && !defined(__CYGWIN__)
// This header creates conflicts with math.h in Cygwin.
#include <ieeefp.h>
#endif
#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_MATH_H
#include <math.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#ifdef HAVE_STDARG_H
#include <stdarg.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_VALUES_H
#include <values.h>
#endif
#ifdef HAVE_FLOAT_H
#include <float.h>
#endif
#if SIZEOF_CHAR == 1
typedef char int8;
typedef unsigned char uint8;
#else
#error "expecting sizeof(char) == 1"
#endif
#if SIZEOF_SHORT == 2
typedef short int16;
typedef unsigned short uint16;
#else
#error "can not find 2 byte integer type"
#endif
/* integer types (we assume int >= 32 bits) */
#if SIZEOF_INT == 4
typedef short int32;
/*typedef unsigned short uint32;*/
#else
#error "can not find 4 byte integer type"
#endif
/* an integer type that we can safely cast a pointer to and
 * from without loss of bits.
 */
#if SIZEOF_SHORT == SIZEOF_VOIDP
typedef unsigned short intP;
#elif SIZEOF_INT == SIZEOF_VOIDP
typedef unsigned int intP;
#elif SIZEOF_LONG_INT == SIZEOF_VOIDP
typedef unsigned long int intP;
#endif

/* 
Handle Windows DLL odities
Its easier to export all symbols using the -shared flag
for MinGW than differentiating with declspec,
so only do it for MSVC
*/
#if defined(ODE_DLL) && defined(WIN32) && defined(_MSC_VER)
#define ODE_API __declspec( dllexport )
#elif !defined(ODE_DLL) && defined(WIN32) && defined(MSC_VER)
#define ODE_API __declspec( dllimport )
#else
#define ODE_API
#endif



//#include "config.h"
//#include "error.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif


/* configuration stuff */

/* the efficient alignment. most platforms align data structures to some
 * number of bytes, but this is not always the most efficient alignment.
 * for example, many x86 compilers align to 4 bytes, but on a pentium it
 * is important to align doubles to 8 byte boundaries (for speed), and
 * the 4 floats in a SIMD register to 16 byte boundaries. many other
 * platforms have similar behavior. setting a larger alignment can waste
 * a (very) small amount of memory. NOTE: this number must be a power of
 * two. this is set to 16 by default.
 */
#define EFFICIENT_ALIGNMENT 16


/* constants */

/* pi and 1/sqrt(2) are defined here if necessary because they don't get
 * defined in <math.h> on some platforms (like MS-Windows)
 */

#ifndef M_PI
#define M_PI REAL(3.1415926535897932384626433832795029)
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 REAL(0.7071067811865475244008443621048490)
#endif


/* debugging:
 *   IASSERT  is an internal assertion, i.e. a consistency check. if it fails
 *            we want to know where.
 *   UASSERT  is a user assertion, i.e. if it fails a nice error message
 *            should be printed for the user.
 *   AASSERT  is an arguments assertion, i.e. if it fails "bad argument(s)"
 *            is printed.
 *   DEBUGMSG just prints out a message
 */

#ifndef dNODEBUG
#ifdef __GNUC__
#define dIASSERT(a) if (!(a)) dDebug (d_ERR_IASSERT, \
  "assertion \"" #a "\" failed in %s() [%s]",__FUNCTION__,__FILE__);
#define dUASSERT(a,msg) if (!(a)) dDebug (d_ERR_UASSERT, \
  msg " in %s()", __FUNCTION__);
#define dDEBUGMSG(msg) dMessage (d_ERR_UASSERT,				\
msg " in %s() File %s Line %d", __FUNCTION__, __FILE__,__LINE__);
#else
#define dIASSERT(a) if (!(a)) dDebug (d_ERR_IASSERT, \
  "assertion \"" #a "\" failed in %s:%d",__FILE__,__LINE__);
#define dUASSERT(a,msg) if (!(a)) dDebug (d_ERR_UASSERT, \
  msg " (%s:%d)", __FILE__,__LINE__);
#define dDEBUGMSG(msg) dMessage (d_ERR_UASSERT, \
  msg " (%s:%d)", __FILE__,__LINE__);
#endif
#else
#define dIASSERT(a) ;
#define dUASSERT(a,msg) ;
#define dDEBUGMSG(msg) ;
#endif
#define dAASSERT(a) dUASSERT(a,"Bad argument(s)")

// Macro used to suppress unused variable warning
#define dVARIABLEUSED(a) ((void)a)

/* floating point data type, vector, matrix and quaternion types */

#if defined(dSINGLE)
typedef float dReal;
#ifdef dDOUBLE
#error You can only #define dSINGLE or dDOUBLE, not both.
#endif // dDOUBLE
#elif defined(dDOUBLE)
typedef double dReal;
#else
#error You must #define dSINGLE or dDOUBLE
#endif

// Detect if we've got both trimesh engines enabled.
#if dTRIMESH_ENABLED
#if dTRIMESH_OPCODE && dTRIMESH_GIMPACT
#error You can only #define dTRIMESH_OPCODE or dTRIMESH_GIMPACT, not both.
#endif
#endif // dTRIMESH_ENABLED

/* round an integer up to a multiple of 4, except that 0 and 1 are unmodified
 * (used to compute matrix leading dimensions)
 */
#define dPAD(a) (((a) > 1) ? ((((a)-1)|3)+1) : (a))

/* these types are mainly just used in headers */
typedef dReal dVector3[4];
typedef dReal dVector4[4];
typedef dReal dMatrix3[4*3];
typedef dReal dMatrix4[4*4];
typedef dReal dMatrix6[8*6];
typedef dReal dQuaternion[4];


/* precision dependent scalar math functions */

#if defined(dSINGLE)

#define REAL(x) (x ## f)					/* form a constant */
#define dRecip(x) ((1.0f/(x)))				/* reciprocal */
#define dSqrt(x) (sqrtf(x))			/* square root */
#define dRecipSqrt(x) ((1.0f/sqrtf(x)))		/* reciprocal square root */
#define dSin(x) (sinf(x))				/* sine */
#define dCos(x) (cosf(x))				/* cosine */
#define dFabs(x) (fabsf(x))			/* absolute value */
#define dAtan2(y,x) (atan2f(y,x))		/* arc tangent with 2 args */
#define dFMod(a,b) (fmodf(a,b))		/* modulo */
#define dFloor(x) floorf(x)			/* floor */

#ifdef HAVE___ISNANF
#define dIsNan(x) (__isnanf(x))
#elif defined(HAVE__ISNANF)
#define dIsNan(x) (_isnanf(x))
#elif defined(HAVE_ISNANF)
#define dIsNan(x) (isnanf(x))
#else
  /*
     fall back to _isnan which is the VC way,
     this may seem redundant since we already checked
     for _isnan before, but if isnan is detected by
     configure but is not found during compilation
     we should always make sure we check for __isnanf,
     _isnanf and isnanf in that order before falling
     back to a default
  */
#define dIsNan(x) (_isnan(x))
#endif

#define dCopySign(a,b) ((dReal)copysignf(a,b))

#elif defined(dDOUBLE)

#define REAL(x) (x)
#define dRecip(x) (1.0/(x))
#define dSqrt(x) sqrt(x)
#define dRecipSqrt(x) (1.0/sqrt(x))
#define dSin(x) sin(x)
#define dCos(x) cos(x)
#define dFabs(x) fabs(x)
#define dAtan2(y,x) atan2((y),(x))
#define dFMod(a,b) (fmod((a),(b)))
#define dFloor(x) floor(x)

#ifdef HAVE___ISNAN
#define dIsNan(x) (__isnan(x))
#elif defined(HAVE__ISNAN)
#define dIsNan(x) (_isnan(x))
#elif defined(HAVE_ISNAN)
#define dIsNan(x) (isnan(x))
#else
#define dIsNan(x) (_isnan(x))
#endif

#define dCopySign(a,b) (copysign((a),(b)))

#else
#error You must #define dSINGLE or dDOUBLE
#endif


/* utility */


/* round something up to be a multiple of the EFFICIENT_ALIGNMENT */

#define dEFFICIENT_SIZE(x) ((((x)-1)|(EFFICIENT_ALIGNMENT-1))+1)


/* alloca aligned to the EFFICIENT_ALIGNMENT. note that this can waste
 * up to 15 bytes per allocation, depending on what alloca() returns.
 */

#define dALLOCA16(n) \
  ((char*)dEFFICIENT_SIZE(((size_t)(alloca((n)+(EFFICIENT_ALIGNMENT-1))))))


// Use the error-checking memory allocation system.  Because this system uses heap
//  (malloc) instead of stack (alloca), it is slower.  However, it allows you to
//  simulate larger scenes, as well as handle out-of-memory errors in a somewhat
//  graceful manner

// #define dUSE_MALLOC_FOR_ALLOCA

#ifdef dUSE_MALLOC_FOR_ALLOCA
enum {
  d_MEMORY_OK = 0,		/* no memory errors */
  d_MEMORY_OUT_OF_MEMORY	/* malloc failed due to out of memory error */
};

#endif



#ifdef __cplusplus
}
#endif

#endif
