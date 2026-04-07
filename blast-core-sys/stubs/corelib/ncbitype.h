/* Minimal ncbitype.h stub — delegates to the real NCBI header */
#ifndef CORELIB___NCBITYPE__H
#define CORELIB___NCBITYPE__H

/* The real ncbitype.h includes <ncbiconf.h> and defines Int1..Int8 etc.
   We just forward to it since our stubs/ is on the include path first
   and provides ncbiconf.h. */
#include <ncbiconf.h>

#include <inttypes.h>

typedef          char  Char;
typedef signed   char  Schar;
typedef unsigned char  Uchar;

typedef int8_t   Int1;
typedef uint8_t  Uint1;
typedef int16_t  Int2;
typedef uint16_t Uint2;
typedef int32_t  Int4;
typedef uint32_t Uint4;
typedef int64_t  Int8;
typedef uint64_t Uint8;

#define NCBI_INT8_IS_LONG 1

typedef long double Ncbi_BigScalar;

#define NCBI_CONST_INT8(v)     INT64_C(v)
#define NCBI_CONST_UINT8(v)    UINT64_C(v)
#define NCBI_INT8_FORMAT_SPEC  PRId64
#define NCBI_UINT8_FORMAT_SPEC PRIu64
#define NCBI_CONST_LONGDOUBLE(v) v##L

#endif /* CORELIB___NCBITYPE__H */
