/* Minimal ncbiconf.h stub for standalone BLAST core build */
#ifndef NCBICONF__H
#define NCBICONF__H

#define SIZEOF_CHAR        1
#define SIZEOF_SHORT       2
#define SIZEOF_INT         4
#define SIZEOF_LONG        8
#define SIZEOF_LONG_LONG   8
#define SIZEOF_DOUBLE      8
#define SIZEOF_LONG_DOUBLE 16
#define SIZEOF_VOIDP       8
#define SIZEOF_SIZE_T      8

#define HAVE_INTTYPES_H    1
#define HAVE_INTPTR_T      1
#define HAVE_UINTPTR_T     1
#define STDC_HEADERS       1

#define NCBI_CXX_TOOLKIT   1
#define NCBI_OS_UNIX       1

/* Disable visibility attributes for static linking */
/* #define HAVE_ATTRIBUTE_VISIBILITY_DEFAULT 1 */

#endif /* NCBICONF__H */
