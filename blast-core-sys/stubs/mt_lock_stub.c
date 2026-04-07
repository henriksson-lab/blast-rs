/* Minimal MT_LOCK stubs for standalone BLAST core build.
 * These provide no-op implementations since we handle threading in Rust. */

#include <stdlib.h>

typedef void* MT_LOCK;

/* These symbols are referenced by blast_hspstream.c and blast_diagnostics.c */

MT_LOCK MT_LOCK_Delete(MT_LOCK lock) {
    (void)lock;
    return NULL;
}

int MT_LOCK_DoInternal(MT_LOCK lock, int how) {
    (void)lock;
    (void)how;
    return 1; /* success */
}
