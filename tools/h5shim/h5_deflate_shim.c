/*
 * LD_PRELOAD shim: override the gzip level pkdgrav3 passes to HDF5 (3.5.1 hardcodes 9 in io/fio.c).
 *   PKD_H5_DEFLATE=<0..9>  level to use; 0 removes the gzip filter entirely; unset = unchanged.
 * Lets us reproduce the production setting (level 1) with an unmodified pkdgrav3 binary.
 * Build (Linux):  gcc -O2 -shared -fPIC -o libh5_deflate_shim.so h5_deflate_shim.c -ldl
 * Use:            mpirun -x LD_PRELOAD=$PWD/libh5_deflate_shim.so -x PKD_H5_DEFLATE=1 -np 1 pkdgrav3 ...
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef int herr_t;
typedef int64_t hid_t;

herr_t H5Pset_deflate(hid_t plist, unsigned level) {
    static herr_t (*real)(hid_t, unsigned) = 0;
    static int announced = 0;
    if (!real) real = (herr_t (*)(hid_t, unsigned))dlsym(RTLD_NEXT, "H5Pset_deflate");
    const char* e = getenv("PKD_H5_DEFLATE");
    unsigned lvl = e ? (unsigned)atoi(e) : level;
    if (!announced) {
        fprintf(stderr, "h5_deflate_shim: H5Pset_deflate level %u -> %u%s\n", level, lvl, lvl ? "" : " (filter removed)");
        announced = 1;
    }
    if (lvl == 0) return 0;
    return real(plist, lvl);
}
