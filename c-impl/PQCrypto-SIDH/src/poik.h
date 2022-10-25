
#if POIK_BITS == 434
#include "P434/P434_poik_api.h"
#elif POIK_BITS == 503
#include "P503/P503_poik_api.h"
#elif POIK_BITS == 610
#include "P610/P610_poik_api.h"
#elif POIK_BITS == 751
#include "P751/P751_poik_api.h"
#else
#error "please #define POIK_BITS to {434,503,610,751} prior to including <poik.h>"
#endif

#ifndef POIK_H
#define POIK_H

////////////////////////////////////////////////////////////////

#if POIK_BITS == 434
    #define POIK_REPS       219
    #define POIK_WIDTH      4
    #define POIK_HEIGHT     7
    #define HASH_BYTES      32      // length of hash commitments
    #define HIDE_BYTES      64      // entropy in hash commitments
#elif POIK_BITS == 503
    #define POIK_REPS       219
    #define POIK_WIDTH      4
    #define POIK_HEIGHT     7
    #define HASH_BYTES      32      // length of hash commitments
    #define HIDE_BYTES      64      // entropy in hash commitments
#elif POIK_BITS == 610
    #define POIK_REPS       329
    #define POIK_WIDTH      4
    #define POIK_HEIGHT     7
    #define HASH_BYTES      48      // length of hash commitments
    #define HIDE_BYTES      96      // entropy in hash commitments
#elif POIK_BITS == 751
    #define POIK_REPS       438
    #define POIK_WIDTH      4
    #define POIK_HEIGHT     7
    #define HASH_BYTES      64      // length of hash commitments
    #define HIDE_BYTES      128     // entropy in hash commitments
#else
#error "please #define POIK_BITS to {434,503,610,751} prior to including <poik.h>"
#endif

////////////////////////////////////////////////////////////////

typedef union {
    unsigned char bytes[SIDH_BYTES];
    struct {
        unsigned char real[SIDH_BYTES/2];
        unsigned char imag[SIDH_BYTES/2];
    };
} packed_fp2;
_Static_assert(sizeof(packed_fp2) == SIDH_BYTES, "");

struct secret {
    packed_fp2 E[POIK_WIDTH+1];
    packed_fp2 P[POIK_WIDTH];
};

typedef struct {
    packed_fp2 x[POIK_WIDTH];
} packed_chain_A;

typedef struct {
    packed_fp2 x[POIK_HEIGHT];
} packed_chain_B;

struct proof {
    packed_fp2 E0, E1;
    struct response {
        enum {
            RESP_LEFT   = -1,
            RESP_BOTTOM =  0,
            RESP_RIGHT  = +1,
        } tag;
        union {
            struct {
                unsigned char r2[HIDE_BYTES], h3[HASH_BYTES];
                packed_chain_B chain;
            } left;
            struct {
                unsigned char r2[HIDE_BYTES], r3[HIDE_BYTES];
                packed_fp2 E2;
                packed_chain_A chain;
            } bottom;
            struct {
                unsigned char h2[HASH_BYTES], r3[HIDE_BYTES];
                packed_chain_B chain;
            } right;
        };
    } resp[POIK_REPS];
};

////////////////////////////////////////////////////////////////

#include <stdbool.h>
#include <stdio.h>

extern FILE *Progress_Indicator_File;

bool Generate(packed_fp2 const *E0, struct secret *sec);
void Prove(struct secret const *sec, struct proof *prf);
bool Verify(struct proof *prf);

bool Dump_Fp2(FILE *file, packed_fp2 const *el);
bool Dump(FILE *file, struct proof const *prf);
bool Parse_Fp2(FILE *file, packed_fp2 *el);
bool Parse(FILE *file, struct proof *prf);

////////////////////////////////////////////////////////////////

#endif  // POIK_H

