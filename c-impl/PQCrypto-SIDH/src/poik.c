/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: proof of isogeny knowledge (PoIK)
*********************************************************************************************/

#include "poik.h"

typedef struct {
    f2elm_t A, C;
} curve_proj;

typedef struct {
    f2elm_t B;      // codomain Montgomery coefficient
    f2elm_t r, u2;  // isomorphism maps x to (x-r)/u2
} mont_iso;

typedef struct {
    curve_proj E[POIK_WIDTH+1];
    point_proj P[POIK_WIDTH];
} chain_A;

typedef struct {
    curve_proj E[POIK_HEIGHT+1];
    point_proj P[POIK_HEIGHT];
} chain_B;

struct instance
{
    unsigned char r2[HIDE_BYTES], r3[HIDE_BYTES];
    curve_proj E2, E3;
    chain_B left;
    chain_A bottom;
    chain_B right;
};

// domain-separation strings
#define RO_CHALLENGE    "\x01"  // Fiat-Shamir
#define RO_COMMIT       "\x02"  // hash commitments

struct commitment
{
    packed_fp2 j0, j1;
    struct commitment1 {
        unsigned char h2[HASH_BYTES], h3[HASH_BYTES];
    } hs[POIK_REPS];
} __attribute__((__packed__));
_Static_assert(sizeof(struct commitment) == 2*SIDH_BYTES + 2*POIK_REPS*HASH_BYTES, "");

////////////////////////////////////////////////////////////////

//TODO These utility functions probably shouldn't be here.

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define fail(S) do { fprintf(stderr, "%s\n", (S)); abort(); } while (false)

//XXX for debugging
#define dbg(s, obj) do { \
        fprintf(stderr, "(%u) %s: ", __LINE__, s); \
        packed_fp2 _BUF; \
        fp2_encode((obj), _BUF.bytes); \
        Dump_Fp2(stderr, &_BUF); \
        fputc('\n', stderr); \
    } while (false)
//XXX

static void fp2one(f2elm_t el)  //FIXME should be a constant
{
    fpcopy((digit_t *) Montgomery_one, el[0]);
    fpzero(el[1]);
}

static const f2elm_t zero = {0};

static uint8_t fp_equal(const felm_t a, const felm_t b)
{
    felm_t c;
    fpsub(a, b, c);
    fpcorrection(c);
    return ~ct_compare((void const *) c, (void const *) *zero, sizeof(felm_t));
}
static uint8_t fp2_equal(const f2elm_t a, const f2elm_t b)
{
    return fp_equal(a[0], b[0]) & fp_equal(a[1], b[1]);
}
#define fp2_iszero(A) fp2_equal((A), zero)
#define fp_iszero(A) fp_equal((A), *zero)

static uint8_t proj_equal(const f2elm_t X1, const f2elm_t Z1, const f2elm_t X2, const f2elm_t Z2)
{
    f2elm_t U,V;
    fp2mul_mont(X1, Z2, U);
    fp2mul_mont(X2, Z1, V);
    return fp2_equal(U, V);
}
#define curve_equal(E1,E2) proj_equal((E1)->A, (E1)->C, (E2)->A, (E2)->C)
#define point_equal(P, Q) proj_equal((P)->X, (P)->Z, (Q)->X, (Q)->Z)

static void randombytes_strict(unsigned char *buf, size_t len)
{
    if (randombytes(buf, len))
        fail("randombytes() failed -- could not get entropy");
}

static void fprand(felm_t el)
{
    randombytes_strict((void *) el, sizeof(felm_t));
    _Static_assert(NBITS_FIELD >= RADIX * (NWORDS_FIELD - 1), "");
    el[NWORDS_FIELD-1] &= (digit_t) -1 >> NBITS_FIELD % RADIX;  //FIXME not uniform
}

static void fp2div(const f2elm_t X, const f2elm_t Z, f2elm_t x)
{
    assert(!fp2_iszero(Z));
    f2elm_t inv;
    fp2copy(Z, inv);
    fp2inv_mont(inv);
    fp2mul_mont(X, inv, x);
}
static void fp2div_enc(const f2elm_t X, const f2elm_t Z, packed_fp2 *r)
{
    f2elm_t x;
    fp2div(X, Z, x);
    fp2_encode(x, r->bytes);
}

// NB: is_felm_lt is not constant-time, hence this
static uint8_t fp_lt(const felm_t x, const felm_t y)
{
    digit_t z[NWORDS_FIELD];
    unsigned b = mp_sub(x, y, z, NWORDS_FIELD);
    return - (uint8_t) (b & 1);
}

// adapted from sqrt_Fp2 in fpx.c -- <XXX CURRENTLY NOT constant-time XXX> & works for subfield elements
static bool fp2sqrt(const f2elm_t x, f2elm_t y, bool known_square)
{
    digit_t const *a = x[0], *b = x[1];

    felm_t t0, t1;

    if (fp_iszero(b)) {
        fpcopy(a, t1);
        for (unsigned i = 0; i < OALICE_BITS - 2; i++) {
            fpsqr_mont(t1, t1);
        }
        for (unsigned i = 0; i < OBOB_EXPON; i++) {
            fpsqr_mont(t1, t0);
            fpmul_mont(t1, t0, t1);
        }
        // t1 = a^((p+1)/4)
        fpsqr_mont(t1, t0);
        if (fp_equal(t0, a)) {
            fpcopy(t1, y[0]);
            fpzero(y[1]);
        }
        else {
            fpneg(t0);
            assert(fp_equal(t0, a));
            fpcopy(t1, y[1]);
            fpzero(y[0]);
        }
    }
    else {
        felm_t v, w;
        {
            fpsqr_mont(a, t0);
            fpsqr_mont(b, t1);
            fpadd(t0, t1, t0);
            fpcopy(t0, t1);
            for (unsigned i = 0; i < OALICE_BITS - 2; i++) {
                fpsqr_mont(t1, t1);
            }
            for (unsigned i = 0; i < OBOB_EXPON; i++) {
                fpsqr_mont(t1, t0);
                fpmul_mont(t1, t0, t1);
            }
            // t1 = (a^2+b^2)^((p+1)/4)
            fpadd(a, t1, v);
            fpdiv2(v, v);
            fpcorrection(v);
            // v = (a + sqrt(a^2+b^2)) / 2
        }

        felm_t t2;
        fpcopy(v, t2);
        fpinv_chain_mont(t2);               // t2 = v^((p-3)/4)
        fpmul_mont(v, t2, t1);              // t1 = t2*v
        fpmul_mont(t2, b, t2);              // t2 = t2*b
        fpdiv2(t2, t2);                     // t2 = t2/2
        fpsqr_mont(t1, w);                  // w = t1^2
        fpcorrection(w);

        if (memcmp(v, w, NBITS_TO_NBYTES(NBITS_FIELD)) == 0) {
            fpcopy(t1, y[0]);
            fpcopy(t2, y[1]);
        } else {
            fpneg(t1);
            fpcopy(t2, y[0]);
            fpcopy(t1, y[1]);
        }
    }

#ifndef NDEBUG
    {
        f2elm_t chk;
        fp2sqr_mont(y, chk);
        assert(!known_square || fp2_equal(chk, x));
    }
#endif

    if (known_square)
        return true;

    //FIXME this is a hack; should probably use is_sqr_fp2
    f2elm_t chk;
    fp2sqr_mont(y, chk);
    return fp2_equal(chk, x);
}

////////////////////////////////////////////////////////////////

//TODO is this the best way of doing this? (luckily, 2^8 is close to 3^5)
static size_t sample_ternary(uint8_t const *bytes, size_t num_bytes, int8_t *trits, size_t num_trits)
{
    for (size_t done = 0, idx = 0; done < num_trits; ) {
        if (idx >= num_bytes)
            return done;
        uint8_t v = bytes[idx++];
        if (v >= 243)   // 3^5
            continue;
        for (unsigned j = 0; j < 5 && done < num_trits; ++j) {
            int8_t r = v%3;
            trits[done++] = r > 3/2 ? r-3 : r;
            v /= 3;
        }
    }
    return num_trits;
}

////////////////////////////////////////////////////////////////

static void AC_to_AAC24(curve_proj const *E, f2elm_t A24plus, f2elm_t A24minus, f2elm_t C24)
{
    f2elm_t c;
    felm_t *const CC = C24 ? C24 : c;
    fp2add(E->C, E->C, CC);             // 2C
    fp2add(E->A, CC, A24plus);          // A+2C
    fp2add(CC, CC, CC);                 // 4C
    if (A24minus)
        fp2sub(A24plus, CC, A24minus);  // A-2C
}

static void AC24_to_AC(const f2elm_t A24plus, const f2elm_t C24, curve_proj *E)
{
    fp2add(A24plus, A24plus, E->A);     // 2A+4C
    fp2sub(E->A, C24, E->A);            // 2A
    fp2add(E->A, E->A, E->A);           // 4A
    fp2copy(C24, E->C);                 // 4C
}

static void AA24_to_AC(const f2elm_t A24minus, const f2elm_t A24plus, curve_proj *E)
{
    fp2add(A24plus, A24minus, E->A);    // 2A
    fp2add(A24minus, A24minus, E->C);   // 2A-4C
    fp2sub(E->A, E->C, E->C);           // 4C
    fp2add(E->A, E->A, E->A);           // 4A
}

////////////////////////////////////////////////////////////////

//TODO use Elligator
static void random_point(const f2elm_t A, point_proj *P)
{
    f2elm_t x;
    while (true) {
        fprand(x[0]);
        fprand(x[1]);

        f2elm_t rhs;
        fp2add(A, x, rhs);                      // x + A
        fp2mul_mont(rhs, x, rhs);               // x^2 + A*x
        fpadd(rhs[0], Montgomery_one, rhs[0]);  // x^2 + A*x + 1
        fp2mul_mont(rhs, x, rhs);               // x^3 + A*x^2 + x

        f2elm_t y;
        if (fp2sqrt(rhs, y, false))
            break;
    }

    fp2copy(x, P->X);
    fp2one(P->Z);
}

// output: P a random l^e-torsion point, Q = l^(e-1)P
static void random_point_with_order(const f2elm_t A, unsigned l, point_proj *P, unsigned *e, point_proj *Q)
{
    assert(l == 2 || l == 3);

    f2elm_t A24plus, A24minus, C24;
    {
        fp2one(C24);
        fp2add(C24, C24, C24);          // 2
        fp2add(A, C24, A24plus);        // A + 2
        fp2sub(A, C24, A24minus);       // A - 2
        fp2add(C24, C24, C24);          // 4
    }

    for (unsigned count = 0; ++count; ) {
        if (count >= 9999)
            fail("random_point_with_order() appears to be stuck in an infinite loop");

        random_point(A, P);

        switch (l) {
        case 2: xTPLe(P, P, A24minus, A24plus, OBOB_EXPON); break;
        case 3: xDBLe(P, P, A24plus, C24, OALICE_BITS); break;
        default: __builtin_unreachable();
        }
        if (!fp2_iszero(P->Z))
            break;
    }

    *Q = *P;
    *e = 1;

    for (unsigned count = 0; ++count; ++*e) {
        if (count >= 9999)
            fail("random_point_with_order() appears to be stuck in an infinite loop");

        point_proj T;
        switch (l) {
        case 2: xDBL(Q, &T, A24plus, C24); break;
        case 3: xTPL(Q, &T, A24minus, A24plus); break;
        default: __builtin_unreachable();
        }
        if (fp2_iszero(T.Z))
            break;
        *Q = T;
    }

#ifndef NDEBUG
    {
        point_proj T;
        if (l == 2) xDBLe(P, &T, A24plus, C24, *e-1);
        else        xTPLe(P, &T, A24minus, A24plus, *e-1);
        assert(point_equal(Q, &T));
        if (l == 2) xDBL(&T, &T, A24plus, C24);
        else        xTPL(&T, &T, A24minus, A24plus);
        assert(fp2_iszero(T.Z));
    }
#endif
}

////////////////////////////////////////////////////////////////

// x -> (x-r)/u2
static void eval_iso(mont_iso const *iso, point_proj const *P, point_proj *Q)
{
    f2elm_t rz;
    fp2mul_mont(iso->r, P->Z, rz);      // r*z
    fp2sub(P->X, rz, Q->X);             // x-r*z
    fp2mul_mont(iso->u2, P->Z, Q->Z);   // u2*z
}

// only solving x^2 + A*x + 1 = 0 is implemented
static size_t solve_quadratic(const f2elm_t A, f2elm_t *sols)
{
    f2elm_t sqrt_disc;
    {
        f2elm_t four;               //FIXME should be a constant
        fp2one(four);               // 1
        fp2add(four, four, four);   // 2
        fp2add(four, four, four);   // 4

        f2elm_t disc;
        fp2sqr_mont(A, disc);       // A^2
        fp2sub(disc, four, disc);   // A^2-4
        fp2sqrt(disc, sqrt_disc, true); // sqrt(A^2-4)
    }

    fp2sub(sqrt_disc, A, sols[0]);
    fp2div2(sols[0], sols[0]);
    fp2sub(sols[0], sqrt_disc, sols[1]);

    size_t num = 2 - fp2_equal(sols[0], sols[1]);

#ifndef NDEBUG
    {
        f2elm_t chk, one;
        fp2one(one);
        for (size_t i = 0; i < num; ++i) {
            fp2add(A, sols[i], chk);
            fp2mul_mont(chk, sols[i], chk);
            fp2add(chk, one, chk);
            assert(fp2_iszero(chk));
        }
    }
#endif

    return num;
}

static void all_models(curve_proj const *E, mont_iso iso[6])
{
    f2elm_t A;
    fp2div(E->A, E->C, A);

    // A
    fp2zero(iso->r);
    fp2one(iso->u2);
    fp2copy(A, iso->B);
    ++iso;

    // -A
    *iso = iso[-1];
    fp2neg(iso->u2);
    fp2neg(iso->B);
    ++iso;

    f2elm_t x2[2];
    size_t n2 = solve_quadratic(A, x2);
    assert(n2 == 2);

    for (size_t i = 0; i < n2; ++i) {
        f2elm_t min2x2;
        fp2add(x2[i], x2[i], min2x2);   // 2*x2
        fp2neg(min2x2);                 // -2*x2

        f2elm_t x4[2];
        size_t n4 = solve_quadratic(min2x2, x4);
        assert(n4 == 2);

        for (size_t j = 0; j < n4; ++j) {
            fp2sub(x4[j], x2[i], iso->u2);      // x4-x2

            fp2sub(x2[i], min2x2, iso->B);      // 3*x2
            fp2add(iso->B, A, iso->B);          // A+3*x2
            fp2div(iso->B, iso->u2, iso->B);    // (A+3*x2) / (x4-x2)

            fp2copy(x2[i], iso->r);
            ++iso;
        }
    }

#ifndef NDEBUG
    {
        iso -= 6;
        f2elm_t one, j0, j1;
        fp2one(one);
        j_inv(iso[0].B, one, j0);
        for (size_t i = 0; i < 6; ++i) {
            j_inv(iso[i].B, one, j1);
            assert(fp2_equal(j0, j1));
//            for (size_t ii = 0; ii < i; ++ii)
//                assert(!fp2_equal(iso[ii].B, iso[i].B));

            point_proj P, Q;
            random_point(A, &P);

            f2elm_t A24plus, C24;
            AC_to_AAC24(E, A24plus, NULL, C24);
            xDBLe(&P, &Q, A24plus, C24, 55);

            point_proj imP, imQ, imQ_;
            eval_iso(&iso[i], &P, &imP);
            eval_iso(&iso[i], &Q, &imQ);

            curve_proj EE;
            fp2copy(iso[i].B, EE.A);
            fp2one(EE.C);
            AC_to_AAC24(&EE, A24plus, NULL, C24);
            xDBLe(&imP, &imQ_, A24plus, C24, 55);

            assert(point_equal(&imQ, &imQ_));
        }
    }
#endif
}

//TODO could batch inversions
static void canonical_model(curve_proj const *E, mont_iso *iso)
{
    mont_iso all[6];
    all_models(E, all);
    *iso = *all;
    for (size_t i = 1; i < 6; ++i) {
        felm_t re0, im0, re, im;
        from_mont(iso->B[0], re0);
        from_mont(iso->B[1], im0);
        from_mont(all[i].B[0], re);
        from_mont(all[i].B[1], im);
        uint8_t const lt_re = fp_lt(re, re0);
        uint8_t const gt_re = fp_lt(re0, re);
        uint8_t const lt_im = fp_lt(im, im0);
        uint8_t const eq_re = ~lt_re & ~gt_re;
        uint8_t const lt = lt_re | (eq_re & lt_im);
        ct_cmov((void *) iso, (void const *) &all[i], sizeof(*iso), lt);
    }
}

static void change_model(curve_proj const *E, mont_iso *iso, mont_iso *identity)
{
    mont_iso all[6];
    all_models(E, all);
    // 0 is x, 1 is -x, 2..5 are nontrivial
    *iso = all[2];
    if (identity)
        *identity = *all;
}

////////////////////////////////////////////////////////////////

// adapted from Ladder in ec_isogeny.c
void Ladder(point_proj const *P, digit_t const *m, const f2elm_t A, size_t order_bits, point_proj *R)
{
    f2elm_t A24;
    fp2one(A24);
    fpadd(A24[0], A24[0], A24[0]);
    fp2add(A, A24, A24);
    fp2div2(A24, A24);
    fp2div2(A24, A24);  // A24 = (A+2)/4

    // R0 <- O, R1 <- P
    point_proj R0, R1 = *P;
    fp2one(R0.X);
    fp2zero(R0.Z);

    unsigned prevbit = 0;

    // Main loop
    for (ssize_t i = order_bits-1;  i >= 0; --i) {
        unsigned bit = (m[i / RADIX] >> (i % RADIX)) & 1;
        digit_t mask = - (digit_t) (bit ^ prevbit);
        prevbit = bit;

        swap_points(&R0, &R1, mask);
        xDBLADD(&R0, &R1, P->X, P->Z, A24);  // R1 += R0; R0 += R0
    }
    digit_t mask = - (digit_t) prevbit;
    swap_points(&R0, &R1, mask);

    *R = R0;
}

////////////////////////////////////////////////////////////////

// computes EE = E/<K> and evaluates phi_K at Q
// K must have order 2^e2
// back will be set to the order-2 point on the codomain that leads to backtracking
static bool Isogeny_A(curve_proj const *E, point_proj const *K, curve_proj *EE, point_proj *Q, point_proj *back)
{
    assert(!fp2_iszero(E->C));

    f2elm_t A24plus, C24;
    AC_to_AAC24(E, A24plus, NULL, C24);

    point_proj R = *K;

    // The 2-isogeny formulas fail when the input point lies above
    // the point (0,0). We switch to a different Montgomery model
    // in that case.
    {
        point_proj T;
        xDBLe(&R, &T, A24plus, C24, OALICE_BITS-1);
        uint8_t above_00 = fp2_iszero(T.X);

        mont_iso iso;
        {
            mont_iso other;
            change_model(E, &other, &iso);
            ct_cmov((void *) &iso, (void const *) &other, sizeof(iso), above_00);
        }
        curve_proj EE;
        fp2copy(iso.B, EE.A);
        fp2one(EE.C);
        AC_to_AAC24(&EE, A24plus, NULL, C24);
        eval_iso(&iso, &R, &R);
        if (Q) eval_iso(&iso, Q, Q);
    }

    // Check that R has the correct order.
    {
        point_proj T, O;
        xDBLe(&R, &T, A24plus, C24, OALICE_BITS-1);
        xDBL(&T, &O, A24plus, C24);
        uint8_t bad = fp2_iszero(T.Z) | ~fp2_iszero(O.Z);
        if (bad)
            return false;
    }

#if (OALICE_BITS % 2 == 1)
    {
        point_proj S;

        xDBLe(&R, &S, A24plus, C24, OALICE_BITS-1);
        get_2_isog(&S, A24plus, C24);
        eval_2_isog(&R, &S);
        if (Q) eval_2_isog(Q, &S);
    }
#endif

    point_proj pts[MAX_INT_POINTS_ALICE];
    unsigned pts_index[MAX_INT_POINTS_ALICE];

    unsigned index = 0, npts = 0, ii = 0;

    f2elm_t coeff[3];

    for (unsigned row = 1; row < MAX_Alice; row++) {
        while (index < MAX_Alice-row) {
            fp2copy(R.X, pts[npts].X);
            fp2copy(R.Z, pts[npts].Z);
            pts_index[npts++] = index;
            unsigned m = strat_Alice[ii++];
            xDBLe(&R, &R, A24plus, C24, 2*m);
            index += m;
        }

        get_4_isog(&R, A24plus, C24, coeff);

        for (unsigned i = 0; i < npts; i++)
            eval_4_isog(&pts[i], coeff);
        if (Q) eval_4_isog(Q, coeff);

        fp2copy(pts[npts-1].X, R.X);
        fp2copy(pts[npts-1].Z, R.Z);
        index = pts_index[npts-1];
        npts -= 1;
    }

    get_4_isog(&R, A24plus, C24, coeff);
    if (Q) eval_4_isog(Q, coeff);

    if (EE) AC24_to_AC(A24plus, C24, EE);

    if (back) {
        fp2zero(back->X);
        fp2one(back->Z);
    }

    return true;
}

// computes EE = E/<K> and evaluates phi_K at Q
// K must have order 3^e3
// back will be set to the order-3 point on the codomain that leads to backtracking
static bool Isogeny_B(curve_proj const *E, point_proj const *K, curve_proj *EE, point_proj *Q, point_proj *back)
{
    assert(!fp2_iszero(E->C));

    f2elm_t A24plus, A24minus;
    AC_to_AAC24(E, A24plus, A24minus, NULL);

    point_proj R = *K;

    // Check that R has the correct order.
    {
        point_proj T, O;
        xTPLe(&R, &T, A24minus, A24plus, OBOB_EXPON-1);
        xTPL(&T, &O, A24minus, A24plus);
        uint8_t bad = fp2_iszero(T.Z) | ~fp2_iszero(O.Z);
        if (bad)
            return false;
    }

    point_proj pts[MAX_INT_POINTS_BOB];
    unsigned pts_index[MAX_INT_POINTS_BOB];

    unsigned index = 0, npts = 0, ii = 0;

    f2elm_t coeff[3];

    for (unsigned row = 1; row < MAX_Bob; row++) {
        while (index < MAX_Bob-row) {
            fp2copy(R.X, pts[npts].X);
            fp2copy(R.Z, pts[npts].Z);
            pts_index[npts++] = index;
            unsigned m = strat_Bob[ii++];
            xTPLe(&R, &R, A24minus, A24plus, m);
            index += m;
        }

        get_3_isog(&R, A24minus, A24plus, coeff);

        for (unsigned i = 0; i < npts; i++) {
            eval_3_isog(&pts[i], coeff);
        }
        if (Q) eval_3_isog(Q, coeff);

        fp2copy(pts[npts-1].X, R.X);
        fp2copy(pts[npts-1].Z, R.Z);
        index = pts_index[npts-1];
        npts -= 1;
    }

    if (back) {
        f2elm_t A;
        {
            curve_proj E;
            AA24_to_AC(A24minus, A24plus, &E);
            fp2div(E.A, E.C, A);
        }
        do {
            point_proj P;
            unsigned e;
            random_point_with_order(A, 3, &P, &e, back);
            assert(e >= 1);
        } while (point_equal(back, &R));
    }

    get_3_isog(&R, A24minus, A24plus, coeff);
    if (Q) eval_3_isog(Q, coeff);
    if (back) eval_3_isog(back, coeff);
    assert(!back || !fp2_iszero(back->Z));

    if (EE) AA24_to_AC(A24minus, A24plus, EE);

    return true;
}

////////////////////////////////////////////////////////////////

// random point of order l^el
static void Random_Point_AB(curve_proj const *E, unsigned l, point_proj *P, point_proj const *avoid)
{
    assert(l == 2 || l == 3);

    f2elm_t A;
    fp2div(E->A, E->C, A);  //TODO can reuse in caller?

    unsigned const emax = l == 2 ? OALICE_BITS : l == 3 ? OBOB_EXPON : -1;

    for (unsigned count = 0; count < 9999; ++count) {

        unsigned e;
        point_proj Q;
        random_point_with_order(A, l, P, &e, &Q);

        if (e < emax)
            continue;

        if (avoid && point_equal(&Q, avoid))
            continue;

        return;
    }

    fail("Random_Point_AB() appears to be stuck in an infinite loop");
}

static void Random_Point_A(curve_proj const* E, point_proj *P, point_proj const *avoid) { Random_Point_AB(E, 2, P, avoid); }
static void Random_Point_B(curve_proj const* E, point_proj *P, point_proj const *avoid) { Random_Point_AB(E, 3, P, avoid); }

////////////////////////////////////////////////////////////////

#ifndef NDEBUG
#include <math.h>
#endif
static bool proves_order(unsigned long order_2part, unsigned long order_3part)
{
    const unsigned long LOG2_3_N = 1054, LOG2_3_D = 665;
    assert((LOG2_3_N + 0.) / LOG2_3_D <= log2(3));
    assert((LOG2_3_N + 1.) / LOG2_3_D >= log2(3));

    // need a divisor > 4p of (p+1)^2
    const unsigned long bound = LOG2_3_D * (2 + OALICE_BITS) + (LOG2_3_N+1) * OBOB_EXPON,
                        value = LOG2_3_D * order_2part       + (LOG2_3_N+0) * order_3part;
    return value > bound;
}

static bool Check_Curve(curve_proj const *E)
{
#ifndef NDEBUG
    f2elm_t one;
    fp2one(one);
    assert(fp2_equal(E->C, one));
#endif

    f2elm_t A24plus, A24minus, C24;
    AC_to_AAC24(E, A24plus, A24minus, C24);

    struct torsion_data {
        point_proj T;
        unsigned h;
    } max2 = {0}, max3 = {0};

    for (unsigned count = 0; count < 9999; ++count) {
        point_proj P, Q;
        random_point(E->A, &P);

        struct torsion_data cur2, cur3;

        // determine 2-valuation of order
        xTPLe(&P, &Q, A24minus, A24plus, OBOB_EXPON);
        for (cur2.h = 0; cur2.h < OALICE_BITS; ++cur2.h) {
            if (fp2_iszero(Q.Z))
                break;
            cur2.T = Q;
            xDBL(&Q, &Q, A24plus, C24);
        }
        if (!fp2_iszero(Q.Z))
            return false;

        // (only) independent components add up
        unsigned ord2 = max2.h;
        if (!point_equal(&cur2.T, &max2.T))
            ord2 += cur2.h;
        if (proves_order(ord2, max3.h))
            return true;

        // determine 3-valuation of order
        xDBLe(&P, &Q, A24plus, C24, cur2.h);    // NB: uses known 2-valuation
        for (cur3.h = 0; cur3.h < OBOB_EXPON; ++cur3.h) {
            if (fp2_iszero(Q.Z))
                break;
            cur3.T = Q;
            xTPL(&Q, &Q, A24minus, A24plus);
        }
        if (!fp2_iszero(Q.Z))
            return false;

        // (only) independent components add up
        unsigned ord3 = max3.h;
        if (!point_equal(&cur3.T, &max3.T))
            ord3 += cur3.h;
        if (proves_order(ord2, ord3))
            return true;

        if (cur2.h > max2.h)
            max2 = cur2;
        if (cur3.h > max3.h)
            max3 = cur3;
    }

    fail("Check_Curve() appears to be stuck in an infinite loop");
}

////////////////////////////////////////////////////////////////

//FIXME some code for 2-/3-torsion is needlessly duplicated

static void Random_Isogeny_Chain_A(curve_proj const *E0, chain_A *chain)
{
    point_proj const Pback = {{0}, {1}};

    chain->E[0] = *E0;
    for (size_t i = 0; i < POIK_WIDTH; ++i) {
        Random_Point_A(&chain->E[i], &chain->P[i], i ? &Pback : NULL);
#ifndef NDEBUG
        point_proj Pback2;
        Isogeny_A(&chain->E[i], &chain->P[i], &chain->E[i+1], NULL, &Pback2);
        assert(point_equal(&Pback, &Pback2));
#else
        Isogeny_A(&chain->E[i], &chain->P[i], &chain->E[i+1], NULL, NULL);
#endif
    }

#ifndef NDEBUG  // ensure there is no backtracking
    point_proj P, Q;
    Random_Point_A(&chain->E[0], &P, NULL);
    {
        f2elm_t A24plus, C24;
        AC_to_AAC24(&chain->E[0], A24plus, NULL, C24);
        point_proj U, V, T;
        xDBLe(&P, &U, A24plus, C24, OALICE_BITS-1);
        xDBL(&U, &T, A24plus, C24); assert(fp2_iszero(T.Z));
        Random_Point_A(&chain->E[0], &Q, &U);
        xDBLe(&Q, &V, A24plus, C24, OALICE_BITS-1); assert(!point_equal(&U, &V));
        xDBL(&V, &T, A24plus, C24); assert(fp2_iszero(T.Z));
    }
    for (size_t i = 0; i < POIK_WIDTH; ++i) {
        {
            curve_proj E;
            assert(Isogeny_A(&chain->E[i], &chain->P[i], &E, &P, NULL)); assert(curve_equal(&E, &chain->E[i+1]));
            assert(Isogeny_A(&chain->E[i], &chain->P[i], &E, &Q, NULL)); assert(curve_equal(&E, &chain->E[i+1]));
        }
    }
    f2elm_t A24plus, C24;
    AC_to_AAC24(&chain->E[POIK_WIDTH], A24plus, NULL, C24);
    size_t ordP, ordQ;
    for (ordP = 0; !fp2_iszero(P.Z); ++ordP) xDBL(&P, &P, A24plus, C24); assert(fp2_iszero(P.Z));
    for (ordQ = 0; !fp2_iszero(Q.Z); ++ordQ) xDBL(&Q, &Q, A24plus, C24); assert(fp2_iszero(Q.Z));
    assert(ordP == OALICE_BITS || ordQ == OALICE_BITS);  // kernel does not contain 2-torsion
#endif
}

static void Random_Isogeny_Chain_B(curve_proj const *E0, chain_B *chain)
{
    point_proj Pback;

    chain->E[0] = *E0;
    for (size_t i = 0; i < POIK_HEIGHT; ++i) {
        Random_Point_B(&chain->E[i], &chain->P[i], i ? &Pback : NULL);
        Isogeny_B(&chain->E[i], &chain->P[i], &chain->E[i+1], NULL, &Pback);
    }

#ifndef NDEBUG  // ensure there is no backtracking
    point_proj P, Q;
    Random_Point_B(&chain->E[0], &P, NULL);
    {
        f2elm_t A24plus, A24minus;
        AC_to_AAC24(&chain->E[0], A24plus, A24minus, NULL);
        point_proj U, V, T;
        xTPLe(&P, &U, A24minus, A24plus, OBOB_EXPON-1);
        xTPL(&U, &T, A24minus, A24plus); assert(fp2_iszero(T.Z));
        Random_Point_B(&chain->E[0], &Q, &U);
        xTPLe(&Q, &V, A24minus, A24plus, OBOB_EXPON-1); assert(!point_equal(&U, &V));
        xTPL(&V, &T, A24minus, A24plus); assert(fp2_iszero(T.Z));
    }
    for (size_t i = 0; i < POIK_HEIGHT; ++i) {
        {
            curve_proj E;
            assert(Isogeny_B(&chain->E[i], &chain->P[i], &E, &P, NULL)); assert(curve_equal(&E, &chain->E[i+1]));
            assert(Isogeny_B(&chain->E[i], &chain->P[i], &E, &Q, NULL)); assert(curve_equal(&E, &chain->E[i+1]));
        }
    }
    f2elm_t A24plus, A24minus;
    AC_to_AAC24(&chain->E[POIK_HEIGHT], A24plus, A24minus, NULL);
    size_t ordP, ordQ;
    for (ordP = 0; !fp2_iszero(P.Z); ++ordP) xTPL(&P, &P, A24minus, A24plus); assert(fp2_iszero(P.Z));
    for (ordQ = 0; !fp2_iszero(Q.Z); ++ordQ) xTPL(&Q, &Q, A24minus, A24plus); assert(fp2_iszero(Q.Z));
    assert(ordP == OBOB_EXPON || ordQ == OBOB_EXPON);  // kernel does not contain 3-torsion
#endif
}

////////////////////////////////////////////////////////////////

static inline void pack_curve_and_point(curve_proj const *E, point_proj const *P, packed_fp2 *A, packed_fp2 *x)
{
    assert(E);
    mont_iso iso;
    canonical_model(E, &iso);
    if (A)
        fp2_encode(iso.B, A->bytes);
    if (!x)
        return;
    assert(P);
    point_proj Q;
    eval_iso(&iso, P, &Q);
    fp2div_enc(Q.X, Q.Z, x);
}
static void pack_curve(curve_proj const *E, packed_fp2 *a) { pack_curve_and_point(E, NULL, a, NULL); }
static void pack_point(curve_proj const *E, point_proj const *P, packed_fp2 *x) { pack_curve_and_point(E, P, NULL, x); }
static bool unpack_curve(packed_fp2 const *a, curve_proj *E, bool check) { fp2_decode(a->bytes, E->A); fp2one(E->C); return !check || Check_Curve(E); }
static void unpack_point(packed_fp2 const *x, point_proj *P) { fp2_decode(x->bytes, P->X); fp2one(P->Z); }

static void Pack_Chain_A(chain_A const *in, packed_chain_A *out)
{
    for (size_t i = 0; i < POIK_WIDTH; ++i) {
        digit_t scalar[NWORDS_ORDER];
        {
            unsigned char buf[SIDH_SECRETKEYBYTES_A];
            random_mod_order_A(buf);
            buf[0] |= 1;    // make it odd
            decode_to_digits(buf, scalar, SECRETKEY_A_BYTES, NWORDS_ORDER);
        }

        f2elm_t A;
        fp2div((&in->E[i])->A, (&in->E[i])->C, A);

        point_proj Q;
        Ladder(&in->P[i], scalar, A, OALICE_BITS, &Q);

        pack_point(&in->E[i], &Q, &out->x[i]);

#ifndef NDEBUG
        mont_iso iso;
        canonical_model(&in->E[i], &iso);

        curve_proj E;
        fp2copy(iso.B, E.A);
        fp2one(E.C);

        point_proj K;
        unpack_point(&out->x[i], &K);

        {
            f2elm_t A24plus, C24;
            AC_to_AAC24(&E, A24plus, NULL, C24);
            point_proj Q;
            xDBLe(&K, &Q, A24plus, C24, OALICE_BITS);
            assert(fp2_iszero(Q.Z));
        }

        curve_proj EE;
        bool ret = Isogeny_A(&E, &K, &EE, NULL, NULL);
        assert(ret); (void) ret;

        f2elm_t j0, j1;
        j_inv(EE.A, EE.C, j0);
        j_inv(in->E[i+1].A, in->E[i+1].C, j1);
        assert(fp2_equal(j0, j1));
#endif
    }
}

static void Pack_Chain_B(chain_B const *in, packed_chain_B *out)
{
    for (size_t i = 0; i < POIK_HEIGHT; ++i) {
        pack_point(&in->E[i], &in->P[i], &out->x[i]);

#ifndef NDEBUG
        mont_iso iso;
        canonical_model(&in->E[i], &iso);

        curve_proj E;
        fp2copy(iso.B, E.A);
        fp2one(E.C);

        point_proj K;
        unpack_point(&out->x[i], &K);

        {
            f2elm_t A24plus, A24minus;
            AC_to_AAC24(&E, A24plus, A24minus, NULL);
            point_proj Q;
            xTPLe(&K, &Q, A24minus, A24plus, OBOB_EXPON);
            assert(fp2_iszero(Q.Z));
        }

        curve_proj EE;
        bool ret = Isogeny_B(&E, &K, &EE, NULL, NULL);
        assert(ret); (void) ret;

        f2elm_t j0, j1;
        j_inv(EE.A, EE.C, j0);
        j_inv(in->E[i+1].A, in->E[i+1].C, j1);
        assert(fp2_equal(j0, j1));
#endif
    }
}

static bool Unpack_Chain_A(packed_fp2 const *E0, packed_chain_A const *in, curve_proj *E)
{
    unpack_curve(E0, E, true);

    point_proj back;

    for (size_t i = 0; i < POIK_WIDTH; ++i) {
        if (i) {
            mont_iso iso;
            canonical_model(E, &iso);
            fp2copy(iso.B, E->A);
            fp2one(E->C);
            eval_iso(&iso, &back, &back);
        }

        point_proj K;
        unpack_point(&in->x[i], &K);
        if (!Isogeny_A(E, &K, E, i ? &back : NULL, i ? NULL : &back))
            return false;
    }

    return !fp2_iszero(back.Z);
}

static bool Unpack_Chain_B(packed_fp2 const *E0, packed_chain_B const *in, curve_proj *E)
{
    unpack_curve(E0, E, true);

    point_proj back = {0};

    for (size_t i = 0; i < POIK_HEIGHT; ++i) {
        if (i) {
            mont_iso iso;
            canonical_model(E, &iso);
            fp2copy(iso.B, E->A);
            fp2one(E->C);
            eval_iso(&iso, &back, &back);
        }

        point_proj K;
        unpack_point(&in->x[i], &K);
        if (!Isogeny_B(E, &K, E, i ? &back : NULL, i ? NULL : &back))
            return false;
    }

    return !fp2_iszero(back.Z);
}
////////////////////////////////////////////////////////////////

static void enc_j(curve_proj const *E, packed_fp2 *j)
{
    f2elm_t jinv;
    j_inv(E->A, E->C, jinv);
    fp2_encode(jinv, j->bytes);
}

static void hash_hiding(unsigned char const r[HIDE_BYTES], packed_fp2 const *v, unsigned char h[HASH_BYTES])
{
    unsigned char buf[strlen(RO_COMMIT) + HIDE_BYTES + sizeof(*v)];
    memcpy(buf, RO_COMMIT, strlen(RO_COMMIT));
    memcpy(buf + strlen(RO_COMMIT), r, HIDE_BYTES);
    memcpy(buf + strlen(RO_COMMIT) + HIDE_BYTES, v, sizeof(*v));
    shake256(h, HASH_BYTES, buf, sizeof(buf));
}

static void Instance(struct secret const *sec, struct instance *v, struct commitment1 *comm)
{
    chain_A bottom;
    for (size_t i = 0; i <= POIK_WIDTH; ++i)
        unpack_curve(&sec->E[i], &bottom.E[i], false);
    for (size_t i = 0; i <  POIK_WIDTH; ++i)
        unpack_point(&sec->P[i], &bottom.P[i]);

    chain_B left, right;
    Random_Isogeny_Chain_B(&bottom.E[0], &left);

    for (unsigned row = 0; row < POIK_HEIGHT; ++row) {
        curve_proj El = left.E[row];
        point_proj Kl = left.P[row];
        for (unsigned col = 0; col < POIK_WIDTH; ++col) {
            curve_proj Et = bottom.E[col];
            point_proj Kt = bottom.P[col];
            assert(curve_equal(&Et, &El));
            curve_proj Eb, Er;
            point_proj Kb = Kt, Kr = Kl;
            Isogeny_A(&Et, &Kt, &Er, &Kr, NULL);
            Isogeny_B(&El, &Kl, &Eb, &Kb, NULL);
            bottom.E[col] = Eb;
            bottom.P[col] = Kb;
            El = Er;
            Kl = Kr;
        }
        right.E[row] = El;
        right.P[row] = Kl;
    }

    Isogeny_A(&bottom.E[POIK_WIDTH-1], &bottom.P[POIK_WIDTH-1], &bottom.E[POIK_WIDTH], NULL, NULL);
    Isogeny_B(&right.E[POIK_HEIGHT-1], &right.P[POIK_HEIGHT-1], &right.E[POIK_HEIGHT], NULL, NULL);
    assert(curve_equal(&bottom.E[POIK_WIDTH], &right.E[POIK_HEIGHT]));

    // proof instance
    v->E2 = left.E[POIK_HEIGHT];
    v->E3 = right.E[POIK_HEIGHT];
    v->left = left;
    v->bottom = bottom;
    v->right = right;

    // commitment
    packed_fp2 j2, j3;
    enc_j(&v->E2, &j2);
    enc_j(&v->E3, &j3);
    randombytes_strict(v->r2, sizeof(v->r2));
    randombytes_strict(v->r3, sizeof(v->r3));
    hash_hiding(v->r2, &j2, comm->h2);
    hash_hiding(v->r3, &j3, comm->h3);
}

static void Challenge(struct commitment const *comm, int8_t challs[POIK_REPS])
{
    uint64_t state[25] = {0};
    {
        char buf[strlen(RO_CHALLENGE) + sizeof(*comm)];
        memcpy(buf, RO_CHALLENGE, strlen(RO_CHALLENGE));
        memcpy(buf + strlen(RO_CHALLENGE), comm, sizeof(*comm));
        shake256_absorb(state, (void *) buf, sizeof(buf));
    }

    for (size_t done = 0, n; done < POIK_REPS; done += n) {
        unsigned char blk[SHAKE256_RATE];
        shake256_squeezeblocks(blk, 1, state);
        n = sample_ternary(blk, sizeof(blk), challs + done, POIK_REPS - done);
    }
}

//TODO could batch inversions between multiple proof instances
static void Pack_Response(int8_t which, struct instance const *v, struct commitment1 *c, struct response *r)
{
    switch (r->tag = which) {

    case RESP_LEFT:
        memcpy(r->left.r2, v->r2, HIDE_BYTES);
        memcpy(r->left.h3, c->h3, HASH_BYTES);
        Pack_Chain_B(&v->left, &r->left.chain);
        break;

    case RESP_BOTTOM:
        memcpy(r->bottom.r2, v->r2, HIDE_BYTES);
        memcpy(r->bottom.r3, v->r3, HIDE_BYTES);
        pack_curve(&v->E2, &r->bottom.E2);
        Pack_Chain_A(&v->bottom, &r->bottom.chain);
        break;

    case RESP_RIGHT:
        memcpy(r->right.h2, c->h2, HASH_BYTES);
        memcpy(r->right.r3, v->r3, HIDE_BYTES);
        Pack_Chain_B(&v->right, &r->right.chain);
        break;

    default:
        __builtin_unreachable();
    }
}

static bool Unpack_Response(struct proof const *prf, struct response const *r, int8_t *which, struct commitment1 *comm)
{
    curve_proj E2, E3;
    packed_fp2 j2, j3;

    switch (*which = r->tag) {

    case RESP_LEFT:
        if (!Unpack_Chain_B(&prf->E0, &r->left.chain, &E2))
            return false;
        enc_j(&E2, &j2);
        hash_hiding(r->left.r2, &j2, comm->h2);
        memcpy(comm->h3, r->left.h3, HASH_BYTES);
        break;

    case RESP_BOTTOM:
        if (!unpack_curve(&r->bottom.E2, &E2, true))
            return false;
        if (!Unpack_Chain_A(&r->bottom.E2, &r->bottom.chain, &E3))
            return false;
        enc_j(&E2, &j2);
        enc_j(&E3, &j3);
        hash_hiding(r->bottom.r2, &j2, comm->h2);
        hash_hiding(r->bottom.r3, &j3, comm->h3);
        break;

    case RESP_RIGHT:
        if (!Unpack_Chain_B(&prf->E1, &r->right.chain, &E3))
            return false;
        memcpy(comm->h2, r->right.h2, HASH_BYTES);
        enc_j(&E3, &j3);
        hash_hiding(r->right.r3, &j3, comm->h3);
        break;

    default:
        return false;
    }

    return true;
}

////////////////////////////////////////////////////////////////

FILE *Progress_Indicator_File = NULL;

static void render_progress(size_t idx)
{
    if (!Progress_Indicator_File)
        return;

    if (idx == -1) {
        fputs("\r\x1b[K", Progress_Indicator_File);
        fflush(Progress_Indicator_File);
        return;
    }

    static size_t last = -1;
    unsigned const done = 63 * idx / POIK_REPS;
    if (done == last)
        return;

    char buf[65] = {0};
    memset(buf, '.', 64);
    memset(buf, '=', done);
    buf[done] = '>';
    fprintf(Progress_Indicator_File, "\r\x1b[K%s", buf);
    fflush(Progress_Indicator_File);
    last = done;
}

////////////////////////////////////////////////////////////////

#ifdef THREADING

#define num_threads() 16    // default

#if defined __has_include
    #if __has_include(<unistd.h>)
        #include <unistd.h>
        #undef num_threads
        #define num_threads() sysconf(_SC_NPROCESSORS_ONLN)
    #endif
#endif

//TODO determine number of cores on other systems

#include <pthread.h>
#include <stdatomic.h>

struct thread_state
{
    atomic_size_t idx;
    union {
        struct thread_state_Instance {
            struct secret const *sec;
            struct instance *inst;
            struct commitment1 *comm1;
        } Instance;
        struct thread_state_Pack_Response {
            int8_t *chall;
            struct instance const *inst;
            struct commitment1 *comm1;
            struct response *resp;
        } Pack_Response;
        struct thread_state_Unpack_Response {
            struct proof const *prf;
            int8_t *tags;
            struct commitment1 *comm1;
            atomic_bool *okay;
        } Unpack_Response;
    };
};

static void run_threads(void *(*fun)(void *), struct thread_state *state)
{
    static size_t count = 0;
    if (!count)
        count = num_threads();

    pthread_t threads[count];

    for (size_t i = 0; i < sizeof(threads) / sizeof(*threads); ++i)
        if (pthread_create(&threads[i], NULL, fun, state))
            fail("pthread_create() failed -- could not start worker");

    if (Progress_Indicator_File) {
        for (size_t idx; (idx = atomic_load(&state->idx)) < POIK_REPS; ) {
            render_progress(idx);
            usleep(1000);
        }
        render_progress(-1);
    }

    for (size_t i = 0; i < sizeof(threads) / sizeof(*threads); ++i)
        if (pthread_join(threads[i], NULL))
            fail("pthread_join() failed -- could not wait for worker");
}

static void *Worker_Instance(void *ptr)
{
    struct thread_state *const state = ptr;
    struct thread_state_Instance *const data = &state->Instance;

    while (true) {
        size_t i = atomic_fetch_add(&state->idx, 1);
        if (i >= POIK_REPS)
            return 0;

        Instance(data->sec, &data->inst[i], &data->comm1[i]);
    }
}

static void *Worker_Pack_Response(void *ptr)
{
    struct thread_state *const state = ptr;
    struct thread_state_Pack_Response *const data = &state->Pack_Response;

    while (true) {
        size_t i = atomic_fetch_add(&state->idx, 1);
        if (i >= POIK_REPS)
            return 0;

        Pack_Response(data->chall[i], &data->inst[i], &data->comm1[i], &data->resp[i]);
    }
}

static void *Worker_Unpack_Response(void *ptr)
{
    struct thread_state *const state = ptr;
    struct thread_state_Unpack_Response *const data = &state->Unpack_Response;

    while (true) {
        size_t i = atomic_fetch_add(&state->idx, 1);
        if (i >= POIK_REPS)
            return 0;

        bool ok = Unpack_Response(data->prf, &data->prf->resp[i], &data->tags[i], &data->comm1[i]);
        if (!ok) {
            atomic_store(data->okay, false);
            atomic_store(&state->idx, POIK_REPS);  // signal other threads to stop
            return 0;
        }
    }
}

#endif  // THREADING

////////////////////////////////////////////////////////////////

bool Generate(packed_fp2 const *E0, struct secret *sec)
{
    curve_proj E;
    if (!unpack_curve(E0, &E, true))
        return false;

    chain_A chain;
    Random_Isogeny_Chain_A(&E, &chain);

    for (size_t i = 0; i <= POIK_WIDTH; ++i)
        fp2div_enc(chain.E[i].A, chain.E[i].C, &sec->E[i]);
    for (size_t i = 0; i <  POIK_WIDTH; ++i)
        fp2div_enc(chain.P[i].X, chain.P[i].Z, &sec->P[i]);

    return true;
}

void Prove(struct secret const *sec, struct proof *prf)
{
    curve_proj E0, E1;

    unpack_curve(&sec->E[0], &E0, false);
    unpack_curve(&sec->E[POIK_WIDTH], &E1, false);

    // The size of this array gets dangerously close
    // to the default stack size limit on Linux:
    //     struct instance inst[POIK_REPS];
    // Allocate it on the heap instead to be sure.
    struct instance *inst = malloc(sizeof(*inst) * POIK_REPS);

    if (!inst)
        fail("malloc() failed -- out of memory?");

    struct commitment comm;
    enc_j(&E0, &comm.j0);
    enc_j(&E1, &comm.j1);

    // generate PoIK instances & commitments
#ifdef THREADING
    {
        struct thread_state state = {.Instance = {.sec=sec, .inst=inst, .comm1=comm.hs}};
        run_threads(Worker_Instance, &state);
    }
#else
    for (size_t i = 0; i < POIK_REPS; ++i) {
        render_progress(i+1);
        Instance(sec, &inst[i], &comm.hs[i]);
    }
    render_progress(-1);
#endif

    // derive challenge (Fiat-Shamir)
    int8_t chall[POIK_REPS];
    Challenge(&comm, chall);

    // compute proof
    pack_curve(&E0, &prf->E0);
    pack_curve(&E1, &prf->E1);
#ifdef THREADING
    {
        struct thread_state state = {.Pack_Response = {.chall=chall, .inst=inst, .comm1=comm.hs, .resp=prf->resp}};
        run_threads(Worker_Pack_Response, &state);
    }
#else
    for (size_t i = 0; i < POIK_REPS; ++i) {
        render_progress(i+1);
        Pack_Response(chall[i], &inst[i], &comm.hs[i], &prf->resp[i]);
    }
    render_progress(-1);
#endif

    free(inst);
}

//TODO think carefully about input validation
bool Verify(struct proof *prf)
{
    curve_proj E0, E1;

    if (!unpack_curve(&prf->E0, &E0, true))
        return false;

    if (!unpack_curve(&prf->E1, &E1, true))
        return false;

    struct commitment comm;
    enc_j(&E0, &comm.j0);
    enc_j(&E1, &comm.j1);

    int8_t tags[POIK_REPS];

#ifdef THREADING
    {
        atomic_bool okay = true;

        struct thread_state state = {.Unpack_Response = {.prf=prf, .tags=tags, .comm1=comm.hs, .okay=&okay}};
        run_threads(Worker_Unpack_Response, &state);

        if (!okay)
            return false;
    }
#else
    for (size_t i = 0; i < POIK_REPS; ++i)
        if (!Unpack_Response(prf, &prf->resp[i], &tags[i], &comm.hs[i]))
            return false;
#endif

    // re-derive challenge
    int8_t challenge[POIK_REPS];
    Challenge(&comm, challenge);
    return !memcmp(challenge, tags, sizeof(challenge));
}

////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

static bool Dump_Bytes(FILE *file, unsigned char const *buf, size_t len)
{
    for (size_t i = 0; i < len; ++i)
        if (2 != fprintf(file, "%02hhx", buf[i]))
            return false;
    if (EOF == fputc('\n', file))
        return false;
    return true;
}

#define FP_ENCODED_BYTES (FP2_ENCODED_BYTES/2)

static bool dump_packed_fp(FILE *file, unsigned char const *v)
{
    if (2 != fprintf(file, "0x"))
        return false;
    for (ssize_t i = FP_ENCODED_BYTES-1; i >= 0; --i)
        if (2 != fprintf(file, "%02hhx", v[i]))
            return false;
    return true;
}

bool Dump_Fp2(FILE *file, packed_fp2 const *v)
{
    if (!dump_packed_fp(file, v->real))
        return false;
    if (EOF == fputc(',', file))
        return false;
    if (!dump_packed_fp(file, v->imag))
        return false;
    if (EOF == fputc('\n', file))
        return false;
    return true;
}

bool Dump(FILE *file, struct proof const *prf) {
#define dump(v) do { if (!Dump_Fp2(file, (v))) return false; } while (false)
#define dump_bytes(v) do { if (!Dump_Bytes(file, (v), sizeof(v))) return false; } while (false)
    dump(&prf->E0);
    dump(&prf->E1);

    for (size_t i = 0; i < POIK_REPS; ++i) {
        struct response const *const r = &prf->resp[i];

        fprintf(file, "%d\n", r->tag);

        switch (r->tag) {
        case RESP_LEFT:
            dump_bytes(r->left.r2);
            dump_bytes(r->left.h3);
            for (size_t i = 0; i < POIK_HEIGHT; ++i)
                dump(&r->left.chain.x[i]);
            break;
        case RESP_BOTTOM:
            dump_bytes(r->bottom.r2);
            dump_bytes(r->bottom.r3);
            dump(&r->bottom.E2);
            for (size_t i = 0; i < POIK_WIDTH; ++i)
                dump(&r->bottom.chain.x[i]);
            break;
        case RESP_RIGHT:
            dump_bytes(r->right.h2);
            dump_bytes(r->right.r3);
            for (size_t i = 0; i < POIK_HEIGHT; ++i)
                dump(&r->right.chain.x[i]);
            break;
        default:
            __builtin_unreachable();
        }
    }

    return true;
#undef dump
#undef dump_bytes
}

static int parse_hex_byte(FILE *file)
{
    char a,b;
    if (2 != fscanf(file, "%c%c", &a, &b) || !isxdigit(a) || !isxdigit(b))
        return -1;
    char s[3] = {a, b, 0};
    return strtol(s, NULL, 16);
}

static bool Parse_Bytes(FILE *file, unsigned char *buf, size_t len)
{
    fscanf(file, " ");  // ignore whitespace
    for (size_t i = 0; i < len; ++i) {
        int c = parse_hex_byte(file);
        if (c < 0)
            return false;
        buf[i] = c;
    }
    return true;
}

static bool parse_packed_fp(FILE *file, unsigned char *v)
{
    char a,b;
    fscanf(file, " ");  // ignore whitespace
    if (2 != fscanf(file, "%c%c", &a, &b) || a != '0' || b != 'x')
        return false;
    for (ssize_t i = FP_ENCODED_BYTES-1; i >= 0; --i) {
        int c = parse_hex_byte(file);
        if (c < 0)
            return false;
        v[i] = c;
    }
    return true;
}

bool Parse_Fp2(FILE *file, packed_fp2 *v)
{
    if (!parse_packed_fp(file, v->real))
        return false;
    char c;
    if (1 != fscanf(file, " %c", &c) || c != ',')
        return false;
    if (!parse_packed_fp(file, v->imag))
        return false;
    return true;
}

bool Parse(FILE *file, struct proof *prf)
{
#define parse(v) do { if (!Parse_Fp2(file, (v))) return false; } while (false)
#define parse_bytes(v) do { if (!Parse_Bytes(file, (v), sizeof(v))) return false; } while (false)
    parse(&prf->E0);
    parse(&prf->E1);

    for (size_t i = 0; i < POIK_REPS; ++i) {
        struct response *const r = &prf->resp[i];

        if (1 != fscanf(file, " %d", &r->tag))
            return false;

        switch (r->tag) {
        case RESP_LEFT:
            parse_bytes(r->left.r2);
            parse_bytes(r->left.h3);
            for (size_t i = 0; i < POIK_HEIGHT; ++i)
                parse(&r->left.chain.x[i]);
            break;
        case RESP_BOTTOM:
            parse_bytes(r->bottom.r2);
            parse_bytes(r->bottom.r3);
            parse(&r->bottom.E2);
            for (size_t i = 0; i < POIK_WIDTH; ++i)
                parse(&r->bottom.chain.x[i]);
            break;
        case RESP_RIGHT:
            parse_bytes(r->right.h2);
            parse_bytes(r->right.r3);
            for (size_t i = 0; i < POIK_HEIGHT; ++i)
                parse(&r->right.chain.x[i]);
            break;
        default:
            return false;
        }
    }

    return true;
#undef parse
#undef parse_bytes
}

