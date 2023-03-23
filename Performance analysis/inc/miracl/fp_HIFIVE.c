/*
 * Copyright (c) 2012-2020 MIRACL UK Ltd.
 *
 * This file is part of MIRACL Core
 * (see https://github.com/miracl/core).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* CORE mod p functions */
/* Small Finite Field arithmetic */
/* SU=m, SU is Stack Usage (NOT_SPECIAL Modulus) */

#include "fp_HIFIVE.h"

/* Fast Modular Reduction Methods */

/* r=d mod m */
/* d MUST be normalised */
/* Products must be less than pR in all cases !!! */
/* So when multiplying two numbers, their product *must* be less than MODBITS+BASEBITS*NLEN */
/* Results *may* be one bit bigger than MODBITS */

#if MODTYPE_HIFIVE == PSEUDO_MERSENNE
/* r=d mod m */

/* Converts from BIG integer to residue form mod Modulus */
void FP_HIFIVE_nres(FP_HIFIVE *y, BIG_336_60 x)
{
    BIG_336_60 mdls;
    BIG_336_60_rcopy(mdls, Modulus_HIFIVE);
    BIG_336_60_copy(y->g, x);
    BIG_336_60_mod(y->g,mdls);
    y->XES = 1;
}

/* Converts from residue form back to BIG integer form */
void FP_HIFIVE_redc(BIG_336_60 x, FP_HIFIVE *y)
{
    BIG_336_60_copy(x, y->g);
}

/* reduce a DBIG to a BIG exploiting the special form of the modulus */
void FP_HIFIVE_mod(BIG_336_60 r, DBIG_336_60 d)
{
    BIG_336_60 t, b;
    chunk v, tw;
    BIG_336_60_split(t, b, d, MODBITS_HIFIVE);

    /* Note that all of the excess gets pushed into t. So if squaring a value with a 4-bit excess, this results in
       t getting all 8 bits of the excess product! So products must be less than pR which is Montgomery compatible */

    if (MConst_HIFIVE < NEXCESS_336_60)
    {
        BIG_336_60_imul(t, t, MConst_HIFIVE);
        BIG_336_60_norm(t);
        BIG_336_60_add(r, t, b);
        BIG_336_60_norm(r);
        tw = r[NLEN_336_60 - 1];
        r[NLEN_336_60 - 1] &= TMASK_HIFIVE;
        r[0] += MConst_HIFIVE * ((tw >> TBITS_HIFIVE));
    }
    else
    {
        v = BIG_336_60_pmul(t, t, MConst_HIFIVE);
        BIG_336_60_add(r, t, b);
        BIG_336_60_norm(r);
        tw = r[NLEN_336_60 - 1];
        r[NLEN_336_60 - 1] &= TMASK_HIFIVE;
#if CHUNK == 16
        r[1] += muladd_336_60(MConst_HIFIVE, ((tw >> TBITS_HIFIVE) + (v << (BASEBITS_336_60 - TBITS_HIFIVE))), 0, &r[0]);
#else
        r[0] += MConst_HIFIVE * ((tw >> TBITS_HIFIVE) + (v << (BASEBITS_336_60 - TBITS_HIFIVE)));
#endif
    }
    BIG_336_60_norm(r);
}
#endif

/* This only applies to Curve C448, so specialised (for now) */
#if MODTYPE_HIFIVE == GENERALISED_MERSENNE

void FP_HIFIVE_nres(FP_HIFIVE *y, BIG_336_60 x)
{
    BIG_336_60 mdls;
    BIG_336_60_rcopy(mdls, Modulus_HIFIVE);
    BIG_336_60_copy(y->g, x);
    BIG_336_60_mod(y->g,mdls);
    y->XES = 1;
}

/* Converts from residue form back to BIG integer form */
void FP_HIFIVE_redc(BIG_336_60 x, FP_HIFIVE *y)
{
    BIG_336_60_copy(x, y->g);
}

/* reduce a DBIG to a BIG exploiting the special form of the modulus */
void FP_HIFIVE_mod(BIG_336_60 r, DBIG_336_60 d)
{
    BIG_336_60 t, b;
    chunk carry;
    BIG_336_60_split(t, b, d, MBITS_HIFIVE);

    BIG_336_60_add(r, t, b);

    BIG_336_60_dscopy(d, t);
    BIG_336_60_dshl(d, MBITS_HIFIVE / 2);

    BIG_336_60_split(t, b, d, MBITS_HIFIVE);

    BIG_336_60_add(r, r, t);
    BIG_336_60_add(r, r, b);
    BIG_336_60_norm(r);
    BIG_336_60_shl(t, MBITS_HIFIVE / 2);

    BIG_336_60_add(r, r, t);

    carry = r[NLEN_336_60 - 1] >> TBITS_HIFIVE;

    r[NLEN_336_60 - 1] &= TMASK_HIFIVE;
    r[0] += carry;

    r[224 / BASEBITS_336_60] += carry << (224 % BASEBITS_336_60); /* need to check that this falls mid-word */
    BIG_336_60_norm(r);
}

#endif

#if MODTYPE_HIFIVE == MONTGOMERY_FRIENDLY

/* convert to Montgomery n-residue form */
void FP_HIFIVE_nres(FP_HIFIVE *y, BIG_336_60 x)
{
    DBIG_336_60 d;
    BIG_336_60 r;
    BIG_336_60_rcopy(r, R2modp_HIFIVE);
    BIG_336_60_mul(d, x, r);
    FP_HIFIVE_mod(y->g, d);
    y->XES = 2;
}

/* convert back to regular form */
void FP_HIFIVE_redc(BIG_336_60 x, FP_HIFIVE *y)
{
    DBIG_336_60 d;
    BIG_336_60_dzero(d);
    BIG_336_60_dscopy(d, y->g);
    FP_HIFIVE_mod(x, d);
}

/* fast modular reduction from DBIG to BIG exploiting special form of the modulus */
void FP_HIFIVE_mod(BIG_336_60 a, DBIG_336_60 d)
{
    int i;

    for (i = 0; i < NLEN_336_60; i++)
        d[NLEN_336_60 + i] += muladd_336_60(d[i], MConst_HIFIVE - 1, d[i], &d[NLEN_336_60 + i - 1]);

    BIG_336_60_sducopy(a, d);
    BIG_336_60_norm(a);
}

#endif

#if MODTYPE_HIFIVE == NOT_SPECIAL

/* convert to Montgomery n-residue form */
void FP_HIFIVE_nres(FP_HIFIVE *y, BIG_336_60 x)
{
    DBIG_336_60 d;
    BIG_336_60 r;
    BIG_336_60_rcopy(r, R2modp_HIFIVE);
    BIG_336_60_mul(d, x, r);
    FP_HIFIVE_mod(y->g, d);
    y->XES = 2;
}

/* convert back to regular form */
void FP_HIFIVE_redc(BIG_336_60 x, FP_HIFIVE *y)
{
    DBIG_336_60 d;
    BIG_336_60_dzero(d);
    BIG_336_60_dscopy(d, y->g);
    FP_HIFIVE_mod(x, d);
}


/* reduce a DBIG to a BIG using Montgomery's no trial division method */
/* d is expected to be dnormed before entry */
/* SU= 112 */
void FP_HIFIVE_mod(BIG_336_60 a, DBIG_336_60 d)
{
    BIG_336_60 mdls;
    BIG_336_60_rcopy(mdls, Modulus_HIFIVE);
    BIG_336_60_monty(a, mdls, MConst_HIFIVE, d);
}

#endif

void FP_HIFIVE_from_int(FP_HIFIVE *x,int a)
{
    BIG_336_60 w;
    if (a<0) BIG_336_60_rcopy(w, Modulus_HIFIVE);
    else BIG_336_60_zero(w); 
    BIG_336_60_inc(w,a); BIG_336_60_norm(w); 
    FP_HIFIVE_nres(x,w);
}

/* test x==0 ? */
/* SU= 48 */
int FP_HIFIVE_iszilch(FP_HIFIVE *x)
{
    BIG_336_60 m;
    FP_HIFIVE y;
    FP_HIFIVE_copy(&y,x);
    FP_HIFIVE_reduce(&y);
    FP_HIFIVE_redc(m,&y);
    return BIG_336_60_iszilch(m);
}

int FP_HIFIVE_isunity(FP_HIFIVE *x)
{
    BIG_336_60 m;
    FP_HIFIVE y;
    FP_HIFIVE_copy(&y,x);
    FP_HIFIVE_reduce(&y);
    FP_HIFIVE_redc(m,&y);
    return BIG_336_60_isunity(m);
}


void FP_HIFIVE_copy(FP_HIFIVE *y, FP_HIFIVE *x)
{
    BIG_336_60_copy(y->g, x->g);
    y->XES = x->XES;
}

void FP_HIFIVE_rcopy(FP_HIFIVE *y, const BIG_336_60 c)
{
    BIG_336_60 b;
    BIG_336_60_rcopy(b, c);
    FP_HIFIVE_nres(y, b);
}

/* Swap a and b if d=1 */
void FP_HIFIVE_cswap(FP_HIFIVE *a, FP_HIFIVE *b, int d)
{
    sign32 t, c = d;
    BIG_336_60_cswap(a->g, b->g, d);

    c = ~(c - 1);
    t = c & ((a->XES) ^ (b->XES));
    a->XES ^= t;
    b->XES ^= t;

}

/* Move b to a if d=1 */
void FP_HIFIVE_cmove(FP_HIFIVE *a, FP_HIFIVE *b, int d)
{
    sign32 c = -d;

    BIG_336_60_cmove(a->g, b->g, d);
    a->XES ^= (a->XES ^ b->XES)&c;
}

void FP_HIFIVE_zero(FP_HIFIVE *x)
{
    BIG_336_60_zero(x->g);
    x->XES = 1;
}

int FP_HIFIVE_equals(FP_HIFIVE *x, FP_HIFIVE *y)
{
    FP_HIFIVE xg, yg;
    FP_HIFIVE_copy(&xg, x);
    FP_HIFIVE_copy(&yg, y);
    FP_HIFIVE_reduce(&xg);
    FP_HIFIVE_reduce(&yg);
    if (BIG_336_60_comp(xg.g, yg.g) == 0) return 1;
    return 0;
}

// Is x lexically larger than p-x?
// return -1 for no, 0 if x=0, 1 for yes
int FP_HIFIVE_islarger(FP_HIFIVE *x)
{
    BIG_336_60 p,fx,sx;
    if (FP_HIFIVE_iszilch(x)) return 0;
    BIG_336_60_rcopy(p,Modulus_HIFIVE);
    FP_HIFIVE_redc(fx,x);
    BIG_336_60_sub(sx,p,fx);  BIG_336_60_norm(sx); 
    return BIG_336_60_comp(fx,sx);
}

void FP_HIFIVE_toBytes(char *b,FP_HIFIVE *x)
{
    BIG_336_60 t;
    FP_HIFIVE_redc(t, x);
    BIG_336_60_toBytes(b, t);
}

void FP_HIFIVE_fromBytes(FP_HIFIVE *x,char *b)
{
    BIG_336_60 t;
    BIG_336_60_fromBytes(t, b);
    FP_HIFIVE_nres(x, t);
}

/* output FP */
/* SU= 48 */
void FP_HIFIVE_output(FP_HIFIVE *r)
{
    BIG_336_60 c;
    FP_HIFIVE_reduce(r);
    FP_HIFIVE_redc(c, r);
    BIG_336_60_output(c);
}

void FP_HIFIVE_rawoutput(FP_HIFIVE *r)
{
    BIG_336_60_rawoutput(r->g);
}

#ifdef GET_STATS
int tsqr = 0, rsqr = 0, tmul = 0, rmul = 0;
int tadd = 0, radd = 0, tneg = 0, rneg = 0;
int tdadd = 0, rdadd = 0, tdneg = 0, rdneg = 0;
#endif

#ifdef FUSED_MODMUL

/* Insert fastest code here */

#endif

/* r=a*b mod Modulus */
/* product must be less that p.R - and we need to know this in advance! */
/* SU= 88 */
void FP_HIFIVE_mul(FP_HIFIVE *r, FP_HIFIVE *a, FP_HIFIVE *b)
{
    DBIG_336_60 d;

    if ((sign64)a->XES * b->XES > (sign64)FEXCESS_HIFIVE)
    {
#ifdef DEBUG_REDUCE
        printf("Product too large - reducing it\n");
#endif
        FP_HIFIVE_reduce(a);  /* it is sufficient to fully reduce just one of them < p */
    }

#ifdef FUSED_MODMUL
    FP_HIFIVE_modmul(r->g, a->g, b->g);
#else
    BIG_336_60_mul(d, a->g, b->g);
    FP_HIFIVE_mod(r->g, d);
#endif
    r->XES = 2;
}


/* multiplication by an integer, r=a*c */
/* SU= 136 */
void FP_HIFIVE_imul(FP_HIFIVE *r, FP_HIFIVE *a, int c)
{
    int s = 0;

    if (c < 0)
    {
        c = -c;
        s = 1;
    }

#if MODTYPE_HIFIVE==PSEUDO_MERSENNE || MODTYPE_HIFIVE==GENERALISED_MERSENNE
    DBIG_336_60 d;
    BIG_336_60_pxmul(d, a->g, c);
    FP_HIFIVE_mod(r->g, d);
    r->XES = 2;

#else
    //Montgomery
    BIG_336_60 k;
    FP_HIFIVE f;
    if (a->XES * c <= FEXCESS_HIFIVE)
    {
        BIG_336_60_pmul(r->g, a->g, c);
        r->XES = a->XES * c; // careful here - XES jumps!
    }
    else
    {
        // don't want to do this - only a problem for Montgomery modulus and larger constants
        BIG_336_60_zero(k);
        BIG_336_60_inc(k, c);
        BIG_336_60_norm(k);
        FP_HIFIVE_nres(&f, k);
        FP_HIFIVE_mul(r, a, &f);
    }
#endif

    if (s)
    {
        FP_HIFIVE_neg(r, r);
        FP_HIFIVE_norm(r);
    }
}

/* Set r=a^2 mod m */
/* SU= 88 */
void FP_HIFIVE_sqr(FP_HIFIVE *r, FP_HIFIVE *a)
{
    DBIG_336_60 d;

    if ((sign64)a->XES * a->XES > (sign64)FEXCESS_HIFIVE)
    {
#ifdef DEBUG_REDUCE
        printf("Product too large - reducing it\n");
#endif
        FP_HIFIVE_reduce(a);
    }

    BIG_336_60_sqr(d, a->g);
    FP_HIFIVE_mod(r->g, d);
    r->XES = 2;
}

/* SU= 16 */
/* Set r=a+b */
void FP_HIFIVE_add(FP_HIFIVE *r, FP_HIFIVE *a, FP_HIFIVE *b)
{
    BIG_336_60_add(r->g, a->g, b->g);
    r->XES = a->XES + b->XES;
    if (r->XES > FEXCESS_HIFIVE)
    {
#ifdef DEBUG_REDUCE
        printf("Sum too large - reducing it \n");
#endif
        FP_HIFIVE_reduce(r);
    }
}

/* Set r=a-b mod m */
/* SU= 56 */
void FP_HIFIVE_sub(FP_HIFIVE *r, FP_HIFIVE *a, FP_HIFIVE *b)
{
    FP_HIFIVE n;
    FP_HIFIVE_neg(&n, b);
    FP_HIFIVE_add(r, a, &n);
}

// https://graphics.stanford.edu/~seander/bithacks.html
// constant time log to base 2 (or number of bits in)

static int logb2(unsign32 v)
{
    int r;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;

    v = v - ((v >> 1) & 0x55555555);
    v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
    r = (((v + (v >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
    return r;
}

// find appoximation to quotient of a/m
// Out by at most 2.
// Note that MAXXES is bounded to be 2-bits less than half a word
static int quo(BIG_336_60 n, BIG_336_60 m)
{
    int sh;
    chunk num, den;
    int hb = CHUNK / 2;
    if (TBITS_HIFIVE < hb)
    {
        sh = hb - TBITS_HIFIVE;
        num = (n[NLEN_336_60 - 1] << sh) | (n[NLEN_336_60 - 2] >> (BASEBITS_336_60 - sh));
        den = (m[NLEN_336_60 - 1] << sh) | (m[NLEN_336_60 - 2] >> (BASEBITS_336_60 - sh));
    }
    else
    {
        num = n[NLEN_336_60 - 1];
        den = m[NLEN_336_60 - 1];
    }
    return (int)(num / (den + 1));
}

/* SU= 48 */
/* Fully reduce a mod Modulus */
void FP_HIFIVE_reduce(FP_HIFIVE *a)
{
    BIG_336_60 m, r;
    int sr, sb, q;
    chunk carry;

    BIG_336_60_rcopy(m, Modulus_HIFIVE);

    BIG_336_60_norm(a->g);

    if (a->XES > 16)
    {
        q = quo(a->g, m);
        carry = BIG_336_60_pmul(r, m, q);
        r[NLEN_336_60 - 1] += (carry << BASEBITS_336_60); // correction - put any carry out back in again
        BIG_336_60_sub(a->g, a->g, r);
        BIG_336_60_norm(a->g);
        sb = 2;
    }
    else sb = logb2(a->XES - 1); // sb does not depend on the actual data

    BIG_336_60_fshl(m, sb);

    while (sb > 0)
    {
// constant time...
        sr = BIG_336_60_ssn(r, a->g, m); // optimized combined shift, subtract and norm
        BIG_336_60_cmove(a->g, r, 1 - sr);
        sb--;
    }

    //BIG_336_60_mod(a->g,m);
    a->XES = 1;
}

void FP_HIFIVE_norm(FP_HIFIVE *x)
{
    BIG_336_60_norm(x->g);
}

/* Set r=-a mod Modulus */
/* SU= 64 */
void FP_HIFIVE_neg(FP_HIFIVE *r, FP_HIFIVE *a)
{
    int sb;
    BIG_336_60 m;

    BIG_336_60_rcopy(m, Modulus_HIFIVE);

    sb = logb2(a->XES - 1);
    BIG_336_60_fshl(m, sb);
    BIG_336_60_sub(r->g, m, a->g);
    r->XES = ((sign32)1 << sb) + 1;

    if (r->XES > FEXCESS_HIFIVE)
    {
#ifdef DEBUG_REDUCE
        printf("Negation too large -  reducing it \n");
#endif
        FP_HIFIVE_reduce(r);
    }

}

/* Set r=a/2. */
/* SU= 56 */
void FP_HIFIVE_div2(FP_HIFIVE *r, FP_HIFIVE *a)
{
    BIG_336_60 m;
    BIG_336_60 w;
    BIG_336_60_rcopy(m, Modulus_HIFIVE);
    int pr=BIG_336_60_parity(a->g);

    FP_HIFIVE_copy(r, a);
    BIG_336_60_copy(w,r->g);
    BIG_336_60_fshr(r->g,1);
    BIG_336_60_add(w, w, m);
    BIG_336_60_norm(w);
    BIG_336_60_fshr(w, 1);   
    
    BIG_336_60_cmove(r->g,w,pr);

}

// Could leak size of b
// but not used here with secret exponent b
void FP_HIFIVE_pow(FP_HIFIVE *r, FP_HIFIVE *a, BIG_336_60 b)
{
    sign8 w[1 + (NLEN_336_60 * BASEBITS_336_60 + 3) / 4];
    FP_HIFIVE tb[16];
    BIG_336_60 t;
    int i, nb;

    FP_HIFIVE_copy(r,a);
    FP_HIFIVE_norm(r);
    BIG_336_60_copy(t, b);
    BIG_336_60_norm(t);
    nb = 1 + (BIG_336_60_nbits(t) + 3) / 4;
    /* convert exponent to 4-bit window */
    for (i = 0; i < nb; i++)
    {
        w[i] = BIG_336_60_lastbits(t, 4);
        BIG_336_60_dec(t, w[i]);
        BIG_336_60_norm(t);
        BIG_336_60_fshr(t, 4);
    }

    FP_HIFIVE_one(&tb[0]);
    FP_HIFIVE_copy(&tb[1], r);
    for (i = 2; i < 16; i++)
        FP_HIFIVE_mul(&tb[i], &tb[i - 1], r);

    FP_HIFIVE_copy(r, &tb[w[nb - 1]]);
    for (i = nb - 2; i >= 0; i--)
    {
        FP_HIFIVE_sqr(r, r);
        FP_HIFIVE_sqr(r, r);
        FP_HIFIVE_sqr(r, r);
        FP_HIFIVE_sqr(r, r);
        FP_HIFIVE_mul(r, r, &tb[w[i]]);
    }
    FP_HIFIVE_reduce(r);
}


#if MODTYPE_HIFIVE == PSEUDO_MERSENNE || MODTYPE_HIFIVE==GENERALISED_MERSENNE

// See eprint paper https://eprint.iacr.org/2018/1038
// If p=3 mod 4 r= x^{(p-3)/4}, if p=5 mod 8 r=x^{(p-5)/8}

static void FP_HIFIVE_fpow(FP_HIFIVE *r, FP_HIFIVE *x)
{
    int i, j, k, bw, w, c, nw, lo, m, n, nd, e=PM1D2_HIFIVE;
    FP_HIFIVE xp[11], t, key;
    const int ac[] = {1, 2, 3, 6, 12, 15, 30, 60, 120, 240, 255};
// phase 1
    FP_HIFIVE_copy(&xp[0], x); // 1
    FP_HIFIVE_sqr(&xp[1], x); // 2
    FP_HIFIVE_mul(&xp[2], &xp[1], x); //3
    FP_HIFIVE_sqr(&xp[3], &xp[2]); // 6
    FP_HIFIVE_sqr(&xp[4], &xp[3]); // 12
    FP_HIFIVE_mul(&xp[5], &xp[4], &xp[2]); // 15
    FP_HIFIVE_sqr(&xp[6], &xp[5]); // 30
    FP_HIFIVE_sqr(&xp[7], &xp[6]); // 60
    FP_HIFIVE_sqr(&xp[8], &xp[7]); // 120
    FP_HIFIVE_sqr(&xp[9], &xp[8]); // 240
    FP_HIFIVE_mul(&xp[10], &xp[9], &xp[5]); // 255

#if MODTYPE_HIFIVE==PSEUDO_MERSENNE
    n = MODBITS_HIFIVE;
#endif
#if MODTYPE_HIFIVE==GENERALISED_MERSENNE  // Goldilocks ONLY
    n = MODBITS_HIFIVE / 2;
#endif

    n-=(e+1);
    c=(MConst_HIFIVE+(1<<e)+1)/(1<<(e+1));

// need c to be odd
    nd=0;
    while (c%2==0)
    {
        c/=2;
        n-=1;
        nd++;
    }

    bw = 0; w = 1; while (w < c) {w *= 2; bw += 1;}
    k = w - c;

    if (k != 0)
    {
        i = 10; while (ac[i] > k) i--;
        FP_HIFIVE_copy(&key, &xp[i]);
        k -= ac[i];
    }
    while (k != 0)
    {
        i--;
        if (ac[i] > k) continue;
        FP_HIFIVE_mul(&key, &key, &xp[i]);
        k -= ac[i];
    }

// phase 2
    FP_HIFIVE_copy(&xp[1], &xp[2]);
    FP_HIFIVE_copy(&xp[2], &xp[5]);
    FP_HIFIVE_copy(&xp[3], &xp[10]);

    j = 3; m = 8;
    nw = n - bw;
    while (2 * m < nw)
    {
        FP_HIFIVE_copy(&t, &xp[j++]);
        for (i = 0; i < m; i++)
            FP_HIFIVE_sqr(&t, &t);
        FP_HIFIVE_mul(&xp[j], &xp[j - 1], &t);
        m *= 2;
    }

    lo = nw - m;
    FP_HIFIVE_copy(r, &xp[j]);

    while (lo != 0)
    {
        m /= 2; j--;
        if (lo < m) continue;
        lo -= m;
        FP_HIFIVE_copy(&t, r);
        for (i = 0; i < m; i++)
            FP_HIFIVE_sqr(&t, &t);
        FP_HIFIVE_mul(r, &t, &xp[j]);
    }
// phase 3

    if (bw != 0)
    {
        for (i = 0; i < bw; i++ )
            FP_HIFIVE_sqr(r, r);
        FP_HIFIVE_mul(r, r, &key);
    }
#if MODTYPE_HIFIVE==GENERALISED_MERSENNE  // Goldilocks ONLY
    FP_HIFIVE_copy(&key, r);
    FP_HIFIVE_sqr(&t, &key);
    FP_HIFIVE_mul(r, &t, &xp[0]);
    for (i = 0; i < n + 1; i++)
        FP_HIFIVE_sqr(r, r);
    FP_HIFIVE_mul(r, r, &key);
#endif

    for (i=0;i<nd;i++)
        FP_HIFIVE_sqr(r,r);
}

#endif


// calculates r=x^(p-1-2^e)/2^{e+1) where 2^e|p-1
void FP_HIFIVE_progen(FP_HIFIVE *r,FP_HIFIVE *x)
{
#if MODTYPE_HIFIVE==PSEUDO_MERSENNE  || MODTYPE_HIFIVE==GENERALISED_MERSENNE
    FP_HIFIVE_fpow(r, x);  
#else
    int e=PM1D2_HIFIVE;
    BIG_336_60 m;
    BIG_336_60_rcopy(m, Modulus_HIFIVE);
    BIG_336_60_dec(m,1);
    BIG_336_60_shr(m,e);
    BIG_336_60_dec(m,1);
    BIG_336_60_fshr(m,1);
    FP_HIFIVE_pow(r,x,m);
#endif
}

/* Is x a QR? return optional hint for fast follow-up square root */
int FP_HIFIVE_qr(FP_HIFIVE *x,FP_HIFIVE *h)
{
    FP_HIFIVE r;
    int i,e=PM1D2_HIFIVE;
    FP_HIFIVE_progen(&r,x);
    if (h!=NULL)
        FP_HIFIVE_copy(h,&r);

    FP_HIFIVE_sqr(&r,&r);
    FP_HIFIVE_mul(&r,x,&r);
    for (i=0;i<e-1;i++ )
        FP_HIFIVE_sqr(&r,&r);


//    for (i=0;i<e;i++)
//        FP_HIFIVE_sqr(&r,&r);
//    FP_HIFIVE_copy(&s,x);
//    for (i=0;i<e-1;i++ )
//        FP_HIFIVE_sqr(&s,&s);
//    FP_HIFIVE_mul(&r,&r,&s);
    
    return FP_HIFIVE_isunity(&r);
}

/* Modular inversion */
void FP_HIFIVE_inv(FP_HIFIVE *r,FP_HIFIVE *x,FP_HIFIVE *h)
{
    int i,e=PM1D2_HIFIVE;
    FP_HIFIVE s,t;
    FP_HIFIVE_norm(x);
    FP_HIFIVE_copy(&s,x);

    if (h==NULL)
        FP_HIFIVE_progen(&t,x);
    else
        FP_HIFIVE_copy(&t,h);

    for (i=0;i<e-1;i++)
    {  
        FP_HIFIVE_sqr(&s,&s);
        FP_HIFIVE_mul(&s,&s,x);
    }
  
    for (i=0;i<=e;i++)
        FP_HIFIVE_sqr(&t,&t);
    
    FP_HIFIVE_mul(r,&t,&s);
    FP_HIFIVE_reduce(r);
}

// Tonelli-Shanks in constant time
void FP_HIFIVE_sqrt(FP_HIFIVE *r, FP_HIFIVE *a, FP_HIFIVE *h)
{
    int i,j,k,u,e=PM1D2_HIFIVE;
    FP_HIFIVE v,g,t,b;
    BIG_336_60 m;

    if (h==NULL)
        FP_HIFIVE_progen(&g,a);
    else
        FP_HIFIVE_copy(&g,h);

    BIG_336_60_rcopy(m,ROI_HIFIVE);
    FP_HIFIVE_nres(&v,m);

    FP_HIFIVE_sqr(&t,&g);
    FP_HIFIVE_mul(&t,&t,a);
   
    FP_HIFIVE_mul(r,&g,a);
    FP_HIFIVE_copy(&b,&t);
    for (k=e;k>1;k--)
    {
        for (j=1;j<k-1;j++)
            FP_HIFIVE_sqr(&b,&b);
        u=1-FP_HIFIVE_isunity(&b);
        FP_HIFIVE_mul(&g,r,&v);
        FP_HIFIVE_cmove(r,&g,u);
        FP_HIFIVE_sqr(&v,&v);
        FP_HIFIVE_mul(&g,&t,&v);
        FP_HIFIVE_cmove(&t,&g,u);
        FP_HIFIVE_copy(&b,&t);
    }
// always return +ve square root
    k=FP_HIFIVE_sign(r);
    FP_HIFIVE_neg(&v,r); FP_HIFIVE_norm(&v);
    FP_HIFIVE_cmove(r,&v,k);
}

// Calculate both inverse and square root of x, return QR
int FP_HIFIVE_invsqrt(FP_HIFIVE *i, FP_HIFIVE *s, FP_HIFIVE *x)
{
    FP_HIFIVE h;
    int qr=FP_HIFIVE_qr(x,&h);
    FP_HIFIVE_sqrt(s,x,&h);
    FP_HIFIVE_inv(i,x,&h);
    return qr;
}

// Two for Price of One - See Hamburg https://eprint.iacr.org/2012/309.pdf
// Calculate inverse of i and square root of s, return QR
int FP_HIFIVE_tpo(FP_HIFIVE *i, FP_HIFIVE *s)
{
    int qr;
    FP_HIFIVE w,t;
    FP_HIFIVE_mul(&w,s,i);
    FP_HIFIVE_mul(&t,&w,i);
    qr=FP_HIFIVE_invsqrt(i,s,&t);
    FP_HIFIVE_mul(i,i,&w);
    FP_HIFIVE_mul(s,s,i);
    return qr;
}

/* SU=8 */
/* set n=1 */
void FP_HIFIVE_one(FP_HIFIVE *n)
{
    BIG_336_60 b;
    BIG_336_60_one(b);
    FP_HIFIVE_nres(n, b);
}

int FP_HIFIVE_sign(FP_HIFIVE *x)
{
#ifdef BIG_ENDIAN_SIGN_HIFIVE
    int cp;
    BIG_336_60 m,pm1d2;
    FP_HIFIVE y;
    BIG_336_60_rcopy(pm1d2, Modulus_HIFIVE);
    BIG_336_60_dec(pm1d2,1);
    BIG_336_60_fshr(pm1d2,1); //(p-1)/2
     
    FP_HIFIVE_copy(&y,x);
    FP_HIFIVE_reduce(&y);
    FP_HIFIVE_redc(m,&y);
    cp=BIG_336_60_comp(m,pm1d2);
    return ((cp+1)&2)>>1;

#else
    BIG_336_60 m;
    FP_HIFIVE y;
    FP_HIFIVE_copy(&y,x);
    FP_HIFIVE_reduce(&y);
    FP_HIFIVE_redc(m,&y);
    return BIG_336_60_parity(m);
#endif
}

void FP_HIFIVE_rand(FP_HIFIVE *x,csprng *rng)
{
    BIG_336_60 w,m;
    BIG_336_60_rcopy(m,Modulus_HIFIVE);
    BIG_336_60_randomnum(w,m,rng);
    FP_HIFIVE_nres(x,w);
}


