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

/**
 * @file fp.h
 * @author Mike Scott
 * @brief FP Header File
 *
 */

#ifndef FP_BRAINPOOL_H
#define FP_BRAINPOOL_H

#include "big_256_56.h"
#include "config_field_BRAINPOOL.h"


/**
	@brief FP Structure - quadratic extension field
*/

typedef struct
{
    BIG_256_56 g;	/**< Big representation of field element */
    sign32 XES;	/**< Excess */
} FP_BRAINPOOL;


/* Field Params - see rom.c */
extern const BIG_256_56 Modulus_BRAINPOOL;	/**< Actual Modulus set in rom_field_yyy.c */
extern const BIG_256_56 ROI_BRAINPOOL;	    /**< Root of unity set in rom_field_yyy.c */
extern const BIG_256_56 R2modp_BRAINPOOL;	/**< Montgomery constant */
extern const BIG_256_56 CRu_BRAINPOOL;       /**< Cube Root of Unity */
extern const BIG_256_56 SQRTm3_BRAINPOOL; /**< Square root of -3 */
extern const BIG_256_56 TWK_BRAINPOOL; /**< Tweak for square roots, pre-calculated from field norm */
extern const chunk MConst_BRAINPOOL;		/**< Constant associated with Modulus - for Montgomery = 1/p mod 2^BASEBITS */


#define MODBITS_BRAINPOOL MBITS_BRAINPOOL                        /**< Number of bits in Modulus for selected curve */
#define TBITS_BRAINPOOL (MBITS_BRAINPOOL%BASEBITS_256_56)           /**< Number of active bits in top word */
#define TMASK_BRAINPOOL (((chunk)1<<TBITS_BRAINPOOL)-1)          /**< Mask for active bits in top word */
#define FEXCESS_BRAINPOOL (((sign32)1<<MAXXES_BRAINPOOL)-1)	     /**< 2^(BASEBITS*NLEN-MODBITS)-1 - normalised BIG can be multiplied by less than this before reduction */
#define OMASK_BRAINPOOL (-((chunk)(1)<<TBITS_BRAINPOOL))         /**<  for masking out overflow bits */

//#define BIG_ENDIAN_SIGN_BRAINPOOL 

//#define FUSED_MODMUL
//#define DEBUG_REDUCE

/* FP prototypes */

/**	@brief Create FP from integer
 *
	@param x FP to be initialised
	@param a integer
 */
extern void FP_BRAINPOOL_from_int(FP_BRAINPOOL *x,int a);

/**	@brief Tests for FP equal to zero mod Modulus
 *
	@param x FP number to be tested
	@return 1 if zero, else returns 0
 */
extern int FP_BRAINPOOL_iszilch(FP_BRAINPOOL *x);


/**	@brief Tests for lexically largest 
 *
	@param x FP number to be tested if larger than -x
	@return 1 if larger, else returns 0
 */
extern int FP_BRAINPOOL_islarger(FP_BRAINPOOL *x);

/**	@brief Serialize out FP  
 *
    @param b buffer for output
	@param x FP number to be serialized
 */
extern void FP_BRAINPOOL_toBytes(char *b,FP_BRAINPOOL *x);

/**	@brief Serialize in FP  
 *
	@param x FP number to be serialized
    @param b buffer for input
 */
extern void FP_BRAINPOOL_fromBytes(FP_BRAINPOOL *x,char *b);



/**	@brief Tests for FP equal to one mod Modulus
 *
	@param x FP number to be tested
	@return 1 if one, else returns 0
 */
extern int FP_BRAINPOOL_isunity(FP_BRAINPOOL *x);


/**	@brief Set FP to zero
 *
	@param x FP number to be set to 0
 */
extern void FP_BRAINPOOL_zero(FP_BRAINPOOL *x);

/**	@brief Copy an FP
 *
	@param y FP number to be copied to
	@param x FP to be copied from
 */
extern void FP_BRAINPOOL_copy(FP_BRAINPOOL *y, FP_BRAINPOOL *x);

/**	@brief Copy from ROM to an FP
 *
	@param y FP number to be copied to
	@param x BIG to be copied from ROM
 */
extern void FP_BRAINPOOL_rcopy(FP_BRAINPOOL *y, const BIG_256_56 x);


/**	@brief Compares two FPs
 *
	@param x FP number
	@param y FP number
	@return 1 if equal, else returns 0
 */
extern int FP_BRAINPOOL_equals(FP_BRAINPOOL *x, FP_BRAINPOOL *y);


/**	@brief Conditional constant time swap of two FP numbers
 *
	Conditionally swaps parameters in constant time (without branching)
	@param x an FP number
	@param y another FP number
	@param s swap takes place if not equal to 0
 */
extern void FP_BRAINPOOL_cswap(FP_BRAINPOOL *x, FP_BRAINPOOL *y, int s);
/**	@brief Conditional copy of FP number
 *
	Conditionally copies second parameter to the first (without branching)
	@param x an FP number
	@param y another FP number
	@param s copy takes place if not equal to 0
 */
extern void FP_BRAINPOOL_cmove(FP_BRAINPOOL *x, FP_BRAINPOOL *y, int s);
/**	@brief Converts from BIG integer to residue form mod Modulus
 *
	@param x BIG number to be converted
	@param y FP result
 */
extern void FP_BRAINPOOL_nres(FP_BRAINPOOL *y, BIG_256_56 x);
/**	@brief Converts from residue form back to BIG integer form
 *
	@param y FP number to be converted to BIG
	@param x BIG result
 */
extern void FP_BRAINPOOL_redc(BIG_256_56 x, FP_BRAINPOOL *y);
/**	@brief Sets FP to representation of unity in residue form
 *
	@param x FP number to be set equal to unity.
 */
extern void FP_BRAINPOOL_one(FP_BRAINPOOL *x);


/**	@brief returns "sign" of an FP
 *
	@param x FP number
    @return 0 for positive, 1 for negative
 */
extern int FP_BRAINPOOL_sign(FP_BRAINPOOL *x);


/**	@brief Reduces DBIG to BIG exploiting special form of the modulus
 *
	This function comes in different flavours depending on the form of Modulus that is currently in use.
	@param r BIG number, on exit = d mod Modulus
	@param d DBIG number to be reduced
 */
extern void FP_BRAINPOOL_mod(BIG_256_56 r, DBIG_256_56 d);

#ifdef FUSED_MODMUL
extern void FP_BRAINPOOL_modmul(BIG_256_56, BIG_256_56, BIG_256_56);
#endif

/**	@brief Fast Modular multiplication of two FPs, mod Modulus
 *
	Uses appropriate fast modular reduction method
	@param x FP number, on exit the modular product = y*z mod Modulus
	@param y FP number, the multiplicand
	@param z FP number, the multiplier
 */
extern void FP_BRAINPOOL_mul(FP_BRAINPOOL *x, FP_BRAINPOOL *y, FP_BRAINPOOL *z);
/**	@brief Fast Modular multiplication of an FP, by a small integer, mod Modulus
 *
	@param x FP number, on exit the modular product = y*i mod Modulus
	@param y FP number, the multiplicand
	@param i a small number, the multiplier
 */
extern void FP_BRAINPOOL_imul(FP_BRAINPOOL *x, FP_BRAINPOOL *y, int i);
/**	@brief Fast Modular squaring of an FP, mod Modulus
 *
	Uses appropriate fast modular reduction method
	@param x FP number, on exit the modular product = y^2 mod Modulus
	@param y FP number, the number to be squared

 */
extern void FP_BRAINPOOL_sqr(FP_BRAINPOOL *x, FP_BRAINPOOL *y);
/**	@brief Modular addition of two FPs, mod Modulus
 *
	@param x FP number, on exit the modular sum = y+z mod Modulus
	@param y FP number
	@param z FP number
 */
extern void FP_BRAINPOOL_add(FP_BRAINPOOL *x, FP_BRAINPOOL *y, FP_BRAINPOOL *z);
/**	@brief Modular subtraction of two FPs, mod Modulus
 *
	@param x FP number, on exit the modular difference = y-z mod Modulus
	@param y FP number
	@param z FP number
 */
extern void FP_BRAINPOOL_sub(FP_BRAINPOOL *x, FP_BRAINPOOL *y, FP_BRAINPOOL *z);
/**	@brief Modular division by 2 of an FP, mod Modulus
 *
	@param x FP number, on exit =y/2 mod Modulus
	@param y FP number
 */
extern void FP_BRAINPOOL_div2(FP_BRAINPOOL *x, FP_BRAINPOOL *y);
/**	@brief Fast Modular exponentiation of an FP, to the power of a BIG, mod Modulus
 *
	@param x FP number, on exit  = y^z mod Modulus
	@param y FP number
	@param z BIG number exponent
 */
extern void FP_BRAINPOOL_pow(FP_BRAINPOOL *x, FP_BRAINPOOL *y, BIG_256_56 z);


/**	@brief Inverse square root precalculation
 *
	@param r FP number, on exit  = x^(p-2*e-1)/2^(e+1) mod Modulus
	@param x FP number
 */
extern void FP_BRAINPOOL_progen(FP_BRAINPOOL *r,FP_BRAINPOOL *x);

/**	@brief Fast Modular square root of a an FP, mod Modulus
 *
	@param x FP number, on exit  = sqrt(y) mod Modulus
	@param y FP number, the number whose square root is calculated
    @param h an optional precalculation
 */
extern void FP_BRAINPOOL_sqrt(FP_BRAINPOOL *x, FP_BRAINPOOL *y, FP_BRAINPOOL *h);

/**	@brief Modular negation of a an FP, mod Modulus
 *
	@param x FP number, on exit = -y mod Modulus
	@param y FP number
 */
extern void FP_BRAINPOOL_neg(FP_BRAINPOOL *x, FP_BRAINPOOL *y);
/**	@brief Outputs an FP number to the console
 *
	Converts from residue form before output
	@param x an FP number
 */
extern void FP_BRAINPOOL_output(FP_BRAINPOOL *x);
/**	@brief Outputs an FP number to the console, in raw form
 *
	@param x a BIG number
 */
extern void FP_BRAINPOOL_rawoutput(FP_BRAINPOOL *x);
/**	@brief Reduces possibly unreduced FP mod Modulus
 *
	@param x FP number, on exit reduced mod Modulus
 */
extern void FP_BRAINPOOL_reduce(FP_BRAINPOOL *x);
/**	@brief normalizes FP
 *
	@param x FP number, on exit normalized
 */
extern void FP_BRAINPOOL_norm(FP_BRAINPOOL *x);
/**	@brief Tests for FP a quadratic residue mod Modulus
 *
	@param x FP number to be tested
    @param h an optional precalculation
	@return 1 if quadratic residue, else returns 0 if quadratic non-residue
 */
extern int FP_BRAINPOOL_qr(FP_BRAINPOOL *x,FP_BRAINPOOL *h);

/**	@brief Simultaneous Inverse and Square root
 *
	@param i FP number, on exit = 1/x mod Modulus
	@param s FP number, on exit = sqrt(x) mod Modulus
	@param x FP number
	@return 1 if quadratic residue, else returns 0 if quadratic non-residue
 */
extern int FP_BRAINPOOL_invsqrt(FP_BRAINPOOL *i,FP_BRAINPOOL *s,FP_BRAINPOOL *x);

/**	@brief Simultaneous Inverse and Square root of different numbers
 *
	@param i FP number, on exit = 1/i mod Modulus
	@param s FP number, on exit = sqrt(s) mod Modulus
	@return 1 if quadratic residue, else returns 0 if quadratic non-residue
 */
extern int FP_BRAINPOOL_tpo(FP_BRAINPOOL *i, FP_BRAINPOOL *s);



/**	@brief Modular inverse of a an FP, mod Modulus
 *
	@param x FP number, on exit = 1/y mod Modulus
	@param y FP number
    @param h an optional input precalculation
 */
extern void FP_BRAINPOOL_inv(FP_BRAINPOOL *x, FP_BRAINPOOL *y, FP_BRAINPOOL *h);

/**	@brief Generate random FP
 *
	@param x random FP number
	@param rng random number generator
 */
extern void FP_BRAINPOOL_rand(FP_BRAINPOOL *x, csprng *rng);



#endif