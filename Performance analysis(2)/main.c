#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "pair_BLS12381.h"
#include "test/cpucycles.h"

#define ATT_SIZE 15
#define MLEN 59

typedef struct {
	ECP_BLS12381 cmt;
	BIG_384_58 rsp;
}prf;

typedef struct {
	ECP2_BLS12381 cmt;
	BIG_384_58 rsp[ATT_SIZE + 2];
}prf2;

static int cmp_llu(const void *a, const void *b)
{
	if (*(unsigned long long *)a < *(unsigned long long *)b) return -1;
	if (*(unsigned long long *)a > *(unsigned long long *)b) return 1;
	return 0;
}

static unsigned long long median(unsigned long long *l, size_t llen)
{
	qsort(l, llen, sizeof(unsigned long long), cmp_llu);

	if (llen % 2) return l[llen / 2];
	else return (l[llen / 2 - 1] + l[llen / 2]) / 2;
}

static unsigned long long average(unsigned long long *t, size_t tlen)
{
	unsigned long long acc = 0;
	size_t i;
	for (i = 0; i < tlen; i++)
		acc += t[i];
	return acc / (tlen);
}

static void print_results(const char *s, unsigned long long *t, size_t tlen)
{
	printf("%s", s);
	printf("\n");
	printf("median:  %llu ", median(t, tlen));  print_unit; printf("\n");
	printf("average: %llu ", average(t, tlen - 1));  print_unit; printf("\n");
	printf("\n");
}

void HASH_CORE(octet *out, octet *in, unsigned char *domain);
void HASH_H1(BIG_384_58 out, ECP_BLS12381 *cmt, int len);
void HASH_H2(BIG_384_58 out, ECP2_BLS12381 *cmt, int len);

void PoK(ECP_BLS12381 *pk1, BIG_384_58 usk, ECP_BLS12381 upk, prf *p);
int Vrf_PoK(ECP_BLS12381 *pk1, ECP_BLS12381 upk, prf p);

void NIZK(ECP2_BLS12381 *pk2, BIG_384_58 usk, BIG_384_58 t, int *I, char m[][MLEN], prf2 *p);
int Vrf_NIZK(ECP2_BLS12381 *pk2, ECP2_BLS12381 usig3, char rm[][MLEN], int I[], prf2 p);

void OrgKeyGen(BIG_384_58 *sk, FP12_BLS12381 *pk, BIG_384_58 *sk1, ECP_BLS12381 *pk1, ECP2_BLS12381 *pk2);

void UserKeyGen(ECP_BLS12381 *pk1, BIG_384_58 *usk, ECP_BLS12381 *upk);

void Obt_and_Iss(BIG_384_58 usk, ECP_BLS12381 upk, BIG_384_58 sk, ECP_BLS12381 *pk1, prf *p, char m[][MLEN], ECP_BLS12381 *sig1, ECP_BLS12381 *sig2);
void Obtain(ECP_BLS12381 *pk1, BIG_384_58 usk, ECP_BLS12381 upk, prf *p, char m[][MLEN]);
void Issue(ECP_BLS12381 upk, BIG_384_58 sk, ECP_BLS12381 *pk1, prf p, char m[][MLEN], ECP_BLS12381 *sig1, ECP_BLS12381 *sig2);

void Show_and_Vrf(ECP_BLS12381 *pk1, ECP2_BLS12381 *pk2, FP12_BLS12381 pk, BIG_384_58 usk, ECP_BLS12381 upk, char m[][MLEN], int *I, prf2 *p, ECP_BLS12381 sig1, ECP_BLS12381 sig2, ECP_BLS12381 *usig1, ECP_BLS12381 *usig2, ECP2_BLS12381 *usig3);
void Show(ECP_BLS12381 *pk1, ECP2_BLS12381 *pk2, BIG_384_58 usk, ECP_BLS12381 upk, char m[][MLEN], int *I, prf2 *p, ECP_BLS12381 sig1, ECP_BLS12381 sig2, ECP_BLS12381 *usig1, ECP_BLS12381 *usig2, ECP2_BLS12381 *usig3);
void Verify(FP12_BLS12381 pk, ECP2_BLS12381 *pk2, char rm[][MLEN], int I[], prf2 p, ECP_BLS12381 usig1, ECP_BLS12381 usig2, ECP2_BLS12381 usig3);

int main(void)
{
	srand(time(NULL));
	unsigned int i, j;

	//performance
	struct timespec begin, end;
	double time = 0.0;

	//Message and initialization
	char str[10] = "ATTRIBUTE";
	char m[ATT_SIZE][MLEN];
	for (j = 0; j < ATT_SIZE; j++)
	{
		strcpy(m[j], str);
	}

	//index of message to be opened (random selection)
	int I[ATT_SIZE];
	for (j = 0; j < ATT_SIZE; j++)
	{
		I[j] = rand() % 2;
	}

	//proof
	prf p1;
	prf2 p2;

	//Organization Key
	BIG_384_58 sk;						//sk
	FP12_BLS12381 pk;					//e(g, g')^sk
	BIG_384_58 sk1[ATT_SIZE + 2];		//y0, y1, ..., yn, h
	ECP_BLS12381 pk1[ATT_SIZE + 3];		//g0, g1, ..., gn, h, g
	ECP2_BLS12381 pk2[ATT_SIZE + 3];	//g0', g1', ..., gn', h', g'

	//User Key
	BIG_384_58 usk;
	ECP_BLS12381 upk;

	//Signature from Org
	ECP_BLS12381 sig1, sig2;

	//Signature from User
	ECP_BLS12381 usig1, usig2;
	ECP2_BLS12381 usig3;
	
	
	OrgKeyGen(&sk, &pk, sk1, pk1, pk2);
	
	clock_gettime(CLOCK_MONOTONIC, &begin); //performance
	UserKeyGen(pk1, &usk, &upk);
	clock_gettime(CLOCK_MONOTONIC, &end); //performance

	Obt_and_Iss(usk, upk, sk, pk1, &p1, m, &sig1, &sig2);

	Show_and_Vrf(pk1, pk2, pk, usk, upk, m, I, &p2, sig1, sig2, &usig1, &usig2, &usig3);

	printf("\n");
	printf("...VALID CREDENTIALS...\n\n");

	time = (end.tv_sec - begin.tv_sec) + 1e-9 * (end.tv_nsec - begin.tv_nsec);
	time = time * 1000.0;

	printf("Time cost per one voter (milli unit): %.6lf ms\n", time);

	return 0;
}

void HASH_CORE(octet *out, octet *in, unsigned char *domain)
{
	unsigned char *m = malloc(in->len + 8);
	octet M = { 8, in->len + 8, m };

	memcpy(m, domain, 8);
	OCT_joctet(&M, in);

	SPhash(MC_SHA2, 48, out, &M);

	free(m);
}

void HASH_H1(BIG_384_58 out, ECP_BLS12381 *cmt, int len)
{
	BIG_384_58 p;

	unsigned char domain[8] = { 'H','A','S','H','2', 0, 0, 0 };
	unsigned char element[49] = { 0 };

	unsigned char *message;
	unsigned char digest[48] = { 0 };

	octet Element = { 0, 49, element };
	octet Message;
	octet Digest = { 0, 48, digest };

	//Initialize Message
	message = (unsigned char *)malloc(4 + 49 * len);
	Message.len = 4 + 49 * len;
	Message.val = message;
	Message.max = 4 + 49 * len;

	//Fill Message
	memcpy(message, &len, sizeof(len));
	for (int i = 0; i < len; i++)
	{
		ECP_BLS12381_toOctet(&Element, cmt + i, 1);
		memcpy(message + 4 + i * 49, element, 49);
	}

	HASH_CORE(&Digest, &Message, domain);
	BIG_384_58_fromBytesLen(out, Digest.val, Digest.len);

	BIG_384_58_rcopy(p, CURVE_Order_BLS12381);

	BIG_384_58_mod(out, p);

	free(message);
}

void HASH_H2(BIG_384_58 out, ECP2_BLS12381 *cmt, int len)
{
	BIG_384_58 p;

	unsigned char domain[8] = { 'H','A','S','H','2', 0, 0, 0 };
	unsigned char element[97] = { 0 };

	unsigned char *message;
	unsigned char digest[48] = { 0 };

	octet Element = { 0, 97, element };
	octet Message;
	octet Digest = { 0, 48, digest };

	//Initialize Message
	message = (unsigned char *)malloc(4 + 97 * len);
	Message.len = 4 + 97 * len;
	Message.val = message;
	Message.max = 4 + 97 * len;

	//Fill Message
	memcpy(message, &len, sizeof(len));
	for (int i = 0; i < len; i++)
	{
		ECP2_BLS12381_toOctet(&Element, cmt + i, 1);
		memcpy(message + 4 + i * 97, element, 97);
	}

	HASH_CORE(&Digest, &Message, domain);
	BIG_384_58_fromBytesLen(out, Digest.val, Digest.len);

	BIG_384_58_rcopy(p, CURVE_Order_BLS12381);

	BIG_384_58_mod(out, p);

	free(message);
}

void PoK(ECP_BLS12381 *pk1, BIG_384_58 usk, ECP_BLS12381 upk, prf *p)
{
	unsigned char seed[32] = { 0 };
	csprng CSPRNG;
	BIG_384_58 order;

	RAND_seed(&CSPRNG, 32, seed);
	BIG_384_58_rcopy(order, CURVE_Order_BLS12381);

	BIG_384_58 out;
	ECP_BLS12381 in[3];

	BIG_384_58 r;
	ECP_BLS12381 cmt_r;

	BIG_384_58_randomnum(r, order, &CSPRNG); //random r 추출

	ECP_BLS12381_copy(&cmt_r, &pk1[ATT_SIZE + 1]); //h
	ECP_BLS12381_mul(&cmt_r, r); //r*h

	ECP_BLS12381_generator(&in[0]);
	ECP_BLS12381_copy(&in[1], &cmt_r);
	ECP_BLS12381_copy(&in[2], &upk);

	//HASH
	HASH_H1(out, in, 3);;

	//compute response
	BIG_384_58_modmul(p->rsp, usk, out, order);		//p->rsp = x*c
	BIG_384_58_modneg(p->rsp, p->rsp, order);		//p->rsp = -x*c
	BIG_384_58_modadd(p->rsp, r, p->rsp, order);	//p->rsp = r-x*c

	ECP_BLS12381_copy(&p->cmt, &cmt_r);				//p->cmt = r*h
}

int Vrf_PoK(ECP_BLS12381 *pk1, ECP_BLS12381 upk, prf p)
{
	BIG_384_58 out;
	ECP_BLS12381 in[3];
	ECP_BLS12381 tmp;

	ECP_BLS12381_copy(&tmp, &pk1[ATT_SIZE + 1]);	//tmp = h

	ECP_BLS12381_generator(&in[0]);
	ECP_BLS12381_copy(&in[1], &p.cmt);
	ECP_BLS12381_copy(&in[2], &upk);

	//HASH
	HASH_H1(out, in, 3);

	//compare
	ECP_BLS12381_mul2(&tmp, &upk, p.rsp, out);

	if (ECP_BLS12381_equals(&tmp, &p.cmt))
		return 1;
	else
		return 0;
}

void NIZK(ECP2_BLS12381 *pk2, BIG_384_58 usk, BIG_384_58 t, int *I, char m[][MLEN], prf2 *p)
{
	int i;

	unsigned char seed[32] = { 0 };
	csprng CSPRNG;
	BIG_384_58 order;

	RAND_seed(&CSPRNG, 32, seed);
	BIG_384_58_rcopy(order, CURVE_Order_BLS12381);

	BIG_384_58 bm[ATT_SIZE];

	BIG_384_58 out;
	ECP2_BLS12381 tmp;
	ECP2_BLS12381 in[ATT_SIZE + 4];	//g0', g1', ..., gn', h', g', cmt

	BIG_384_58 r[ATT_SIZE];
	BIG_384_58 usk_r, t_r;
	ECP2_BLS12381 cmt;

	//random r 추출
	BIG_384_58_randomnum(usk_r, order, &CSPRNG);
	BIG_384_58_randomnum(t_r, order, &CSPRNG);

	for (i = 0; i < ATT_SIZE; i++)
	{
		if (I[i] == 0)
			BIG_384_58_randomnum(r[i], order, &CSPRNG);
	}

	//commitment
	ECP2_BLS12381_copy(&cmt, &pk2[0]);		//cmt = g0'

	for (i = 0; i < ATT_SIZE; i++)
	{
		BIG_384_58_fromBytesLen(bm[i], m[i], MLEN);
	}

	for (i = 0; i < ATT_SIZE; i++)
	{
		if (I[i])
		{
			ECP2_BLS12381_copy(&tmp, &pk2[i + 1]);
			ECP2_BLS12381_mul(&tmp, bm[i]);				//tmp = m[i]*gi'

			ECP2_BLS12381_add(&cmt, &tmp);				//cmt = cmt + m[i]&gi'
		}

		else
		{
			ECP2_BLS12381_copy(&tmp, &pk2[i + 1]);
			ECP2_BLS12381_mul(&tmp, r[i]);				//tmp = r[i]*gi'

			ECP2_BLS12381_add(&cmt, &tmp);				//cmt = cmt + r[i]*gi'
		}
	}

	ECP2_BLS12381_copy(&tmp, &pk2[ATT_SIZE + 1]);
	ECP2_BLS12381_mul(&tmp, usk_r);
	ECP2_BLS12381_add(&cmt, &tmp);				//cmt = cmt + usk*h'

	ECP2_BLS12381_copy(&tmp, &pk2[ATT_SIZE + 2]);
	ECP2_BLS12381_mul(&tmp, t_r);
	ECP2_BLS12381_add(&cmt, &tmp);				//cmt = cmt + t*g'

	//copy hash input
	for (i = 0; i < ATT_SIZE + 3; i++)
	{
		ECP2_BLS12381_copy(&in[i], &pk2[i]);
	}

	ECP2_BLS12381_copy(&in[ATT_SIZE + 3], &cmt);

	//HASH
	HASH_H2(out, in, ATT_SIZE + 4);

	//compute response		p->[i] for all i\in[ATT_SIZE + 2]에 대해 [0, ATT_SIZE-1]에 r1~rn을, p->[ATT_SIZE], p->[ATT_SIZE+1]에는 각각 usk, t에 대한 response 저장
	for (i = 0; i < ATT_SIZE; i++)
	{
		if (I[i] == 0)
		{
			BIG_384_58_modmul(p->rsp[i], bm[i], out, order);		//p->rsp[i] = m[i]*c
			BIG_384_58_modneg(p->rsp[i], p->rsp[i], order);			//p->rsp[i] = -m[i]*c
			BIG_384_58_modadd(p->rsp[i], r[i], p->rsp[i], order);	//p->rsp[i] = r[i]-m[i]*c
		}
	}

	BIG_384_58_modmul(p->rsp[ATT_SIZE], usk, out, order);					//p->rsp[i] = usk*c
	BIG_384_58_modneg(p->rsp[ATT_SIZE], p->rsp[ATT_SIZE], order);			//p->rsp[i] = -usk*c
	BIG_384_58_modadd(p->rsp[ATT_SIZE], usk_r, p->rsp[ATT_SIZE], order);	//p->rsp[i] = usk_r - usk*c

	BIG_384_58_modmul(p->rsp[ATT_SIZE + 1], t, out, order);						//p->rsp[i] = t*c
	BIG_384_58_modneg(p->rsp[ATT_SIZE + 1], p->rsp[ATT_SIZE + 1], order);		//p->rsp[i] = -t*c
	BIG_384_58_modadd(p->rsp[ATT_SIZE + 1], t_r, p->rsp[ATT_SIZE + 1], order);	//p->rsp[i] = t_r - t*c

	ECP2_BLS12381_copy(&p->cmt, &cmt);		//p->cmt = cmt
}

int Vrf_NIZK(ECP2_BLS12381 *pk2, ECP2_BLS12381 usig3, char rm[][MLEN], int I[], prf2 p)
{
	int i;
	BIG_384_58 order;
	BIG_384_58 bm[ATT_SIZE];
	BIG_384_58 out;
	ECP2_BLS12381 in[ATT_SIZE + 4];	//g0', g1', ..., gn', h', g', cmt
	ECP2_BLS12381 tmp1, tmp2, tmp3;

	BIG_384_58_rcopy(order, CURVE_Order_BLS12381);

	//copy hash input
	for (i = 0; i < ATT_SIZE + 3; i++)
	{
		ECP2_BLS12381_copy(&in[i], &pk2[i]);
	}

	ECP2_BLS12381_copy(&in[ATT_SIZE + 3], &p.cmt);

	//HASH
	HASH_H2(out, in, ATT_SIZE + 4);

	//compare
	ECP2_BLS12381_copy(&tmp1, &usig3);

	//usig3에서 공개된 m[i]는 역원을 곱하여 소거->그 결과를 tmp1에 저장
	ECP2_BLS12381_copy(&tmp2, &pk2[0]);
	ECP2_BLS12381_sub(&tmp1, &tmp2);

	for (i = 0; i < ATT_SIZE; i++)
	{
		if (I[i] == 1)
		{
			BIG_384_58_fromBytesLen(bm[i], rm[i], MLEN);

			ECP2_BLS12381_copy(&tmp2, &pk2[i + 1]);
			ECP2_BLS12381_mul(&tmp2, bm[i]);
			ECP2_BLS12381_sub(&tmp1, &tmp2);
		}
	}

	//tmp3에 stmt에 해당하는 base 각각에 대응하는 response를 곱한 결과 저장
	ECP2_BLS12381_copy(&tmp3, &pk2[0]);		//cmt = g0'

	for (i = 0; i < ATT_SIZE; i++)
	{
		if (I[i])
		{
			ECP2_BLS12381_copy(&tmp2, &pk2[i + 1]);
			ECP2_BLS12381_mul(&tmp2, bm[i]);				//tmp2 = m[i]*gi'
			ECP2_BLS12381_add(&tmp3, &tmp2);				//tmp3 = tmp3 + m[i]&gi'
		}

		else
		{
			ECP2_BLS12381_copy(&tmp2, &pk2[i + 1]);
			ECP2_BLS12381_mul(&tmp2, p.rsp[i]);				//tmp2 = p.rsp[i]*gi'
			ECP2_BLS12381_add(&tmp3, &tmp2);				//tmp3 = tmp3 + p.rsp[i]*gi'
		}
	}

	ECP2_BLS12381_copy(&tmp2, &pk2[ATT_SIZE + 1]);
	ECP2_BLS12381_mul(&tmp2, p.rsp[ATT_SIZE]);
	ECP2_BLS12381_add(&tmp3, &tmp2);				//tmp3 = tmp3 + usk*h'

	ECP2_BLS12381_copy(&tmp2, &pk2[ATT_SIZE + 2]);
	ECP2_BLS12381_mul(&tmp2, p.rsp[ATT_SIZE + 1]);
	ECP2_BLS12381_add(&tmp3, &tmp2);				//tmp3 = tmp3 + t*g'

	ECP2_BLS12381_mul(&tmp1, out);		//tmp1 = c*tmp1
	ECP2_BLS12381_add(&tmp1, &tmp3);	//tmp1 = tmp3 + tmp1

	//compare
	if (ECP2_BLS12381_equals(&tmp1, &p.cmt))
		return 1;
	else
		return 0;
}

void OrgKeyGen(BIG_384_58 *sk, FP12_BLS12381 *pk, BIG_384_58 *sk1, ECP_BLS12381 *pk1, ECP2_BLS12381 *pk2)
{
	unsigned char seed[32] = { 0 };
	csprng CSPRNG;
	BIG_384_58 order;
	int i;

	RAND_seed(&CSPRNG, 32, seed);
	BIG_384_58_rcopy(order, CURVE_Order_BLS12381);

	BIG_384_58_randomnum(*sk, order, &CSPRNG);

	ECP_BLS12381_generator(&pk1[ATT_SIZE + 2]);		//g
	ECP2_BLS12381_generator(&pk2[ATT_SIZE + 2]);	//g'

	for (i = 0; i < ATT_SIZE + 2; i++)
	{
		ECP_BLS12381_copy(&pk1[i], &pk1[ATT_SIZE + 2]);
		ECP2_BLS12381_copy(&pk2[i], &pk2[ATT_SIZE + 2]);

		BIG_384_58_randomnum(sk1[i], order, &CSPRNG);

		ECP_BLS12381_mul(&pk1[i], sk1[i]);
		ECP2_BLS12381_mul(&pk2[i], sk1[i]);
	}

	PAIR_BLS12381_ate(pk, &pk2[ATT_SIZE + 2], &pk1[ATT_SIZE + 2]);		//pk = e(g, g')
	PAIR_BLS12381_fexp(pk);
	PAIR_BLS12381_GTpow(pk, *sk);										//pk = e(g, g')^delta
}

void UserKeyGen(ECP_BLS12381 *pk1, BIG_384_58 *usk, ECP_BLS12381 *upk)
{
	unsigned char seed[32] = { 0 };
	csprng CSPRNG;
	BIG_384_58 order;

	RAND_seed(&CSPRNG, 32, seed);
	BIG_384_58_rcopy(order, CURVE_Order_BLS12381);

	BIG_384_58_randomnum(*usk, order, &CSPRNG);

	ECP_BLS12381_copy(upk, &pk1[ATT_SIZE + 1]);	//upk = h
	ECP_BLS12381_mul(upk, *usk);				//upk = h^usk
}

void Obt_and_Iss(BIG_384_58 usk, ECP_BLS12381 upk, BIG_384_58 sk, ECP_BLS12381 *pk1, prf *p, char m[][MLEN], ECP_BLS12381 *sig1, ECP_BLS12381 *sig2)
{


	
	Obtain(pk1, usk, upk, p, m);
	
	
	
	Issue(upk, sk, pk1, *p, m, sig1, sig2);
	

	
}

void Obtain(ECP_BLS12381 *pk1, BIG_384_58 usk, ECP_BLS12381 upk, prf *p, char m[][MLEN])
{
	PoK(pk1, usk, upk, p);
}

void Issue(ECP_BLS12381 upk, BIG_384_58 sk, ECP_BLS12381 *pk1, prf p, char m[][MLEN], ECP_BLS12381 *sig1, ECP_BLS12381 *sig2)
{
	int i;

	if (Vrf_PoK(pk1, upk, p) == 0)
	{
		printf("INVALID USER\n");
		exit(0);
	}

	unsigned char seed[32] = { 0 };
	csprng CSPRNG;
	BIG_384_58 order;

	RAND_seed(&CSPRNG, 32, seed);
	BIG_384_58_rcopy(order, CURVE_Order_BLS12381);

	BIG_384_58 r;
	BIG_384_58 bm[ATT_SIZE];
	ECP_BLS12381 ym[ATT_SIZE];
	ECP_BLS12381 tmp;

	//convert to BIG number from byte array
	for (i = 0; i < ATT_SIZE; i++)
	{
		BIG_384_58_fromBytesLen(bm[i], m[i], MLEN);
	}

	BIG_384_58_randomnum(r, order, &CSPRNG); //random r 추출

	ECP_BLS12381_generator(sig1);	//sig1 = G
	ECP_BLS12381_mul(sig1, r);		//sig1 = r*G

	for (i = 0; i < ATT_SIZE; i++)
	{
		ECP_BLS12381_copy(&ym[i], &pk1[i + 1]);
		ECP_BLS12381_mul(&ym[i], bm[i]);
	}

	ECP_BLS12381_copy(&tmp, &upk);		//tmp = (h*usk)*G
	ECP_BLS12381_add(&tmp, &pk1[0]);	//tmp = (h*usk + y0)*G

	for (i = 0; i < ATT_SIZE; i++)
	{
		ECP_BLS12381_add(&tmp, &ym[i]);	//tmp = (h*usk + y0 + y1*m1 + ... + yn*mn)*G
	}

	ECP_BLS12381_mul(&tmp, r);			//tmp = r*(h*usk + y0 + y1*m1 + ... + yn*mn)*G

	ECP_BLS12381_generator(sig2);		//sig2 = G
	ECP_BLS12381_mul(sig2, sk);			//sig2 = sk*G

	ECP_BLS12381_add(sig2, &tmp);		//sig2 = (r*(h*usk + y0 + y1*m1 + ... + yn*mn) + sk)*G
}

void Show_and_Vrf(ECP_BLS12381 *pk1, ECP2_BLS12381 *pk2, FP12_BLS12381 pk, BIG_384_58 usk, ECP_BLS12381 upk, char m[][MLEN], int *I, prf2 *p, ECP_BLS12381 sig1, ECP_BLS12381 sig2, ECP_BLS12381 *usig1, ECP_BLS12381 *usig2, ECP2_BLS12381 *usig3)
{	
	Show(pk1, pk2, usk, upk, m, I, p, sig1, sig2, usig1, usig2, usig3);

	Verify(pk, pk2, m, I, *p, *usig1, *usig2, *usig3);
}

void Show(ECP_BLS12381 *pk1, ECP2_BLS12381 *pk2, BIG_384_58 usk, ECP_BLS12381 upk, char m[][MLEN], int *I, prf2 *p, ECP_BLS12381 sig1, ECP_BLS12381 sig2, ECP_BLS12381 *usig1, ECP_BLS12381 *usig2, ECP2_BLS12381 *usig3)
{
	int i;

	unsigned char seed[32] = { 0 };
	csprng CSPRNG;
	BIG_384_58 order;
	char str[9] = "REDACTED";

	RAND_seed(&CSPRNG, 32, seed);
	BIG_384_58_rcopy(order, CURVE_Order_BLS12381);

	BIG_384_58 r, t;
	BIG_384_58 bm[ATT_SIZE];
	ECP_BLS12381 ym[ATT_SIZE];
	ECP2_BLS12381 ym2[ATT_SIZE];
	ECP_BLS12381 tmp1;
	ECP2_BLS12381 tmp2;

	BIG_384_58_randomnum(r, order, &CSPRNG); //random r' 추출 (편의상 r'이 아닌 r로 표기)
	BIG_384_58_randomnum(t, order, &CSPRNG); //random t 추출

	//usig1
	ECP_BLS12381_generator(usig1);		//usig1 = G
	ECP_BLS12381_mul(usig1, r);			//usig1 = r'*G
	ECP_BLS12381_add(usig1, &sig1);		//usig1 = (r+r')*G

	//usig2
	//convert to BIG number from byte array
	for (i = 0; i < ATT_SIZE; i++)
	{
		BIG_384_58_fromBytesLen(bm[i], m[i], MLEN);
	}

	ECP_BLS12381_copy(&tmp1, &upk);		//tmp1 = (h*usk)*G
	ECP_BLS12381_add(&tmp1, &pk1[0]);	//tmp1 = (h*usk + y0)*G

	for (i = 0; i < ATT_SIZE; i++)
	{
		ECP_BLS12381_copy(&ym[i], &pk1[i + 1]);
		ECP_BLS12381_mul(&ym[i], bm[i]);
		ECP_BLS12381_add(&tmp1, &ym[i]);	//tmp1 = (h*usk + y0 + y1*m1 + ... + yn*mn)*G
	}

	ECP_BLS12381_mul(&tmp1, r);			//tmp1 = r*(h*usk + y0 + y1*m1 + ... + yn*mn)*G

	ECP_BLS12381_add(&tmp1, &sig2);		//tmp1 = (sk + (h*usk + y0 + y1*m1 + ... + yn*mn)*(r+r'))*G

	ECP_BLS12381_copy(usig2, usig1);	//usig2 = usig1 = (r+r')*G
	ECP_BLS12381_mul(usig2, t);			//usig2 = t*(r+r')*G

	ECP_BLS12381_add(usig2, &tmp1);		//usig2 = (sk + (h*usk + y0 + y1*m1 + ... + yn*mn + t)*(r+r'))*G

	//usig3
	ECP2_BLS12381_copy(usig3, &pk2[ATT_SIZE + 1]);		//usig3 = h'*G'
	ECP2_BLS12381_mul(usig3, usk);						//usig3 = (h'*usk)*G'
	ECP2_BLS12381_add(usig3, &pk2[0]);					//usig3 = (h'*usk + y0')*G'

	for (i = 0; i < ATT_SIZE; i++)
	{
		ECP2_BLS12381_copy(&ym2[i], &pk2[i + 1]);
		ECP2_BLS12381_mul(&ym2[i], bm[i]);
		ECP2_BLS12381_add(usig3, &ym2[i]);				//usig3 = (h'*usk + y0' + y1'*m1 + ... + yn*mn)*G'
	}

	ECP2_BLS12381_generator(&tmp2);			//tmp2 = G'
	ECP2_BLS12381_mul(&tmp2, t);			//tmp2 = t*G'
	ECP2_BLS12381_add(usig3, &tmp2);		//usig3 = (h'*usk + y0' + y1'*m1 + ... + yn*mn + t)*G'

	//NIZK proof
	NIZK(pk2, usk, t, I, m, p);

	//redact messages
	for (i = 0; i < ATT_SIZE; i++)
	{
		if (I[i] == 0)
			strcpy(m[i], str);
	}
}

void Verify(FP12_BLS12381 pk, ECP2_BLS12381 *pk2, char rm[][MLEN], int I[], prf2 p, ECP_BLS12381 usig1, ECP_BLS12381 usig2, ECP2_BLS12381 usig3)
{
	ECP2_BLS12381 G2;
	FP12_BLS12381 tmp1, tmp2;

	ECP2_BLS12381_generator(&G2);

	//check 1
	PAIR_BLS12381_ate(&tmp1, &G2, &usig2);
	PAIR_BLS12381_fexp(&tmp1);
	
	PAIR_BLS12381_ate(&tmp2, &usig3, &usig1);
	PAIR_BLS12381_fexp(&tmp2);
	FP12_BLS12381_mul(&tmp2, &pk);
	
	if (FP12_BLS12381_equals(&tmp1, &tmp2) == 0)
	{
		printf("(Type1) INVALID CREDENTIALS\n");
		exit(0);
	}

	//check 2
	if (Vrf_NIZK(pk2, usig3, rm, I, p) == 0)
	{
		printf("(Type2) INVALID CREDENTIALS\n");
		exit(0);
	}
}