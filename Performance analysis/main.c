#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ecp_NIST256.h"
#include "big_256_56.h"

typedef struct {
	ECP_NIST256 E[8];
	BIG_256_56 B[8];
}rsp;

void HASH_CORE(octet *out, octet *in, unsigned char *domain);
void HASH_H1(BIG_256_56 out, ECP_NIST256 *cmt, int len);

BIG_256_56 **SETUP_SK();
ECP_NIST256 **SETUP_PK();
ECP_NIST256 **SETUP_Z();
ECP_NIST256 **SETUP_Y();
rsp **SETUP_PROOF();

void ROUND1(BIG_256_56 **sk, ECP_NIST256 **pk, rsp **pi1);
void GENERATE_SK(BIG_256_56 **sk);
void GENERATE_PK(BIG_256_56 **sk, ECP_NIST256 **pk);
void NIZK1(ECP_NIST256 **pk, BIG_256_56 **sk, rsp **pi1);

void COMPUTE_Z(ECP_NIST256 **Z, ECP_NIST256 **pk);
void COMPUTE_Y(ECP_NIST256 **Y, ECP_NIST256 **pk);

void ROUND2(ECP_NIST256 **H_Y, ECP_NIST256 **H_Z, ECP_NIST256 **Y, ECP_NIST256 **Z, ECP_NIST256 **pk, BIG_256_56 **sk, int *rnd, rsp **pi2, rsp **pi3);
void COMPUTE_HZ(ECP_NIST256 **H_Z, ECP_NIST256 **Z, BIG_256_56 **sk, int *rnd);
void NIZK2(int *rnd, ECP_NIST256 **H_Z, ECP_NIST256 **Z, ECP_NIST256 **pk, BIG_256_56 **sk, rsp **pi2);
void COMPUTE_HY(ECP_NIST256 **H_Y, ECP_NIST256 **Y, BIG_256_56 **sk, int *rnd);
void NIZK3(int *rnd, ECP_NIST256 **H_Z, ECP_NIST256 **H_Y, ECP_NIST256 **Z, ECP_NIST256 **Y, BIG_256_56 **sk, rsp **pi3);

void SELF_TALLY(ECP_NIST256 **H_Y);

void CLEAR(ECP_NIST256 **E1, ECP_NIST256 **E2, ECP_NIST256 **E3, BIG_256_56 **B4, rsp **r1, rsp **r2, rsp **r3);

int VERIFY_NIZK1(ECP_NIST256 **pk, rsp **pi1);
int VERIFY_NIZK2(ECP_NIST256 **pk, ECP_NIST256 **Z, ECP_NIST256 **H_Z, rsp **pi2);
int VERIFY_NIZK3(ECP_NIST256 **Z, ECP_NIST256 **Y, ECP_NIST256 **H_Z, ECP_NIST256 **H_Y, rsp **pi3);

unsigned int num_vt;		//The number of voters
unsigned int num_cnd;		//The number of candidates

int main()
{
	srand(time(NULL));
	int i, j, k;

	num_vt = 40;

	for (k = 1; k <= 1; k++)
	{
		num_cnd = 2;

		int *rnd = (int *)malloc(sizeof(int) * num_vt);
		int rslt[3] = { 0 };

		// performance
		struct timespec begin, end;
		double time_vote = 0.0;
		double time_tally = 0.0;
		double time = 0.0;
		//



		ECP_NIST256 **H_Y = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);
		ECP_NIST256 **H_Z = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);

		for (i = 0; i < num_vt; i++)
		{
			H_Y[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
			H_Z[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);

			//rnd[i] = rand() % num_cnd;
		}

		for (i = 0; i < num_vt; i++)
		{
			j = i % num_cnd;
			rnd[i] = j;
		}

		BIG_256_56 **sk = SETUP_SK();
		ECP_NIST256 **pk = SETUP_PK();
		ECP_NIST256 **Z = SETUP_Z();
		ECP_NIST256 **Y = SETUP_Y();

		rsp **pi1 = SETUP_PROOF();
		rsp **pi2 = SETUP_PROOF();
		rsp **pi3 = SETUP_PROOF();

		clock_gettime(CLOCK_MONOTONIC, &begin); //performance
		ROUND1(sk, pk, pi1);

		COMPUTE_Z(Z, pk);
		COMPUTE_Y(Y, pk);

		ROUND2(H_Y, H_Z, Y, Z, pk, sk, rnd, pi2, pi3);

		//rslt[0] = VERIFY_NIZK1(pk, pi1);
		////printf("***%d\n", rslt[0]);

		//rslt[1] = VERIFY_NIZK2(pk, Z, H_Z, pi2);
		////printf("***%d\n", rslt[1]);

		//rslt[2] = VERIFY_NIZK3(Z, Y, H_Z, H_Y, pi3);
		////printf("***%d\n", rslt[2]);

		clock_gettime(CLOCK_MONOTONIC, &end); //performance

		time_vote = (end.tv_sec - begin.tv_sec) + 1e-9 * (end.tv_nsec - begin.tv_nsec);
		time_vote = time_vote / (double)num_vt;
		time_vote = time_vote * 1000.0;

		clock_gettime(CLOCK_MONOTONIC, &begin); //performance

		SELF_TALLY(H_Y);

		clock_gettime(CLOCK_MONOTONIC, &end); //performance

		//free
		CLEAR(Z, Y, pk, sk, pi1, pi2, pi3);

		for (i = 0; i < num_vt; i++)
		{
			free(H_Y[i]);
			free(H_Z[i]);
		}

		free(rnd);

		free(H_Y);
		free(H_Z);

		//performance
		time_tally = (end.tv_sec - begin.tv_sec) + 1e-9 * (end.tv_nsec - begin.tv_nsec);
		//time = time / (double)num_vt;
		time_tally = time_tally * 1000.0;

		time = time_vote + time_tally;

		printf("n: %d    k: %d\n", num_vt, num_cnd);
		printf("Time cost per voter (micro unit): %.6lf ms\n\n", time);
	}

	return 0;
}

void ROUND1(BIG_256_56 **sk, ECP_NIST256 **pk, rsp **pi1)
{
	GENERATE_SK(sk);
	GENERATE_PK(sk, pk);
	NIZK1(pk, sk, pi1);
}

void ROUND2(ECP_NIST256 **H_Y, ECP_NIST256 **H_Z, ECP_NIST256 **Y, ECP_NIST256 **Z, ECP_NIST256 **pk, BIG_256_56 **sk, int *rnd, rsp **pi2, rsp **pi3)
{
	COMPUTE_HZ(H_Z, Z, sk, rnd);

	NIZK2(rnd, H_Z, Z, pk, sk, pi2);

	COMPUTE_HY(H_Y, Y, sk, rnd);

	NIZK3(rnd, H_Z, H_Y, Z, Y, sk, pi3);
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

void HASH_H1(BIG_256_56 out, ECP_NIST256 *cmt, int len)
{
	BIG_256_56 p;

	unsigned char domain[8] = { 'H','A','S','H','2', 0, 0, 0 };
	unsigned char element[33] = { 0 };

	unsigned char *message;
	unsigned char digest[48] = { 0 };

	octet Element = { 0, 33, element };
	octet Message;
	octet Digest = { 0, 48, digest };

	//Initialize Message
	message = (unsigned char *)malloc(4 + 33 * len);
	Message.len = 4 + 33 * len;
	Message.val = message;
	Message.max = 4 + 33 * len;

	//Fill Message
	memcpy(message, &len, sizeof(len));
	for (int i = 0; i < len; i++)
	{
		ECP_NIST256_toOctet(&Element, cmt + i, 1);
		memcpy(message + 4 + i * 33, element, 33);
	}

	HASH_CORE(&Digest, &Message, domain);
	BIG_256_56_fromBytesLen(out, Digest.val, Digest.len);

	BIG_256_56_rcopy(p, CURVE_Order_NIST256);

	BIG_256_56_mod(out, p);

	free(message);
}

rsp **SETUP_PROOF()
{
	int i, j;
	rsp **new_list;
	new_list = (rsp **)malloc(sizeof(rsp *) * num_vt);

	for (i = 0; i < num_vt; i++)
	{
		new_list[i] = (rsp *)malloc(sizeof(rsp) * num_cnd);
	}

	return new_list;
}

void NIZK1(ECP_NIST256 **pk, BIG_256_56 **sk, rsp **pi1)
{
	int i, j;

	BIG_256_56 **out = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt); //2-dim
	ECP_NIST256 *cmt = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * 3); //G, r*G, x*G

	BIG_256_56 **r = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt);
	ECP_NIST256 **cmt_r = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);
	BIG_256_56 order;
	ECP_NIST256 G;
	unsigned char seed[32] = { 0 };
	csprng CSPRNG;

	for (i = 0; i < num_vt; i++)
	{
		out[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);
		r[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);
		cmt_r[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
	}

	ECP_NIST256_generator(&G);
	RAND_seed(&CSPRNG, 32, seed);
	BIG_256_56_rcopy(order, CURVE_Order_NIST256);

	//랜덤 r 추출 -> r*G (cmt) 계산
	for (i = 0; i < num_vt; i++)
	{
		for (j = 0; j < num_cnd; j++)
		{
			/*pi1[i][j].B = (BIG_256_56 *)malloc(sizeof(BIG_256_56));
			pi1[i][j].E = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256));*/

			BIG_256_56_randomnum(r[i][j], order, &CSPRNG); //random r 추출

			ECP_NIST256_generator(&cmt_r[i][j]); //G
			ECP_NIST256_mul(&cmt_r[i][j], r[i][j]); //r*G

			ECP_NIST256_copy(&cmt[0], &G); //copy G (stmt)
			ECP_NIST256_copy(&cmt[1], &cmt_r[i][j]); //copy r*G (cmt)
			ECP_NIST256_copy(&cmt[2], &pk[i][j]); //copy x*G (stmt)

			//HASH
			HASH_H1(out[i][j], cmt, 3);

			//computes responses
			BIG_256_56_modmul(pi1[i][j].B[0], sk[i][j], out[i][j], order); //B[0] = x*c
			BIG_256_56_modneg(pi1[i][j].B[0], pi1[i][j].B[0], order); //B[0] = -x*c
			BIG_256_56_modadd(pi1[i][j].B[0], r[i][j], pi1[i][j].B[0], order); //B[0] = r - x*c

			ECP_NIST256_copy(&pi1[i][j].E[0], &cmt_r[i][j]); //E[0] = r*G
		}
	}

	//free
	for (i = 0; i < num_vt; i++)
	{
		free(r[i]);
		free(cmt_r[i]);
		free(out[i]);
	}

	free(r);
	free(cmt_r);
	free(cmt);
	free(out);
}

void NIZK2(int *rnd, ECP_NIST256 **H_Z, ECP_NIST256 **Z, ECP_NIST256 **pk, BIG_256_56 **sk, rsp **pi2)
{
	int i, j;

	BIG_256_56 **out = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt); //2-dim
	ECP_NIST256 *cmt = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * 7); //G, x, y, a1, a2, b1, b2

	ECP_NIST256 G;
	unsigned char seed[32] = { 0 };
	csprng CSPRNG;

	BIG_256_56 order;
	BIG_256_56 **r = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt);
	BIG_256_56 **s1 = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt);
	BIG_256_56 **s2 = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt);
	BIG_256_56 **d1 = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt);
	BIG_256_56 **d2 = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt);

	ECP_NIST256 **cmt_x = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt); //pk로 사용해도 무방
	ECP_NIST256 **cmt_y = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);

	ECP_NIST256 **cmt_a1 = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);
	ECP_NIST256 **cmt_a2 = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);
	ECP_NIST256 **cmt_b1 = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);
	ECP_NIST256 **cmt_b2 = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);

	for (i = 0; i < num_vt; i++)
	{
		out[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);

		r[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);
		s1[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);
		s2[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);
		d1[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);
		d2[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);

		cmt_x[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
		cmt_y[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);

		cmt_a1[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
		cmt_a2[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
		cmt_b1[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
		cmt_b2[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
	}

	for (i = 0; i < num_vt; i++)
	{
		memcpy(cmt_x[i], pk[i], sizeof(ECP_NIST256) * num_cnd); //x는 pk, 즉 x*G로 초기화
		memcpy(cmt_y[i], H_Z[i], sizeof(ECP_NIST256) * num_cnd); //y는 H_Z로 초기화
		memcpy(cmt_b1[i], Z[i], sizeof(ECP_NIST256) * num_cnd);
		memcpy(cmt_b2[i], Z[i], sizeof(ECP_NIST256) * num_cnd); //b1, b2는 for문 들어가기 전에 h로 초기화
	}

	ECP_NIST256_generator(&G);
	BIG_256_56_rcopy(order, CURVE_Order_NIST256);

	RAND_seed(&CSPRNG, 32, seed);

	for (i = 0; i < num_vt; i++)
	{
		for (j = 0; j < num_cnd; j++)
		{
			/*pi2[i][j].B = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * 4);
			pi2[i][j].E = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * 4);*/

			BIG_256_56_randomnum(r[i][j], order, &CSPRNG); //random r 추출

			if (rnd[i] != j)
			{
				BIG_256_56_randomnum(s2[i][j], order, &CSPRNG); //random s2 추출
				BIG_256_56_randomnum(d2[i][j], order, &CSPRNG); //random d2 추출

				ECP_NIST256_generator(&cmt_a1[i][j]); //cmt_a1 = G
				ECP_NIST256_mul(&cmt_a1[i][j], r[i][j]); //cmt_a1 = r*G

				ECP_NIST256_mul(&cmt_b1[i][j], r[i][j]); //cmt_b1 = r*h

				ECP_NIST256_generator(&cmt_a2[i][j]); //cmt_a2 = G
				ECP_NIST256_mul2(&cmt_a2[i][j], &cmt_x[i][j], s2[i][j], d2[i][j]); //cmt_a2 = s2*G + d2*x

				ECP_NIST256_sub(&cmt_y[i][j], &G); //cmt_y = y-g
				ECP_NIST256_mul2(&cmt_b2[i][j], &cmt_y[i][j], s2[i][j], d2[i][j]); //cmt_b2 = s2*h + d2*(y-g)

				ECP_NIST256_add(&cmt_y[i][j], &G); //cmt_y = y (HASH에 넣기 위해 H_Z로 다시 복구)
			}

			else
			{
				BIG_256_56_randomnum(s1[i][j], order, &CSPRNG); //random s1 추출
				BIG_256_56_randomnum(d1[i][j], order, &CSPRNG); //random d1 추출

				ECP_NIST256_generator(&cmt_a1[i][j]); //cmt_a1 = G
				ECP_NIST256_mul2(&cmt_a1[i][j], &cmt_x[i][j], s1[i][j], d1[i][j]); //cmt_a1 = s1*G + d1*x

				ECP_NIST256_mul2(&cmt_b1[i][j], &cmt_y[i][j], s1[i][j], d1[i][j]); //cmt_b1 = s1*h + d1*y

				ECP_NIST256_generator(&cmt_a2[i][j]); //cmt_a2 = G
				ECP_NIST256_mul(&cmt_a2[i][j], r[i][j]); //cmt_a2 = r*G
				ECP_NIST256_mul(&cmt_b2[i][j], r[i][j]); //cmt_b2 = r*h
			}

			ECP_NIST256_copy(&cmt[0], &G); //copy G (stmt)
			ECP_NIST256_copy(&cmt[1], &cmt_x[i][j]); //copy x (stmt)
			ECP_NIST256_copy(&cmt[2], &cmt_y[i][j]); //copy y (stmt)
			ECP_NIST256_copy(&cmt[3], &cmt_a1[i][j]); //copy a1 (cmt)
			ECP_NIST256_copy(&cmt[4], &cmt_a2[i][j]); //copy a2 (cmt)
			ECP_NIST256_copy(&cmt[5], &cmt_b1[i][j]); //copy b1 (cmt)
			ECP_NIST256_copy(&cmt[6], &cmt_b2[i][j]); //copy b2 (cmt)

			//HASH
			HASH_H1(out[i][j], cmt, 7);

			//computes responses d1, d2, s1, s2 순서, a1, b1, a2, b2 순서
			if (rnd[i] != j)
			{
				BIG_256_56_modneg(pi2[i][j].B[0], d2[i][j], order); //d1 = -d2
				BIG_256_56_modadd(pi2[i][j].B[0], out[i][j], pi2[i][j].B[0], order); //d1 = c - d2
				BIG_256_56_rcopy(pi2[i][j].B[1], d2[i][j]); //d2

				BIG_256_56_modmul(pi2[i][j].B[2], sk[i][j], pi2[i][j].B[0], order); //s1 = x*d1
				BIG_256_56_modneg(pi2[i][j].B[2], pi2[i][j].B[2], order); //s1 = -x*d1
				BIG_256_56_modadd(pi2[i][j].B[2], r[i][j], pi2[i][j].B[2], order); //s1 = r - x*d1
				BIG_256_56_rcopy(pi2[i][j].B[3], s2[i][j]); //s2
			}

			else
			{
				BIG_256_56_modneg(pi2[i][j].B[1], d1[i][j], order); //d2 = -d1
				BIG_256_56_modadd(pi2[i][j].B[1], out[i][j], pi2[i][j].B[1], order); //d2 = c - d1
				BIG_256_56_rcopy(pi2[i][j].B[0], d1[i][j]); //d1

				BIG_256_56_modmul(pi2[i][j].B[3], sk[i][j], pi2[i][j].B[1], order); //s2 = x*d2
				BIG_256_56_modneg(pi2[i][j].B[3], pi2[i][j].B[3], order); //s2 = -x*d2
				BIG_256_56_modadd(pi2[i][j].B[3], r[i][j], pi2[i][j].B[3], order); //s2 = r - x*d2
				BIG_256_56_rcopy(pi2[i][j].B[2], s1[i][j]); //s1
			}

			ECP_NIST256_copy(&pi2[i][j].E[0], &cmt_a1[i][j]);
			ECP_NIST256_copy(&pi2[i][j].E[1], &cmt_a2[i][j]);
			ECP_NIST256_copy(&pi2[i][j].E[2], &cmt_b1[i][j]);
			ECP_NIST256_copy(&pi2[i][j].E[3], &cmt_b2[i][j]);
		}
	}

	//free
	for (i = 0; i < num_vt; i++)
	{
		free(r[i]);
		free(s1[i]);
		free(s2[i]);
		free(d1[i]);
		free(d2[i]);

		free(cmt_x[i]);
		free(cmt_y[i]);

		free(cmt_a1[i]);
		free(cmt_a2[i]);
		free(cmt_b1[i]);
		free(cmt_b2[i]);

		free(out[i]);
	}

	free(r);
	free(s1);
	free(s2);
	free(d1);
	free(d2);

	free(cmt_x);
	free(cmt_y);

	free(cmt_a1);
	free(cmt_a2);
	free(cmt_b1);
	free(cmt_b2);

	free(out);

	free(cmt);
}

void NIZK3(int *rnd, ECP_NIST256 **H_Z, ECP_NIST256 **H_Y, ECP_NIST256 **Z, ECP_NIST256 **Y, BIG_256_56 **sk, rsp **pi3)
{
	int i, j;

	BIG_256_56 **out = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt); //2-dim
	ECP_NIST256 *cmt = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * 7); //G, h1, h2, H_Y, H_Z, cmt_a, cmt_b

	ECP_NIST256 G;
	unsigned char seed[32] = { 0 };
	csprng CSPRNG;

	BIG_256_56 order;

	BIG_256_56 **r1 = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt);
	BIG_256_56 **r2 = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt);

	ECP_NIST256 **cmt_a = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);
	ECP_NIST256 **cmt_b = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);

	for (i = 0; i < num_vt; i++)
	{
		out[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);

		r1[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);
		r2[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);

		cmt_a[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
		cmt_b[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
	}

	for (i = 0; i < num_vt; i++)
	{
		memcpy(cmt_a[i], Z[i], sizeof(ECP_NIST256) * num_cnd); //a = h1으로 초기화
		memcpy(cmt_b[i], Y[i], sizeof(ECP_NIST256) * num_cnd); //b = h2로 초기화
	}

	ECP_NIST256_generator(&G);
	BIG_256_56_rcopy(order, CURVE_Order_NIST256);

	RAND_seed(&CSPRNG, 32, seed);


	for (i = 0; i < num_vt; i++)
	{
		for (j = 0; j < num_cnd; j++)
		{
			/*pi3[i][j].B = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * 2);
			pi3[i][j].E = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * 2);*/

			BIG_256_56_randomnum(r1[i][j], order, &CSPRNG); //random r1 추출
			BIG_256_56_randomnum(r2[i][j], order, &CSPRNG); //random r2 추출

			ECP_NIST256_mul2(&cmt_a[i][j], &G, r1[i][j], r2[i][j]); // a = r1*h1 + r2*G
			ECP_NIST256_mul2(&cmt_b[i][j], &G, r1[i][j], r2[i][j]); // a = r1*h2 + r2*G

			ECP_NIST256_copy(&cmt[0], &G);
			ECP_NIST256_copy(&cmt[1], &Z[i][j]);
			ECP_NIST256_copy(&cmt[2], &Y[i][j]);
			ECP_NIST256_copy(&cmt[3], &H_Z[i][j]);
			ECP_NIST256_copy(&cmt[4], &H_Y[i][j]);
			ECP_NIST256_copy(&cmt[5], &cmt_a[i][j]);
			ECP_NIST256_copy(&cmt[6], &cmt_b[i][j]);

			//HASH
			HASH_H1(out[i][j], cmt, 7);

			//computes responses s1, s2, a, b 순서
			BIG_256_56_modmul(pi3[i][j].B[0], sk[i][j], out[i][j], order); //s1 = c*x
			BIG_256_56_modneg(pi3[i][j].B[0], pi3[i][j].B[0], order); //s1 = -c*x
			BIG_256_56_modadd(pi3[i][j].B[0], r1[i][j], pi3[i][j].B[0], order); //s1 = r1 - c*x

			if (rnd[i] != j)
			{
				BIG_256_56_rcopy(pi3[i][j].B[1], r2[i][j]); //s2 = r2
			}
			else
			{
				BIG_256_56_rcopy(pi3[i][j].B[1], out[i][j]); //s2 = c
				BIG_256_56_modneg(pi3[i][j].B[1], pi3[i][j].B[1], order); //s2 = -c
				BIG_256_56_modadd(pi3[i][j].B[1], r2[i][j], pi3[i][j].B[1], order); //s2 = r2 - c
			}

			ECP_NIST256_copy(&pi3[i][j].E[0], &cmt_a[i][j]);
			ECP_NIST256_copy(&pi3[i][j].E[1], &cmt_b[i][j]);
		}
	}

	//free
	for (i = 0; i < num_vt; i++)
	{
		free(r1[i]);
		free(r2[i]);
		free(cmt_a[i]);
		free(cmt_b[i]);
		free(out[i]);
	}

	free(r1);
	free(r2);
	free(cmt_a);
	free(cmt_b);
	free(out);
	free(cmt);
}

BIG_256_56 **SETUP_SK()
{
	int i = 0;
	BIG_256_56 **new_list;
	new_list = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt);

	for (i = 0; i < num_vt; i++)
	{
		new_list[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);
	}

	return new_list;
}

ECP_NIST256 **SETUP_PK()
{
	int i;

	ECP_NIST256 **new_list;
	new_list = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);

	for (i = 0; i < num_vt; i++)
	{
		new_list[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
	}

	return new_list;
}

void GENERATE_SK(BIG_256_56 **sk)
{
	int i, j;
	BIG_256_56 order;

	unsigned char seed[32] = { 0 };
	csprng CSPRNG;

	RAND_seed(&CSPRNG, 32, seed);
	BIG_256_56_rcopy(order, CURVE_Order_NIST256);

	for (i = 0; i < num_vt; i++)
	{
		for (j = 0; j < num_cnd; j++)
		{
			BIG_256_56_randomnum(sk[i][j], order, &CSPRNG);
		}
	}
}

void GENERATE_PK(BIG_256_56 **sk, ECP_NIST256 **pk)
{
	int i, j;

	for (i = 0; i < num_vt; i++)
	{
		for (j = 0; j < num_cnd; j++)
		{
			ECP_NIST256_generator(&pk[i][j]);

			ECP_NIST256_mul(&pk[i][j], sk[i][j]);
		}
	}
}

ECP_NIST256 **SETUP_Y()
{
	int i;

	ECP_NIST256 **new_list;
	new_list = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);

	for (i = 0; i < num_vt; i++)
	{
		new_list[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
	}

	return new_list;
}

ECP_NIST256 **SETUP_Z()
{
	int i;

	ECP_NIST256 **new_list;
	new_list = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);

	for (i = 0; i < num_vt; i++)
	{
		new_list[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
	}

	return new_list;
}

void COMPUTE_Z(ECP_NIST256 **Z, ECP_NIST256 **pk)
{
	int i, j, k;
	ECP_NIST256 *ar = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
	ECP_NIST256 tmp;

	for (i = 0; i < num_vt; i++)
	{
		for (j = 0; j < num_cnd; j++)
		{
			for (k = 0; k < num_cnd; k++)
			{
				ECP_NIST256_copy(&ar[k], &pk[i][k]);
			}

			ECP_NIST256_copy(&tmp, &ar[j]);

			for (k = 0; k < j; k++)
			{
				ECP_NIST256_neg(&ar[k]); //j 번째 앞으로는 -, 뒤로는 +
			}

			for (k = 0; k < num_cnd; k++)
			{
				if (k == j)
				{
					ECP_NIST256_sub(&ar[j], &tmp); //j 번째 자기 자신을 빼주는 역할

					continue;
				}

				ECP_NIST256_add(&ar[j], &ar[k]);
			}

			//ECP_NIST256_mul(&ar[j], sk[i][j]);
			ECP_NIST256_copy(&Z[i][j], &ar[j]);
		}
	}

	free(ar);
}

void COMPUTE_Y(ECP_NIST256 **Y, ECP_NIST256 **pk)
{
	int i, j, k;
	ECP_NIST256 *ar = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_vt);
	ECP_NIST256 tmp;

	for (i = 0; i < num_cnd; i++) //for candidates
	{
		for (j = 0; j < num_vt; j++) //for voters
		{
			for (k = 0; k < num_vt; k++)
			{
				ECP_NIST256_copy(&ar[k], &pk[k][i]);
			}

			ECP_NIST256_copy(&tmp, &ar[j]);

			for (k = 0; k < j; k++)
			{
				ECP_NIST256_neg(&ar[k]);
			}

			for (k = 0; k < num_vt; k++)
			{
				if (k == j)
				{
					ECP_NIST256_sub(&ar[j], &tmp); //j 번째 자기 자신을 빼주는 역할

					continue;
				}

				ECP_NIST256_add(&ar[j], &ar[k]);
			}

			//ECP_NIST256_mul(&ar[j], sk[j][i]);
			ECP_NIST256_copy(&Y[j][i], &ar[j]);
		}
	}

	free(ar);
}

void COMPUTE_HZ(ECP_NIST256 **H_Z, ECP_NIST256 **Z, BIG_256_56 **sk, int *rnd)
{
	int i, j;
	ECP_NIST256 G;
	ECP_NIST256_generator(&G);

	for (i = 0; i < num_vt; i++)
	{
		memcpy(H_Z[i], Z[i], sizeof(ECP_NIST256) * num_cnd);
	}

	for (i = 0; i < num_vt; i++)
	{
		for (j = 0; j < num_cnd; j++)
		{
			ECP_NIST256_mul(&H_Z[i][j], sk[i][j]);

			//투표를 수행(본 코드에서는 첫 번째 유권자(index가 0이므로)에게 투표하는 것으로 설정)
			if (j == rnd[i])
			{
				ECP_NIST256_add(&H_Z[i][j], &G);
			}
		}
	}
}

void COMPUTE_HY(ECP_NIST256 **H_Y, ECP_NIST256 **Y, BIG_256_56 **sk, int *rnd)
{
	int i, j;
	ECP_NIST256 G;
	ECP_NIST256_generator(&G);

	for (i = 0; i < num_vt; i++)
	{
		memcpy(H_Y[i], Y[i], sizeof(ECP_NIST256) * num_cnd);
	}

	for (i = 0; i < num_vt; i++)
	{
		for (j = 0; j < num_cnd; j++)
		{
			ECP_NIST256_mul(&H_Y[i][j], sk[i][j]);

			//투표를 수행(본 코드에서는 첫 번째 유권자(index가 0이므로)에게 투표하는 것으로 설정)
			if (j == rnd[i])
			{
				ECP_NIST256_add(&H_Y[i][j], &G);
			}
		}
	}
}

void SELF_TALLY(ECP_NIST256 **H_Y)
{
	int i, j, cnt, ttl;
	ECP_NIST256 G;
	ECP_NIST256 I;
	ECP_NIST256 tmp;
	ECP_NIST256 *Rslt;

	Rslt = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);

	ECP_NIST256_generator(&G); //generator
	ECP_NIST256_inf(&I); //identity value

	ttl = 0;

	//소거한 결과를 Rslt에 저장
	for (i = 0; i < num_cnd; i++)
	{
		for (j = 1; j < num_vt; j++)
		{
			ECP_NIST256_add(&H_Y[0][i], &H_Y[j][i]);
		}

		ECP_NIST256_copy(&Rslt[i], &H_Y[0][i]);
	}

	for (i = 0; i < num_cnd; i++)
	{
		ECP_NIST256_copy(&tmp, &I);
		cnt = 0;

		while (!(ECP_NIST256_equals(&tmp, &Rslt[i])))
		{
			ECP_NIST256_add(&tmp, &G);
			cnt += 1;
		}
		//printf("Candidate %d: %d\n", i + 1, cnt);
		ttl += cnt;
	}

	free(Rslt);
	//printf("Total: %d\n", ttl);
}

void CLEAR(ECP_NIST256 **E1, ECP_NIST256 **E2, ECP_NIST256 **E3, BIG_256_56 **B4, rsp **r1, rsp **r2, rsp **r3)
{
	int i, j;

	for (i = 0; i < num_vt; i++)
	{
		free(E1[i]);
		free(E2[i]);
		free(E3[i]);
		free(B4[i]);

		free(r1[i]);
		free(r2[i]);
		free(r3[i]);
	}

	free(E1);
	free(E2);
	free(E3);
	free(B4);

	free(r1);
	free(r2);
	free(r3);
}

int VERIFY_NIZK1(ECP_NIST256 **pk, rsp **pi1)
{
	int i, j;

	int rslt_cnt = 0;
	int rslt = (int)(num_vt * num_cnd);
	BIG_256_56 **out = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt); //2-dim
	ECP_NIST256 *cmt = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * 3); //G, r*G, x*G
	ECP_NIST256 G;
	ECP_NIST256 tmp;

	for (i = 0; i < num_vt; i++)
	{
		out[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);
	}

	ECP_NIST256_generator(&G);

	for (i = 0; i < num_vt; i++)
	{
		for (j = 0; j < num_cnd; j++)
		{
			ECP_NIST256_generator(&tmp);

			ECP_NIST256_copy(&cmt[0], &G); //copy G (stmt)
			ECP_NIST256_copy(&cmt[1], &pi1[i][j].E[0]); //copy r*G (cmt)
			ECP_NIST256_copy(&cmt[2], &pk[i][j]); //copy x*G (stmt)

			//HASH
			HASH_H1(out[i][j], cmt, 3);

			//compare
			ECP_NIST256_mul2(&tmp, &pk[i][j], pi1[i][j].B[0], out[i][j]);

			if (ECP_NIST256_equals(&tmp, &pi1[i][j].E[0]))
			{
				rslt_cnt += 1;
			}
		}
	}

	//free
	for (i = 0; i < num_vt; i++)
	{
		free(out[i]);
	}

	free(out);
	free(cmt);

	if (rslt_cnt == rslt)
		return 1;
	else
		return 0;
}

int VERIFY_NIZK2(ECP_NIST256 **pk, ECP_NIST256 **Z, ECP_NIST256 **H_Z, rsp **pi2)
{
	int i, j;

	int rslt = (int)(num_vt * num_cnd);
	int rslt1;
	int rslt2 = 0;

	BIG_256_56 **out = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt); //2-dim
	ECP_NIST256 *cmt = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * 7); //G, x, y, a1, a2, b1, b2
	ECP_NIST256 G;
	BIG_256_56 c;
	ECP_NIST256 a1, a2, tmp;
	ECP_NIST256 **b1 = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);
	ECP_NIST256 **b2 = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);
	BIG_256_56 order;

	for (i = 0; i < num_vt; i++)
	{
		out[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);

		b1[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
		b2[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);

		memcpy(b1[i], Z[i], sizeof(ECP_NIST256) * num_cnd);
		memcpy(b2[i], Z[i], sizeof(ECP_NIST256) * num_cnd);
	}

	ECP_NIST256_generator(&G);
	BIG_256_56_rcopy(order, CURVE_Order_NIST256);

	for (i = 0; i < num_vt; i++)
	{
		for (j = 0; j < num_cnd; j++)
		{
			rslt1 = 0;

			ECP_NIST256_generator(&a1);
			ECP_NIST256_generator(&a2);

			ECP_NIST256_copy(&cmt[0], &G); //copy G (stmt)
			ECP_NIST256_copy(&cmt[1], &pk[i][j]); //copy x (stmt)
			ECP_NIST256_copy(&cmt[2], &H_Z[i][j]); //copy y (stmt)
			ECP_NIST256_copy(&cmt[3], &pi2[i][j].E[0]); //copy a1 (cmt)
			ECP_NIST256_copy(&cmt[4], &pi2[i][j].E[1]); //copy a2 (cmt)
			ECP_NIST256_copy(&cmt[5], &pi2[i][j].E[2]); //copy b1 (cmt)
			ECP_NIST256_copy(&cmt[6], &pi2[i][j].E[3]); //copy b2 (cmt)

			//HASH
			HASH_H1(out[i][j], cmt, 7);

			//compare
			BIG_256_56_modadd(c, pi2[i][j].B[0], pi2[i][j].B[1], order); //c =? d1+d2
			if (!(BIG_256_56_comp(c, out[i][j])))
			{
				rslt1 += 1;
			}

			ECP_NIST256_mul2(&a1, &pk[i][j], pi2[i][j].B[2], pi2[i][j].B[0]); //a1 =? s1*G+d1*x
			if (ECP_NIST256_equals(&a1, &pi2[i][j].E[0]))
			{
				rslt1 += 1;
			}

			ECP_NIST256_mul2(&b1[i][j], &H_Z[i][j], pi2[i][j].B[2], pi2[i][j].B[0]); //a1 =? s1*h+d1*y
			if (ECP_NIST256_equals(&b1[i][j], &pi2[i][j].E[2]))
			{
				rslt1 += 1;
			}

			ECP_NIST256_mul2(&a2, &pk[i][j], pi2[i][j].B[3], pi2[i][j].B[1]); //a2 =? s2*G+d2*x
			if (ECP_NIST256_equals(&a2, &pi2[i][j].E[1]))
			{
				rslt1 += 1;
			}

			ECP_NIST256_copy(&tmp, &H_Z[i][j]);
			ECP_NIST256_sub(&tmp, &G);
			ECP_NIST256_mul2(&b2[i][j], &tmp, pi2[i][j].B[3], pi2[i][j].B[1]); //b2 =? s2*h+d2*(y-g)
			if (ECP_NIST256_equals(&b2[i][j], &pi2[i][j].E[3]))
			{
				rslt1 += 1;
			}

			if (rslt1 == 5)
			{
				rslt2 += 1;
			}
		}
	}

	//free
	for (i = 0; i < num_vt; i++)
	{
		free(out[i]);
		free(b1[i]);
		free(b2[i]);
	}

	free(out);
	free(b1);
	free(b2);

	free(cmt);

	if (rslt2 == rslt)
		return 1;
	else
		return 0;
}

int VERIFY_NIZK3(ECP_NIST256 **Z, ECP_NIST256 **Y, ECP_NIST256 **H_Z, ECP_NIST256 **H_Y, rsp **pi3)
{
	int i, j;

	BIG_256_56 **out = (BIG_256_56 **)malloc(sizeof(BIG_256_56 *) * num_vt); //2-dim
	ECP_NIST256 *cmt = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * 7); //G, h1, h2, H_Y, H_Z, cmt_a, cmt_b
	ECP_NIST256 G;
	ECP_NIST256 **a1 = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);
	ECP_NIST256 **a2 = (ECP_NIST256 **)malloc(sizeof(ECP_NIST256 *) * num_vt);
	ECP_NIST256 t1, t2;

	int rslt = (int)(num_vt * num_cnd);
	int rslt1;
	int rslt2 = 0;

	for (i = 0; i < num_vt; i++)
	{
		out[i] = (BIG_256_56 *)malloc(sizeof(BIG_256_56) * num_cnd);

		a1[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);
		a2[i] = (ECP_NIST256 *)malloc(sizeof(ECP_NIST256) * num_cnd);

		memcpy(a1[i], Z[i], sizeof(ECP_NIST256) * num_cnd); //a1 = h1
		memcpy(a2[i], Y[i], sizeof(ECP_NIST256) * num_cnd); //a2 = h2
	}

	ECP_NIST256_generator(&G);

	for (i = 0; i < num_vt; i++)
	{
		for (j = 0; j < num_cnd; j++)
		{
			rslt1 = 0;

			ECP_NIST256_copy(&cmt[0], &G);
			ECP_NIST256_copy(&cmt[1], &Z[i][j]);
			ECP_NIST256_copy(&cmt[2], &Y[i][j]);
			ECP_NIST256_copy(&cmt[3], &H_Z[i][j]);
			ECP_NIST256_copy(&cmt[4], &H_Y[i][j]);
			ECP_NIST256_copy(&cmt[5], &pi3[i][j].E[0]);
			ECP_NIST256_copy(&cmt[6], &pi3[i][j].E[1]);

			//HASH
			HASH_H1(out[i][j], cmt, 7);

			//compare
			ECP_NIST256_mul2(&a1[i][j], &G, pi3[i][j].B[0], pi3[i][j].B[1]); //a1 = s1*h1 + s2*G
			ECP_NIST256_copy(&t1, &H_Z[i][j]); //t1 = x*h1 + v*G
			ECP_NIST256_mul(&t1, out[i][j]); //t1 = c*t1
			ECP_NIST256_add(&a1[i][j], &t1);

			if (ECP_NIST256_equals(&a1[i][j], &pi3[i][j].E[0]))
			{
				rslt1 += 1;
			}

			ECP_NIST256_mul2(&a2[i][j], &G, pi3[i][j].B[0], pi3[i][j].B[1]); //a2 = s1*h2 + s2*G
			ECP_NIST256_copy(&t2, &H_Y[i][j]); //t2 = x*h2 + v*G
			ECP_NIST256_mul(&t2, out[i][j]); //t2 = c*t2
			ECP_NIST256_add(&a2[i][j], &t2);

			if (ECP_NIST256_equals(&a2[i][j], &pi3[i][j].E[1]))
			{
				rslt1 += 1;
			}

			if (rslt1 == 2)
				rslt2 += 1;
		}
	}

	//free
	for (i = 0; i < num_vt; i++)
	{
		free(out[i]);
		free(a1[i]);
		free(a2[i]);
	}

	free(out);
	free(a1);
	free(a2);
	free(cmt);

	if (rslt2 == rslt)
		return 1;
	else
		return 0;
}