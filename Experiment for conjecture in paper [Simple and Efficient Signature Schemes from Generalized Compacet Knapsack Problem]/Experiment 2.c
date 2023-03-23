#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define REP 1000
#define PARAM_Q 19
#define PARAM_M 4
#define PARAM_N 4
#define PARAM_BETA 1

void Test_for_conjecture(int q, int m, int n, int beta);
void init(int q, int n);
int *cal_poly(int **a, int q, int m, int n, int beta, uint32_t limit);
int *del_dup(int *arr, int q, int m, int n, uint32_t limit);
int check_a(int **a, int m, int n);
void q_sort(int *arr, int l, int r);
int partition(int *arr, int l, int r);
void swap(int *p, int *q);

int *cnt_img = NULL;


int main()
{
    uint32_t q, m, n, beta;
    q = PARAM_Q;
    m = PARAM_M;
    n = PARAM_N;
    beta = PARAM_BETA;

    Test_for_conjecture(q, m, n, beta);

    return 0;
}


void Test_for_conjecture(int q, int m, int n, int beta)
{
    uint32_t limit, mn;
    int **a = NULL;
    int h, i, j, o_i;
    double k;
    int nonce = 0;
    int *arr = NULL;
    int max;
    double avg_img;
    int cnt = 0;

    srand(time(NULL));

    limit = 2 * beta + 1;
    mn = m * n;
    limit = (uint32_t)pow((float)limit, (float)mn);

    a = (int **)malloc(sizeof(int *) * m);
    for (i = 0; i < m; i++)
    {
        a[i] = (int *)malloc(sizeof(int) * n);
        memset(a[i], 0, sizeof(int) * n);
    }

    max = (int)pow((float)q, (float)n);

    avg_img = 0.0;

    init(q, n);

    for (o_i = 0; o_i < 1; o_i++)
    {

        while (1)
        {
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    a[i][j] = rand() % q;
                }
            }
            if (check_a(a, m, n))
                break;
        }

        arr = cal_poly(a, q, m, n, beta, limit);

        q_sort(arr, 0, limit - 1);

        arr = del_dup(arr, q, m, n, limit);

        avg_img += (double)arr[0];

        j = (int)pow((float)q, (float)n);

        if (arr[0] == j)
        {
            cnt += 1;
        }
        else
        {
            printf("Not ideal case: %d\n", arr[0]);
            for (h = 0; h < m; h++)
            {
                for (i = 0; i < n; i++)
                {
                    printf("%d ", a[h][i]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }

    avg_img /= (float)o_i;

    j = (int)pow((float)q, (float)n);

    printf("Average of the number of each images is as below\n");

    for (i = 0; i < j; i++)
    {
        k = (float)cnt_img[i];
        k /= (float)o_i;
        printf(" %d: %.6f\n", i, k);
    }

    for (i = 0; i < m; i++)
        free(a[i]);
    free(a);

    printf("Averga of images: %.2f / %d\n", avg_img, max);
    printf("The number of ideal case: %d / %d\n", cnt, o_i);
}


int *cal_poly(int **a, int q, int m, int n, int beta, uint32_t limit)
{
    int **tmp = NULL;
    int **r = NULL;
    int **x = NULL;
    int i, j, k;
    int o_i;
    int *arr = NULL;

    arr = (int *)malloc(sizeof(int) * limit);
    memset(arr, 0, sizeof(int) * limit);

    tmp = (int **)malloc(sizeof(int *) * m);
    x = (int **)malloc(sizeof(int *) * m);
    r = (int **)malloc(sizeof(int *) * m);
    for (i = 0; i < m; i++)
    {
        tmp[i] = (int *)malloc(sizeof(int) * (2 * n));
        x[i] = (int *)malloc(sizeof(int) * n);
		r[i] = (int *)malloc(sizeof(int) * n);
	}

	for (o_i = 0; o_i < (int)limit; o_i++)
	{
		for (i = 0; i < m; i++)
		{
			memset(tmp[i], 0, sizeof(int) * (2 * n));
			memset(r[i], 0, sizeof(int) * n);
		}

		for (i = 0; i < m; i++)
		{
			for (j = 0; j < n; j++)
            {
                x[i][j] = (rand() % (2 * beta + 1)) - beta;
            }
        }

      //multiplication
        for (k = 0; k < m; k++)
        {
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    tmp[k][i + j] += a[k][i] * x[k][j];
                }
            }
        }

        //reduction
        for (k = 0; k < m; k++)
        {
            for (i = 0; i < n; i++)
            {
                r[k][i] = tmp[k][i] - tmp[k][i + n];
            }
        }

        //addition
        for (i = 0; i < n; i++)
        {
            for (k = 1; k < m; k++)
            {
                r[0][i] += r[k][i];
            }
        }

        //modulo
        for (i = 0; i < n; i++)
        {
            r[0][i] = (r[0][i] + 2 * beta * n * m * q) % q;
        }

        for (i = 0; i < n; i++)
        {
            arr[o_i] += (r[0][i] * (int)pow((float)q, (float)i));
        }

        k = arr[o_i];
        cnt_img[k] += 1;
    }

    for (i = 0; i < m; i++)
    {
        free(tmp[i]);
        free(x[i]);
        free(r[i]);
    }
    free(tmp);
    free(r);
    free(x);

    return arr;
}


int check_a(int **a, int m, int n)
{
    int t[4];
    int i, j;
    int r = 0;

    memset(t, 0, sizeof(int) * m);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            t[i] += (a[i][j] * a[i][j]);
        }

        if (t[i] == 0)
            return 0;
    }

    return 1;
}


void init(int q, int n)
{
    int ran;

    ran = (int)pow((float)q, (float)n);

    cnt_img = (int *)malloc(sizeof(int) * ran);
    memset(cnt_img, 0, sizeof(int) * ran);
}


int *del_dup(int *arr, int q, int m, int n, uint32_t limit)
{
    int i, j;
    int *buff = NULL;
    int *n_arr = NULL;


    buff = (int *)malloc(sizeof(int) * limit);
    memset(buff, 0, sizeof(int) * limit);

    j = 0;

    for (i = 0; i < (int)limit; i++)
    {

        if (i == 0)
        {
            buff[j] = arr[i];
            j += 1;
        }

        else
        {
            buff[j] = arr[i];
            j += 1;
        }

        while (arr[i] == arr[i + 1] && i < (int)limit - 1)
        {
            i += 1;
        }
    }

    n_arr = (int *)malloc(sizeof(int) * (j + 1)); //n_arr[0] (first element) is uesd to store the size of n_arr
    memset(n_arr, 0, sizeof(int) * (j + 1));
    memcpy(&n_arr[1], buff, sizeof(int) * j);

    n_arr[0] = j;

    free(buff);
    free(arr);

    return n_arr;
}


void q_sort(int *arr, int l, int r)
{
    int  p = 0;

    if (l < r)
    {
        p = partition(arr, l, r);

        q_sort(arr, l, p - 1);
        q_sort(arr, p + 1, r);
    }
}


int partition(int *arr, int l, int r)
{
    int pivot = arr[r];
    int i = l;

    for (int j = l; j < r; j++)
    {
        if (arr[j] <= pivot)
        {
            swap(&arr[i], &arr[j]);
            i += 1;
        }
    }

    swap(&arr[i], &arr[r]);

    return i;
}


void swap(int *p, int *q)
{
    int tmp;

    tmp = *p;
    *p = *q;
    *q = tmp;
}