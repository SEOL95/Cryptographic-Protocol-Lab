#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define REP 5000
#define PARAM_Q 5
#define PARAM_M 2
#define PARAM_N 4
#define PARAM_BETA 1

void Test_for_conjecture(int q, int m, int n, int beta);
void init(int q, int n);
void cal_poly(int **a, int q, int m, int n, int beta, uint32_t limit);
int *del_dup(int *arr, int q, int m, int n, uint32_t limit);
int check_a(int **a, int m, int n);
void q_sort(int *arr, int l, int r);
int partition(int *arr, int l, int r);
void swap(int *p, int *q);

int *cnt_img = NULL;
int ideal_cnt;
double avg_img;

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
    int i, j, o_i;
    double k;
    int nonce = 0;
    int *arr = NULL;
    int max;

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

    ideal_cnt = 0;

    max = (int)pow((float)q, (float)n);

    avg_img = 0.0;

    init(q, n);

    for (o_i = 0; o_i < REP; o_i++)
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

        cal_poly(a, q, m, n, beta, limit);
    }

    avg_img /= (float)o_i;

    j = (int)pow((float)q, (float)n);

    printf("Average of the number of each images is as below\n");

    for (i = 0; i < j; i++)
    {
        k = (float)cnt_img[i];
        k /= (float)ideal_cnt;
        printf(" %d: %.6f\n", i, k);
    }

    free(cnt_img);

    for (i = 0; i < m; i++)
        free(a[i]);
    free(a);

    printf("Average of images: %.2f / %d\n", avg_img, max);
    printf("The number of ideal case: %d / %d\n", ideal_cnt, o_i);
}


void init(int q, int n)
{
    int ran;

    ran = (int)pow((float)q, (float)n);

    cnt_img = (int *)malloc(sizeof(int) * ran);
    memset(cnt_img, 0, sizeof(int) * ran);
}


void cal_poly(int **a, int q, int m, int n, int beta, uint32_t limit)
{
    int **tmp = NULL;
    int **r = NULL;
    int **tmp_x = NULL;
    int ***x = NULL;
    int i, j, k, o_i;
    int *arr = NULL;
    int *p = NULL;
    int neg_beta = (-1) * beta;
    int i_1, i_2, i_3, i_4, i_5, i_6, i_7, i_8;

    arr = (int *)malloc(sizeof(int) * limit);
    memset(arr, 0, sizeof(int) * limit);

    tmp = (int **)malloc(sizeof(int *) * m);
    for (i = 0; i < m; i++)
    {
        tmp[i] = (int *)malloc(sizeof(int) * (2 * n));
    }

    r = (int **)malloc(sizeof(int *) * m);
    tmp_x = (int **)malloc(sizeof(int *) * m);
    for (i = 0; i < m; i++)
    {
        r[i] = (int *)malloc(sizeof(int) * n);
        tmp_x[i] = (int *)malloc(sizeof(int) * n);

        memset(tmp_x[i], 0, sizeof(int) * n);
    }

    x = (int ***)malloc(sizeof(int **) * limit);
    for (i = 0; i < (int)limit; i++)
        x[i] = (int **)malloc(sizeof(int *) * m);
    for (i = 0; i < (int)limit; i++)
    {
        for (j = 0; j < m; j++)
        {
            x[i][j] = (int *)malloc(sizeof(int) * n);
            memset(x[i][j], 0, sizeof(int) * n);
        }
    }

    o_i = 0;

    for (i_8 = neg_beta; i_8 <= beta; i_8++)
    {
        tmp_x[1][0] = i_8;

        for (i_7 = neg_beta; i_7 <= beta; i_7++)
        {
            tmp_x[1][1] = i_7;

            for (i_6 = neg_beta; i_6 <= beta; i_6++)
            {
                tmp_x[1][2] = i_6;

                for (i_5 = neg_beta; i_5 <= beta; i_5++)
                {
                    tmp_x[1][3] = i_5;

                    for (i_4 = neg_beta; i_4 <= beta; i_4++)
                    {
                        tmp_x[0][0] = i_4;

                        for (i_3 = neg_beta; i_3 <= beta; i_3++)
                        {
                            tmp_x[0][1] = i_3;

                            for (i_2 = neg_beta; i_2 <= beta; i_2++)
                            {
                                tmp_x[0][2] = i_2;

                                for (i_1 = neg_beta; i_1 <= beta; i_1++)
                                {
                                    tmp_x[0][3] = i_1;

                                    for (i = 0; i < m; i++)
                                    {
                                        for (j = 0; j < n; j++)
                                        {
                                            x[o_i][i][j] = tmp_x[i][j];
                                        }
                                    }
                                    o_i += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (o_i = 0; o_i < (int)limit; o_i++)
    {

        for (i = 0; i < m; i++)
        {
            memset(tmp[i], 0, sizeof(int) * (2 * n));
            memset(r[i], 0, sizeof(int) * n);
        }

        //multiplication
        for (k = 0; k < m; k++)
        {
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    tmp[k][i + j] += a[k][i] * x[o_i][k][j];
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
    }

    q_sort(arr, 0, limit - 1);
    p = del_dup(arr, q, m, n, limit);

    avg_img += p[0];

    k = (int)pow((float)q, (float)n);

    if (p[0] == k)
    {
        ideal_cnt += 1;

        for (i = 0; i < (int)limit; i++)
        {
            j = arr[i];
            cnt_img[j] += 1;
        }
    }

    free(p);
    free(arr);
    for (i = 0; i < m; i++)
        free(tmp[i]);
    free(tmp);
    for (i = 0; i < m; i++)
    {
        free(r[i]);
        free(tmp_x[i]);
    }
    free(r);
    free(tmp_x);
    for (i = 0; i < (int)limit; i++)
    {
        for (j = 0; j < m; j++)
        {
            free(x[i][j]);
        }
    }
    for (i = 0; i < (int)limit; i++)
        free(x[i]);
    free(x);
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

    return n_arr;
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
