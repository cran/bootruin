#include <R.h>
void rplimits(double *hl, double *psizero, int *num, double *rp)
{
    double fl[*num], fu[*num];
    double facl = *psizero / (1 - *psizero * hl[0]);

    fu[0] = (1 - *psizero);
    fl[0] = fu[0] / (1 - *psizero * hl[0]);

    *rp = 2;

    for(int i = 1; i <= *num; i++)
    {
        fl[i] = 0;
        fu[i] = 0;
        for(int j = 1; j <= i; j++)
        {
            fl[i] += hl[j] * fl[i-j];
            fu[i] += hl[j-1] * fu[i-j];
        }
        fl[i] *= facl;
        fu[i] *= *psizero;
        *rp -= ((i != *num) ? fl[i-1] : 0) + fu[i-1];
    }

    *rp *= 0.5;
}

/* void hlvec_nonp(double *data, double *hl, double *reserve, double *interval)
{
} */

/* void hlvec_logn(double *data, double *hl, double *reserve, double *interval)
{
} */

/* double moment_estimator(double *data, int k)
{
    int n = sizeof (data) / sizeof *(data);
    double mean = 0;
    for (int i = 0; i <= n; i++) {
        mean += pow(data[i], k);
    }
    return mean / (double) n;
} */

/* double average(double *data)
{
    return moment_estimator(data, 1);
} */

/* double variance(double *data)
{
    return moment_estimator(data, 2) - pow(moment_estimator(data, 1), 2);
} */

/* double stdev(double *data)
{
    return sqrt(variance(data));
} */

/* void hlvec_exp(double *datamean, double *hl, double *reserve, double *interval)
{
    int num = floor(*reserve / *interval) + 1;
    //double mean = average(data);
    double exphl[num];
    exphl[0] = 1;
    for (int i = 0; i < num; i++) {
        exphl[i+1] = exp(- *interval * (double) (i+1) / *datamean);
        hl[i] = exphl[i] - exphl[i+1];
    }
} */
