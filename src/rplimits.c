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
