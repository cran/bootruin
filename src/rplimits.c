#include <R.h>
void rplimits(double *hl, double *hu, double *psizero, int *num, double *rp)
{
    // Define temporary arrays containing the probablities of the bounds
    double fl[*num], fu[*num];

    // Allocate memory and initialize to 0
    memset(fl, 0.0, *num * *num * sizeof(double));
    memset(fu, 0.0, *num * *num * sizeof(double));

    // Define a constant for later re-use
    const double facl = *psizero / (1.0 - *psizero * hl[0]);

    // Initialize the first array elements
    fu[0] = (1.0 - *psizero);
    fl[0] = fu[0] / (1.0 - *psizero * hl[0]);

    // Set the probability of ruin to 2 (sic!)
    // That we can easily compute the cumulative sums and only have to divide
    // once in order to average the bounds.
    *rp = 2.0;

    // Here comes the magic: Recursion algorithm
    //for(int i = 1; i <= *num; i++)
    for(int i = 1; i < *num; i++)
    {
        //fl[i] = 0.0;
        //fu[i] = 0.0;
        //for(int j = 1; j <= i; j++)
        for(int j = 0; j < i; j++)
        {
            //fl[i] += hl[j]   * fl[i-j];
            //fu[i] += hl[j-1] * fu[i-j];
            fl[i] += hl[j] * fl[i-j];
            fu[i] += hu[j] * fu[i-j];
        }
        fl[i] *= facl;
        fu[i] *= *psizero;
        //*rp -= ((i != *num) ? fl[i-1] : 0) + fu[i-1];
        *rp -= fl[i-1] + fu[i-1];
    }

    *rp *= 0.5;
}
