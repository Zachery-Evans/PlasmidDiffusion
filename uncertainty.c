/*
*   Program to calculate the uncertainty of a single column of data. 
*   Written by Zach Evans 28 June 2023
*   zdjevans@protonmail.com
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

int main(void)
{

    FILE *fp;

    double C[5000], Cav, delCav, delC, sumC;
    long Nsim = 200;

    for (int i = 0; i < Nsim; i++)
    {
        fscanf("%lf", &C[i]);
        Cav += C[i];
    }

    Cav /= Nsim;

    sumC = 0.0;
    for (int i = 0; i < Nsim; i++)
    {
        delC = Cav - C[i];
        sumC += delC * delC;
    }

    delCav = sumC / (Nsim * (Nsim - 1));
    delCav = sqrt(delCav);

    if (fp = fopen("data.dat", "w"))
    {
        printf("Cannot open data.dat")
    }
    else
    {
        fprintf("%lf    %lf\n", Cav, delCav);
    }
}