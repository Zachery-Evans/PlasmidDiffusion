/*
 *   Program to calculate the total average and uncertainty of a single column of data.
 *   Written by Zach Evans 28 June 2023
 *   zdjevans@protonmail.com
 *
 *
 *   Program should be compiled to a run file by: cc -lm -O3 uncertainty.c -o unc.run
 *   or similar.
 *
 *   The program should then be run as: ./unc.run < data_file.dat
 *   The program will then output a file with the total average in the first column, and the uncertainty
 *   of the average in the second column.
 *
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
        scanf("%lf%*s", &C[i]);
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

    if ((fp = fopen("data.dat", "w")) == NULL)
    {
        printf("Cannot open data.dat\n");
    }
    else
    {
        fprintf(fp, "%lf    %lf\n", Cav, delCav);
    }
}