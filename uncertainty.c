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

int main(int argc, char *argv[])
{
    FILE *fp;

    double C[5000], Cav=0.0, delCav=0.0, delC=0.0, sumC=0.0;
    int Nsim=0;
    char filename[500];

    // c will store each char in the file as we read them in one at a time,
    // current_line will store the running tally of how many lines we've foun
    char c;
    
    fp = fopen("x2x3cm.dat", "r");
 
    // Check if file exists
    if (fp == NULL)
    {
        printf("Could not open file %s", filename);
        return 0;
    }
 
    // Extract characters from file and store in character c
    do 
    {
    // read the next character from the file
    c = fgetc(fp);

    // if it's a newline, we've found another line
    if (c == '\n') Nsim++;
  
  // continue until the special end of file value is returned from fgetc
    } 
    while (c != EOF);
    // Close the file
    fclose(fp);
    Nsim;

    printf("Number of completed simulations: %ld\n", Nsim);

    fp = fopen("x2x3cm.dat", "r");
    
    for (int i = 0; i < Nsim; i++)
    {
        scanf("%lf", &C[i]);
        //printf("%lf\n", C[i]);
        Cav += C[i];
    }

    fclose(fp);

    Cav /= Nsim;
    //printf("Total Average: %lf\n", Cav);
    sumC = 0.0;
    for (int i = 0; i < Nsim; i++)
    {
        delC = Cav - C[i];
        sumC += delC * delC;
    }

    delCav = sumC / (Nsim * (Nsim - 1));
    delCav = sqrt(delCav);
    
    printf("%lf    %lf\n", Cav, delCav);


    if ((fp = fopen("data.dat", "w")) == NULL)
    {
        printf("Cannot open data.dat\n");
    }
    else
    {
        fprintf(fp, "%lf    %lf\n", Cav, delCav);
	fclose(fp);
    }
}
