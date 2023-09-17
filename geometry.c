#include <stdio.h>
#include <math.h>

#define PI 3.141592653589793
#define BUFFER_SIZE 1000

void input(void);

long nseg1, nseg2, nseg3, nseg4, ncyc, overlap, nacc, kk, itest, iseed;
long neq, nbintot, ibin, ichain, nsamp, nacc_shift, nshift;
long imov, kmaxtest, freq_samp, cmFreqSamp, freq_mon, freq_samp, ncmt, ngridx, ngridy;

double L, H, Ld2, Hd2, rmax, xt, yt, zt, dx, dy, dz, re, dr2, drxy2, dr2min, dr2max;
double qmin, qmax, re2av, re2, drmin, drmax, gridspace, gridspacex_real, gridspacey_real;
double amax, bmin, amax2, bmin2, ecc, Area, rectangleArea, rectangleXYRatio, xBoxMax, yBoxMax, rshift_max;
double xBoxMaxd2, yBoxMaxd2, kappa, xold, yold, zold, delphi_max;

int main(void)
{
    double ii, jj;
    input(); // Read data given in the input file

    FILE *gp;
    FILE *fPtr;
    FILE *fTemp;

    int Nsim = 0, line, count, fgetwc(FILE *stream);
    char path[100] = "geometry.xyz", buffer[BUFFER_SIZE], newline[BUFFER_SIZE], c;

    // c will store each char in the file as we read them in one at a time,
    // current_line will store the running tally of how many lines we've foun

    // --------------------------------------------------------------
    //
    // CREATING THE GEOMETRY OF THE SYSTEM IN XYZ FILE FORMAT FOR VMD
    //
    // --------------------------------------------------------------

    if ((gp = fopen("geometry.xyz", "w")) == NULL)
    { // reading mc.inp
        printf("Cannot open file: geometry.xyz\n");
        return(0);
    }
    else
    {
        if (ecc < 1.0)
        {
            amax = bmin / sqrt(1 - ecc * ecc);
        }
        else
        {
            amax = bmin;
        }
        rectangleArea = Area - PI * amax * bmin;

        if (rectangleArea < 0.0)
        {
            rectangleArea = 0.0;
        }
        yBoxMax = 2.0 * bmin;
        yBoxMaxd2 = bmin; // Width of the rectangle section equivalent to the semi-minor axis
        xBoxMax = rectangleArea / (2.0 * bmin);
        xBoxMaxd2 = xBoxMax / 2.0;
        Hd2 = H / 2.0;

        fprintf(gp, "%d\n", 1);
        fprintf(gp, "Surface:  %d\n", 0);

        for (double ii = 0.0; ii < xBoxMax; ii += 0.5)
        {
            fprintf(gp, "N    %lf  %lf  %lf\n", -xBoxMaxd2 + ii, bmin, Hd2);
            fprintf(gp, "N    %lf  %lf  %lf\n", -xBoxMaxd2 + ii, -bmin, Hd2);
        }

        if (ecc < 1.0)
        {
            for (double ii = -PI / 2; ii < PI / 2; ii += 0.01)
            {
                fprintf(gp, "N    %lf  %lf  %lf\n", xBoxMaxd2 + amax * cos(ii), bmin + bmin * sin(ii) - bmin, Hd2);
                fprintf(gp, "N    %lf  %lf  %lf\n", -xBoxMaxd2 - amax * cos(ii), -bmin + bmin * sin(ii) + bmin, Hd2);
            }
        }

        double amax2 = amax * amax, bmin2 = bmin * bmin;
        double angle = -PI / 2.0;

        for (double jj = -yBoxMaxd2; jj < yBoxMaxd2; jj += 0.5)
        {
            for (double ii = -xBoxMaxd2 - amax; ii < amax + xBoxMaxd2; ii += 0.5)
            {
                if (ii > xBoxMaxd2 && ecc < 1.0)
                { // If the polymer is outside of the rightmost semi-ellipse, write only if inside elliptical cavity
                    if ((ii - xBoxMaxd2) * (ii - xBoxMaxd2) < amax2 * (1 - (jj * jj) / bmin2) && jj < bmin2 * (1 - (ii - xBoxMaxd2) * (ii - xBoxMaxd2) / amax2))
                    {
                        fprintf(gp, "N  %lf  %lf  %lf\n", ii, jj, Hd2);
                    }
                }

                if (ii < -xBoxMaxd2 && ecc < 1.0)
                { // If the polymer is outside of the leftmost semi-ellipse, write only if inside elliptical cavity
                    if ((ii + xBoxMaxd2) * (ii + xBoxMaxd2) < amax2 * (1 - (jj * jj) / bmin2) && jj < bmin2 * (1 - (ii + xBoxMaxd2) * (ii + xBoxMaxd2) / amax2))
                    {
                        fprintf(gp, "N  %lf  %lf  %lf\n", ii, jj, Hd2);
                    }
                }

                if (ii < xBoxMaxd2 && ii > -xBoxMaxd2)
                {
                    fprintf(gp, "N  %lf  %lf  %lf\n", ii, jj, Hd2);
                }
                else if (ecc < 1.0)
                {
                    if (jj < bmin * bmin * (1 - (ii + xBoxMaxd2) * (ii + xBoxMaxd2) / amax * amax) && jj < bmin * bmin * (1 - (ii - xBoxMaxd2) * (ii - xBoxMaxd2) / amax * amax))
                    {
                        fprintf(gp, "N  %lf  %lf  %lf\n", ii, jj, Hd2);
                    }
                }
            }
        }
    }

    // --------------------------------------------------------------
    //
    // COUNTING NUMBER OF LINES IN "geometry.xyz" FILE
    //
    // --------------------------------------------------------------

    gp = fopen("geometry.xyz", "r");
    // Check if file exists
    if (gp == NULL)
    {
        printf("Could not open file %s", path);
        return(0);
    }
    // Extract characters from file and store in character c
    // Extract characters from file and store in character c
    for (c = fgetc(gp); c != EOF; c = fgetc(gp))
        if (c == '\n') // Increment count if this character is newline
            Nsim++;
    // Close the file
    fclose(gp);

    // --------------------------------------------------------------
    //
    // REPLACING NUMBER OF LINES WITH CORRECT NUMBER
    //
    // --------------------------------------------------------------

    line = Nsim;
    if ((fTemp = fopen("replace.tmp", "w")) == NULL)
    { // reading mc.inp
        printf("Cannot open file: replace.tmp\n");
        return(0);
    }
    else
    {
        if (ecc < 1.0)
        {
            amax = bmin / sqrt(1 - ecc * ecc);
        }
        else
        {
            amax = bmin;
        }
        rectangleArea = Area - PI * amax * bmin;

        if (rectangleArea < 0.0)
        {
            rectangleArea = 0.0;
        }
        yBoxMax = 2.0 * bmin;
        yBoxMaxd2 = bmin; // Width of the rectangle section equivalent to the semi-minor axis
        xBoxMax = rectangleArea / (2.0 * bmin);
        xBoxMaxd2 = xBoxMax / 2.0;
        Hd2 = H / 2.0;

        fprintf(gp, "%ld\n", line);
        fprintf(gp, "Surface:  %ld\n", 0);

        for (double ii = 0.0; ii < xBoxMax; ii += 0.5)
        {
            fprintf(gp, "N    %lf  %lf  %lf\n", -xBoxMaxd2 + ii, bmin, Hd2);
            fprintf(gp, "N    %lf  %lf  %lf\n", -xBoxMaxd2 + ii, -bmin, Hd2);
        }

        if (ecc < 1.0)
        {
            for (double ii = -PI / 2; ii < PI / 2; ii += 0.01)
            {
                fprintf(gp, "N    %lf  %lf  %lf\n", xBoxMaxd2 + amax * cos(ii), bmin + bmin * sin(ii) - bmin, Hd2);
                fprintf(gp, "N    %lf  %lf  %lf\n", -xBoxMaxd2 - amax * cos(ii), -bmin + bmin * sin(ii) + bmin, Hd2);
            }
        }

        double amax2 = amax * amax, bmin2 = bmin * bmin;
        double angle = -PI / 2.0;

        for (double jj = -yBoxMaxd2; jj < yBoxMaxd2; jj += 0.5)
        {
            for (double ii = -xBoxMaxd2 - amax; ii < amax + xBoxMaxd2; ii += 0.5)
            {
                if (ii > xBoxMaxd2 && ecc < 1.0)
                { // If the polymer is outside of the rightmost semi-ellipse, write only if inside elliptical cavity
                    if ((ii - xBoxMaxd2) * (ii - xBoxMaxd2) < amax2 * (1 - (jj * jj) / bmin2) && jj < bmin2 * (1 - (ii - xBoxMaxd2) * (ii - xBoxMaxd2) / amax2))
                    {
                        fprintf(gp, "N  %lf  %lf  %lf\n", ii, jj, Hd2);
                    }
                }

                if (ii < -xBoxMaxd2 && ecc < 1.0)
                { // If the polymer is outside of the leftmost semi-ellipse, write only if inside elliptical cavity
                    if ((ii + xBoxMaxd2) * (ii + xBoxMaxd2) < amax2 * (1 - (jj * jj) / bmin2) && jj < bmin2 * (1 - (ii + xBoxMaxd2) * (ii + xBoxMaxd2) / amax2))
                    {
                        fprintf(gp, "N  %lf  %lf  %lf\n", ii, jj, Hd2);
                    }
                }

                if (ii < xBoxMaxd2 && ii > -xBoxMaxd2)
                {
                    fprintf(gp, "N  %lf  %lf  %lf\n", ii, jj, Hd2);
                }
                else if (ecc < 1.0)
                {
                    if (jj < bmin * bmin * (1 - (ii + xBoxMaxd2) * (ii + xBoxMaxd2) / amax * amax) && jj < bmin * bmin * (1 - (ii - xBoxMaxd2) * (ii - xBoxMaxd2) / amax * amax))
                    {
                        fprintf(gp, "N  %lf  %lf  %lf\n", ii, jj, Hd2);
                    }
                }
            }
        }

        /* Delete original source file */
        remove(path);
        /* Rename temporary file as original file */
        rename("replace.tmp", path);
    }
}

// --------------------------------------------------------------
//
// COLLECT INPUT PARAMETERS FOR GRAPHING GEOMETRY
//
// --------------------------------------------------------------
void input(void)
{
    FILE *fp;

    if ((fp = fopen("mc.inp", "r")) == NULL)
    { // reading mc.inp
        printf("Cannot open file: mc.inp\n");
    }
    else
    {
        fscanf(fp, "%ld%*s", &nseg1);
        fscanf(fp, "%ld%*s", &nseg2);
        fscanf(fp, "%ld%*s", &nseg3);
        fscanf(fp, "%ld%*s", &nseg4);
        fscanf(fp, "%lf%*s", &Area);
        fscanf(fp, "%lf%*s", &bmin);
        fscanf(fp, "%lf%*s", &ecc);
        fscanf(fp, "%lf%*s", &H);
        fscanf(fp, "%lf%*s", &kappa);
        fscanf(fp, "%lf%*s", &drmin);
        fscanf(fp, "%lf%*s", &drmax);
        fscanf(fp, "%lf%*s", &gridspace);

        fscanf(fp, "\n%ld%*s", &ncyc);
        fscanf(fp, "%ld%*s", &neq);
        fscanf(fp, "%lf%*s", &rmax);
        fscanf(fp, "%lf%*s", &delphi_max);
        fscanf(fp, "%lf%*s", &rshift_max);
        fscanf(fp, "%ld%*s", &iseed);

        fscanf(fp, "\n%ld%*s", &freq_samp);
        fscanf(fp, "\n%ld%*s", &cmFreqSamp);
        fscanf(fp, "\n%ld%*s", &freq_mon);
        fscanf(fp, "\n%ld%*s", &freq_samp);

        fscanf(fp, "\n%ld%*s", &imov);
    }

    fclose(fp);
}
