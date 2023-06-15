#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.141592653589793

void input(void);

long nseg1, nseg2, nseg3, nbin, i, j, k, ii, ncyc, overlap, nacc, kk, itest, iseed;
long neq, nbintot, ibin, ichain, nsamp, nacc_shift, nshift;
long imov, kmaxtest, freq_samp, cmFreqSamp, freq_mon, freq_mov, ncmt, ngridx, ngridy;

double L, H, Ld2, Hd2, rmax, xt, yt, zt, dx, dy, dz, re, dr2, drxy2, dr2min, dr2max;
double qmin, qmax, re2av, re2, drmin, drmax, gridspace, gridspacex_real, gridspacey_real;
double amax, bmin, amax2, bmin2, ecc, Area, rectangleArea, rectangleXYRatio, xBoxMax, yBoxMax, rshift_max;
double xBoxMaxd2, yBoxMaxd2;
double kappa, xold, yold, zold, delphi_max;
double z1min, z1max, z2min, z2max, zcm1, zcm2, z1bcm, z2bcm;
double **prob1, **prob2, **prob3, **probmon;

double x1, x2, x3, yone, y2, y3, z1, z2, z3;
double vax, vay, vaz, vbx, vby, vbz;
double va_sq, vb_sq;
double va_dot_vb;
double theta_old, theta_new;
double energy_new[3];
double energy_old[3];
double E_new, E_old;
double delta_E;

long ind;

long irep;
double rep_prob;
long nacc_rep;
long nrep;
double phi_prime, costheta_prime, dr_prime;
double dx_prime, dy_prime, dz_prime;
double dx_fixed, dy_fixed, dz_fixed;
double alpha, beta;
double cos_alpha, cos_beta, sin_alpha, sin_beta;
double ux, uy, uz;
double u, uxy;

int main(void)
{
    input();

    FILE *gp;
    if ((gp = fopen("geometry.xyz", "w")) == NULL)
    { // reading mc.inp
        printf("Cannot open file: geometry.dat\n");
    }
    else
    {
        amax = bmin / sqrt(1 - ecc * ecc);
        rectangleArea = Area - PI * amax * bmin;
        yBoxMax = 2.0 * bmin;
        yBoxMaxd2 = bmin; // Width of the rectangle section equivalent to the semi-minor axis
        xBoxMax = rectangleArea / (2.0 * bmin);
        xBoxMaxd2 = xBoxMax / 2.0;
        Hd2 = H / 2.0;

        fprintf(gp, "%ld\n", 50716);
        fprintf(gp, "Surface:  %ld\n", 0);

        for (double i = 0.0; i < xBoxMax; i+=0.5)
        {
            for (double j = 0.0; j < H; j += 0.5)
            {
                fprintf(gp, "N    %lf  %lf  %lf\n", -xBoxMaxd2 + i, bmin, -Hd2 + j);
                fprintf(gp, "N    %lf  %lf  %lf\n", -xBoxMaxd2 + i, -bmin, -Hd2 + j);
            }

            fprintf(gp, "N    %lf  %lf  %lf\n", -xBoxMaxd2 + i, bmin, -Hd2);
            fprintf(gp, "N    %lf  %lf  %lf\n", -xBoxMaxd2 + i, -bmin, -Hd2);
        }

        for (double i = -PI / 2; i < PI / 2; i += 0.01)
        {

            for (double j = 0.0; j < H; j += 0.5)
            {
                fprintf(gp, "N    %lf  %lf  %lf\n", -xBoxMaxd2 - amax * cos(i), bmin + bmin * sin(i) - bmin, -Hd2 + j);
                fprintf(gp, "N    %lf  %lf  %lf\n", xBoxMaxd2 + amax * cos(i), bmin + bmin * sin(i) - bmin, -Hd2 + j);
            }

            for (double i = PI / 2; i > -PI / 2; i -= 0.1)
            {
                fprintf(gp, "N    %lf  %lf  %lf\n", xBoxMaxd2 + amax * cos(i), bmin + bmin * sin(i) - bmin, Hd2);
                fprintf(gp, "N    %lf  %lf  %lf\n", xBoxMaxd2 + amax * cos(i), -bmin + bmin * sin(i) + bmin, -Hd2);
            }
        }
    }
}

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
        fscanf(fp, "\n%ld%*s", &freq_mov);

        fscanf(fp, "\n%ld%*s", &imov);
    }

    fclose(fp);
}