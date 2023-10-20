/*
This code was written by Dr. James Polson, it has been modified by Zach Evans during the summer of 2023

Email: zdjevans@protonmail.com

This program was written to study the equilibrium behaviour of plasmids confined to a dual pit geometry of an elliptically capped rectangle interacting with
T4 polymer (Much longer linear polymer).

This program was built with the Compute Canada clusters, as such it should include an input file whose contents are read from the function "void input(void)".

The input function reads from a file labeled "mc.inp" short for "Monte Carlo Input". The file should be of similar form to:

1000              nseg1   // Number of Monomers in Linear Polymer
25                nseg2   // Number of Monomers in Plasmid 1
25                nseg3   // Number of Monomers in Plasmid 2
0                 nseg4   // Number of Monomers in Plasmid 3
3000.0            Area    // Total area of geometry, note it is possible to make polymer placements impossible
12.0              bmin    // Length of semi minor axis in monomer widths
0.8               ecc     // Eccentricy of elliptical caps

15.0              H       // Lateral height of confinement ( Height in the z direction )
0.0               kappa   // Bending rigidity of polymer & plasmids (physical relevance from 0.0 -> 10.0)
0.9               drmin
1.1               drmax
1.0               gridspace

50000             ncyc    // Number of Monte Carlo Cycles
10000             neq     // Number of Equilibrium Cycles
0.14              rmax
3.14159           delphi_max
0.5               rshift_max
1348378           iseed   // Seed for Random Number Generation

100               freq_samp  // Sample frequency for probability distributions (1 sample for every 100 MC Cycles)
100               cmFreqSamp // Sample frequency for the cm correlation of plasmids

0                 imov       // 1 for writing xyz files for VMD movies
1                 plasRigid  // Plasmid rigidity
0                 xcmPrint   // Write center mass data for singular plasmids (if 1, write. else do not write)
0                 ycmPrint

Unless otherwise edited.

There are some pecularities that were implemented into this program in order to be most effective for file management as well as ease of use:

  1) If the variable "ecc" (eccentricity of the elliptical caps) is greater than or equal to 1.0, then the geometry is changed into a rectangle.

  2) There is a variable called "plasRigid" that controls whether or not the plasmids have ridigity equal to the linear polymer. If plasRigid == 1 then the plasmiids
     have rigidity. They have zero rigidity otherwise.

  3) The variable imov determines whether or not the "xyz" files used for VMD visualizations are printed to a file. It should be noted that this variable needs to be
     set to 0 if doing any work in the Compute Canada clusters, as there is not enough storage in the cloud to support printing these files anywhere but on the local
     computer.

  4) This code was generalized around adding more plasmids to the system, and was optimized such that if either the polymer or any plasmids were removed from the system
     (nseg* = 0) then no modifications would be necessary for the program to compile, run, and produce relevant data.

  5) It is not always the case that measurements for the CM of the plasmids are necessary, to avoid using too much file space on the clusters, input parameters xcmPrint
     and ycmPrint where created in order for user input to determine what particular data should be collected. It should be noted that the x2CM * x3CM etc. data will always
     be collected, as the equilibriumn behaviour of the plasmids is one way of making sure that the data you receive is within expectation.

*/

#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.141592653589793
#define NR_END 1
#define FREE_ARG char *

// function headers
double ran3(void);
void init_pos(void);
void init_pos_circular(void);
void write_log(void);
void write_data(void);
void input(void);
int check_accept(double[], double[], double[], long);
int check_shift_chain(double[], double[], double[], long);

int check_poly_energy(double[], double[], double[], long);
int check_plasmid_energy(double[], double[], double[], long);
double calc_cosine(int, int, int, double[], double[], double[]);

double **dmatrix(long, long, long, long);
void free_dmatrix(double **, long, long, long, long);

void stateCheckSingle(double);
void stateCheckDouble(double, double);
void stateCheckTriple(double, double, double);

void reptation_move_chain1(void);

void shift_move_chain(void);
void shift_move_plasmid(double[], double[], double[], long);

int check_accept_reptation(double[], double[], double[], long, long);
void calc_delta_xyz(void);

void crank_move_polymer(double[], double[], double[]);
void crank_move_plasmid(double[], double[], double[], long);

long nseg1, nseg2, nseg3, nseg4, i, j, k, ii, ncyc, overlap, nacc, kk, iseed;
long neq, ichain, nsamp, nacc_shift, nshift, xcmPrint, ycmPrint;
long imov, plasRigid, kmaxtest, freq_samp, cmFreqSamp, ngridx, ngridy;

double L, H, Ld2, Hd2, rmax, dx, dy, dz, dr2;
double drmin, drmax, gridspace, gridspacex_real, gridspacey_real;
double amax, bmin, amax2, bmin2, ecc, Area, rectangleArea, xBoxMax, yBoxMax, rshift_max;
double xBoxMaxd2, yBoxMaxd2;
double kappa, xold, yold, zold, delphi_max;
double **prob1, **plas, **probmon;

FILE *fpmov;

double r1x[5000], r1y[5000], r1z[5000];
double r2x[5000], r2y[5000], r2z[5000];
double r3x[5000], r3y[5000], r3z[5000];
double r4x[5000], r4y[5000], r4z[5000];
double x2cm[10000], x3cm[10000], x4cm[10000];
double y2cm[10000], y3cm[10000], y4cm[10000];
double plas12[10000], plas23[10000], plas13[10000];

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
long iter = 0, state[6] = {0, 0, 0, 0, 0, 0};
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

// ------------------------------------------------------------------------
// main function
// ------------------------------------------------------------------------
int main()
{
  long indx, indy;
  double xcm1, ycm1, xcm2, ycm2, xcm3, ycm3, xcm4, ycm4;
  clock_t start, end;

  input();

  rep_prob = 0.95;
  // defining the area of the ellipse as the total area subtracted by the area
  // of the rectangle now between two halves of the ellipse.

  // Requires that the y direction of the box is the same height as the semi-minor
  // axis of the ellipse.
  if (ecc < 1.0)
  {
    amax = bmin / sqrt(1 - ecc * ecc);
  }

  amax2 = amax * amax;
  bmin2 = bmin * bmin;
  rectangleArea = Area - PI * amax * bmin;
  if (rectangleArea < 0.0)
  {
    rectangleArea = 0.0;
  }

  yBoxMax = 2.0 * bmin;
  yBoxMaxd2 = bmin; // Width of the rectangle section equivalent to the semi-minor axis
  xBoxMax = rectangleArea / (2.0 * bmin);
  xBoxMaxd2 = xBoxMax / 2.0;

  //  printf("%lf \t %lf\n", amax, bmin);
  //  printf("Length of the box: %lf\n", xBoxMax);
  //  printf("1/2 Length of the box: %lf\n", xBoxMaxd2);
  //  printf("Semi-major axis: %lf\n", amax);
  //  printf("Semi-minor axis: %lf\n", bmin);
  //  printf("Height of box: %lf\n", yBoxMax);
  Hd2 = H / 2.0;

  if (rectangleArea > 0.0)
  {
    ngridx = 2.0 * (amax + xBoxMaxd2) / gridspace + 0.00001;
  }
  else
  {
    ngridx = 2.0 * (amax + xBoxMaxd2) / gridspace + 0.00001;
  }

  ngridy = 2.0 * bmin / gridspace + 0.00001;
  if (rectangleArea > 0.0)
  {
    gridspacex_real = 2.0 * (amax + xBoxMaxd2) / ngridx;
  }
  else
  {
    gridspacex_real = 2.0 * amax / ngridx;
  }
  gridspacey_real = 2.0 * bmin / ngridy;
  // printf("ngridx = %ld, ngridy = %ld, gridspacex_real = %lf, gridspacey_real = %lf\n",
  // ngridx, ngridy, gridspacex_real, gridspacey_real);
  prob1 = dmatrix(0, ngridx - 1, 0, ngridy - 1);
  plas = dmatrix(0, ngridx - 1, 0, ngridy - 1);
  probmon = dmatrix(0, ngridx - 1, 0, ngridy - 1);

  for (i = 0; i < ngridx; i++)
  {
    for (j = 0; j < ngridy; j++)
    {
      prob1[i][j] = 0.0;
      plas[i][j] = 0.0;
      probmon[i][j] = 0.0;
    }
  }

  nsamp = 0;

  write_log();

  // *************************************************************

  if (imov == 1)
  { // Don't include this in cluster
    fpmov = fopen("chain.xyz", "w");
    start = clock();
  }

  if (plasRigid == 1)
  {
    init_pos(); // function call
  }
  else
  {
    init_pos_circular();
  }

  // Don't include if statement below in cluster
  if (imov == 1)
  {
    if (ii % freq_samp == 0 && ii > neq)
    {
      fprintf(fpmov, "%ld\n", nseg1 + nseg2 + nseg3 + nseg4);
      fprintf(fpmov, "Polymer:  %ld\n", ii);

      for (i = 0; i < nseg1; i++)
      {
        fprintf(fpmov, "N    %lf  %lf  %lf\n", r1x[i], r1y[i], r1z[i]);
      }
      for (i = 0; i < nseg2; i++)
      {
        fprintf(fpmov, "O    %lf  %lf  %lf\n", r2x[i], r2y[i], r2z[i]);
      }

      for (i = 0; i < nseg3; i++)
      {
        fprintf(fpmov, "F    %lf  %lf  %lf\n", r3x[i], r3y[i], r3z[i]);
      }

      for (i = 0; i < nseg4; i++)
      {
        fprintf(fpmov, "B    %lf  %lf  %lf\n", r4x[i], r4y[i], r4z[i]);
      }
    }
  }

  nacc = 0; // will hold number of accepted moves
  nacc_rep = 0;
  nacc_shift = 0;
  nshift = 0;
  nrep = 0;

  for (ii = 0; ii < ncyc; ii++)
  {
    // if (ii % 100 == 0) printf("ii = %ld\n",ii);
    for (j = 0; j < nseg1 + nseg2 + nseg3 + nseg4; j++)
    {
      k = (nseg1 + nseg2 + nseg3 + nseg4) * ran3();

      if (k < nseg1 && nseg1 != 0)
      {
        ichain = 1;
        kmaxtest = nseg1 - 2;
      }
      else if (nseg1 < k && k < (nseg1 + nseg2) && nseg2 != 0)
      {
        ichain = 2;
        k -= nseg1;
        kmaxtest = nseg2 - 2;
      }

      else if ((nseg1 + nseg2) < k && k < (nseg1 + nseg2 + nseg3) && nseg3 != 0)
      {
        ichain = 3;
        k -= nseg1;
        k -= nseg2;
        kmaxtest = nseg3 - 2;
      }

      else if ((nseg1 + nseg2 + nseg3) < k && nseg4 != 0)
      {
        ichain = 4;
        k -= nseg1;
        k -= nseg2;
        k -= nseg3;
        kmaxtest = nseg4 - 2;
      }

      if (ran3() >= rep_prob && (k >= 2 || k < kmaxtest))
      {

        if (ichain == 1)
        {
          xold = r1x[k];
          yold = r1y[k];
          zold = r1z[k];
          crank_move_polymer(r1x, r1y, r1z);
          overlap = check_accept(r1x, r1y, r1z, nseg1);
        }
        else if (ichain == 2)
        {
          xold = r2x[k];
          yold = r2y[k];
          zold = r2z[k];
          crank_move_plasmid(r2x, r2y, r2z, nseg2);
          overlap = check_accept(r2x, r2y, r2z, nseg2);
        }
        else if (ichain == 3)
        {
          xold = r3x[k];
          yold = r3y[k];
          zold = r3z[k];
          crank_move_plasmid(r3x, r3y, r3z, nseg3);
          overlap = check_accept(r3x, r3y, r3z, nseg3);
        }
        else if (ichain == 4)
        {
          xold = r4x[k];
          yold = r4y[k];
          zold = r4z[k];
          crank_move_plasmid(r4x, r4y, r4z, nseg4);
          overlap = check_accept(r4x, r4y, r4z, nseg4);
        }

        if (overlap == 0)
        {
          nacc += 1;
        }
        else
        {
          if (ichain == 1)
          {
            r1x[k] = xold;
            r1y[k] = yold;
            r1z[k] = zold;
          }
          else if (ichain == 2)
          {
            r2x[k] = xold;
            r2y[k] = yold;
            r2z[k] = zold;
          }
          else if (ichain == 3)
          {
            r3x[k] = xold;
            r3y[k] = yold;
            r3z[k] = zold;
          }
          else if (ichain == 4)
          {
            r4x[k] = xold;
            r4y[k] = yold;
            r4z[k] = zold;
          }
        }
      }
      else
      {
        nrep++;
        if (ichain == 1)
        {
          reptation_move_chain1();
        }
        else if (ichain == 2)
        {
          shift_move_plasmid(r2x, r2y, r2z, nseg2);
        }
        else if (ichain == 3)
        {
          shift_move_plasmid(r3x, r3y, r3z, nseg3);
        }
        else if (ichain == 4)
        {
          shift_move_plasmid(r4x, r4y, r4z, nseg4);
        }
      }
    }

    if (ii % freq_samp == 0 && ii > neq)
    {

      xcm1 = 0.0;
      ycm1 = 0.0;
      for (i = 0; i < nseg1; i++)
      {
        xcm1 += r1x[i];
        ycm1 += r1y[i];
      }
      xcm1 /= nseg1;
      ycm1 /= nseg1;

      if (rectangleArea > 0.0)
      {
        indx = (xcm1 + amax + xBoxMaxd2) / gridspacex_real;
      }
      else
      {
        (xcm1 + amax) / gridspacex_real;
      }
      indy = (ycm1 + bmin) / gridspacey_real;
      /*
       * printf statements after the following conditions are met:
       * if (indx >= ngridx || indy >= ngridy)
       * Warning for invalid data
       */
      if (indx >= ngridx || indy >= ngridy)
      {
        printf("1:  indx = %ld/%ld, indy = %ld/%ld\n", indx, ngridx, indy, ngridy);
      }

      if (indx >= 0 && indx < ngridx && indy >= 0 && indy < ngridy)
      {
        prob1[indx][indy] += 1.0;
      }

      for (i = 0; i < nseg1; i++)
      {
        if (rectangleArea > 0.0)
        {
          indx = (r1x[i] + amax + xBoxMaxd2) / gridspacex_real;
        }
        else
        {
          indx = (r1x[i] + amax) / gridspacex_real;
        }
        indy = (r1y[i] + bmin) / gridspacey_real;
        if (indx >= 0 && indx < ngridx && indy >= 0 && indy < ngridy)
        {
          probmon[indx][indy] += 1.0;
        }
      }

      xcm2 = 0.0;
      ycm2 = 0.0;
      for (i = 0; i < nseg2; i++)
      {
        xcm2 += r2x[i];
        ycm2 += r2y[i];
      }
      xcm2 /= nseg2;
      ycm2 /= nseg2;

      if (rectangleArea > 0.0)
      {
        indx = (xcm2 + amax + xBoxMaxd2) / gridspacex_real;
      }
      else
      {
        indx = (xcm2 + amax) / gridspacex_real;
      }

      indy = (ycm2 + bmin) / gridspacey_real;

      if (indx >= ngridx || indy >= ngridy)
      {
        printf("2:  indx = %ld/%ld, indy = %ld/%ld\n", indx, ngridx, indy, ngridy);
        printf("    xcm2 = %lf, 2*amax = %lf,  ycm2 = %lf, 2*bmin = %lf\n",
               xcm2 + amax / 2.0, 2 * amax, ycm2 + bmin / 2.0, 2 * bmin);
      }
      if (indx >= 0 && indx < ngridx && indy >= 0 && indy < ngridy)
      {
        plas[indx][indy] += 1.0;
      }

      xcm3 = 0.0;
      ycm3 = 0.0;
      for (i = 0; i < nseg3; i++)
      {
        xcm3 += r3x[i];
        ycm3 += r3y[i];
      }
      xcm3 /= nseg3;
      ycm3 /= nseg3;

      if (rectangleArea > 0.0)
      {
        indx = (xcm3 + amax + xBoxMaxd2) / gridspacex_real;
      }
      else
      {
        indx = (xcm3 + amax) / gridspacex_real;
      }
      indy = (ycm3 + bmin) / gridspacey_real;

      if (indx >= ngridx || indy >= ngridy)
      {
        printf("3:  indx = %ld/%ld, indy = %ld/%ld\n", indx, ngridx, indy, ngridy);
        printf("    xcm3 = %lf, 2*amax = %lf,  ycm3 = %lf, 2*bmin = %lf\n",
               xcm3 + amax / 2.0, 2 * amax, ycm3 + bmin / 2.0, 2 * bmin);
      }

      if (indx >= 0 && indx < ngridx && indy >= 0 && indy < ngridy)
      {
        plas[indx][indy] += 1.0;
      }

      xcm4 = 0.0;
      ycm4 = 0.0;
      for (i = 0; i < nseg4; i++)
      {
        xcm4 += r4x[i];
        ycm4 += r4y[i];
      }
      xcm4 /= nseg4;
      ycm4 /= nseg4;

      if (rectangleArea > 0.0)
      {
        indx = (xcm4 + amax + xBoxMaxd2) / gridspacex_real;
      }
      else
      {
        indx = (xcm4 + amax) / gridspacex_real;
      }
      indy = (ycm4 + bmin) / gridspacey_real;

      if (indx >= ngridx || indy >= ngridy)
      {
        printf("4:  indx = %ld/%ld, indy = %ld/%ld\n", indx, ngridx, indy, ngridy);
        printf("    xcm4 = %lf, 2*amax = %lf,  ycm4 = %lf, 2*bmin = %lf\n",
               xcm4 + amax / 2.0, 2 * amax, ycm4 + bmin / 2.0, 2 * bmin);
      }

      if (indx >= 0 && indx < ngridx && indy >= 0 && indy < ngridy)
      {
        plas[indx][indy] += 1.0;
      }

      nsamp += 1;
    }

    if (ii % cmFreqSamp == 0 && ii > neq)
    {
      long thing = (ii - neq) / cmFreqSamp;
      // printf("%ld   %lf\n", thing, xcm2);
      iter++;
      plas23[thing] = xcm3 * xcm4;
      plas13[thing] = xcm2 * xcm4;
      plas12[thing] = xcm2 * xcm3;

      x2cm[thing] = xcm2;
      y2cm[thing] = ycm2;

      x3cm[thing] = xcm3;
      y3cm[thing] = ycm3;

      x4cm[thing] = xcm4;
      y4cm[thing] = ycm4;

      if (nseg2 != 0 && nseg3 == 0 && nseg4 == 0)
      {
        stateCheckSingle(xcm2);
      }
      else if (nseg2 != 0 && nseg3 != 0 && nseg4 == 0)
      {
        stateCheckDouble(xcm2, xcm3);
      }
      else if (nseg2 != 0 && nseg3 != 0 && nseg4 != 0)
      {
        stateCheckTriple(xcm2, xcm3, xcm4);
      }
    }

    if (imov == 1)
    {
      if (ii % freq_samp == 0 && ii > -1)
      {
        fprintf(fpmov, "%ld\n", nseg1 + nseg2 + nseg3 + nseg4);
        fprintf(fpmov, "Polymer:  %ld\n", ii);

        for (i = 0; i < nseg1; i++)
        {
          fprintf(fpmov, "N    %lf  %lf  %lf\n", r1x[i], r1y[i], r1z[i]);
        }
        for (i = 0; i < nseg2; i++)
        {
          fprintf(fpmov, "O    %lf  %lf  %lf\n", r2x[i], r2y[i], r2z[i]);
        }

        for (i = 0; i < nseg3; i++)
        {
          fprintf(fpmov, "F    %lf  %lf  %lf\n", r3x[i], r3y[i], r3z[i]);
        }

        for (i = 0; i < nseg4; i++)
        {
          fprintf(fpmov, "B    %lf  %lf  %lf\n", r4x[i], r4y[i], r4z[i]);
        }
      }
    }
  }

  printf("Acc. ratio = %lf\n", 1.0 * nacc / ((ncyc * (nseg1 + nseg2 + nseg3 + nseg4)) - nrep));
  printf("Number of reptation attempts = %ld\n", nrep);
  printf("Reptation Acc. ratio = %lf\n", 1.0 * nacc_rep / nrep);
  printf("Shift 2 Acc. ratio = %ld / %ld = %lf\n", nacc_shift, nshift, 1.0 * nacc_shift / nshift);

  write_data();

  // printf("%ld\t%ld\t%ld\n", leftEllipse, rightEllipse, centerBox);

  if (imov == 1)
  {
    fclose(fpmov);
    end = clock();
    double duration = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Program finished in %lf seconds\n", duration);
  }
}

// ----------------------------------------------------------------------
// Initializes several global variables using mc.inp (formatted into 2
// columns). The file, mc.inp is left unchanged
// ----------------------------------------------------------------------
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
    fscanf(fp, "%ld%*s", &cmFreqSamp);

    fscanf(fp, "\n%ld%*s", &imov);
    fscanf(fp, "%ld%*s", &plasRigid);
    fscanf(fp, "%ld%*s", &xcmPrint);
    fscanf(fp, "%ld%*s", &ycmPrint);
  }

  fclose(fp);
}

// ----------------------------------------------------------------------
// This function writes out the parameters which were inputted into the
// program (same as input file)
// ----------------------------------------------------------------------
void write_log(void)
{

  printf("nseg1    %ld\n", nseg1);
  printf("nseg2    %ld\n", nseg2);
  printf("nseg3    %ld\n", nseg3);
  printf("nseg4    %ld\n", nseg4);
  printf("Area     %lf\n", Area);
  printf("Rect     %lf\n", rectangleArea);
  printf("Ltot     %lf\n", amax + xBoxMax);
  printf("ecc      %lf\n", ecc);
  printf("amax     %lf\n", amax);
  printf("bmin     %lf\n", bmin);
  printf("H        %lf\n", H);
  printf("kappa    %lf\n", kappa);
  printf("drmin    %lf\n", drmin);
  printf("drmax    %lf\n", drmax);
  printf("gridspace %lf\n", gridspace);
  printf("gridspacex_real %lf\n", gridspacex_real);
  printf("gridspacey_real %lf\n", gridspacey_real);
  printf("ngridx  %ld\n", ngridx);
  printf("ngridy  %ld\n", ngridy);

  printf("\n");

  printf("ncyc     %ld\n", ncyc);
  printf("neq      %ld\n", neq);
  printf("rmax     %lf\n", rmax);
  printf("delphi_max   %lf\n", delphi_max);
  printf("rshift_max   %lf\n", rshift_max);
  printf("iseed    %ld\n", iseed);
  printf("\n");

  printf("freq_samp  %ld\n", freq_samp);
  printf("cmFreqSamp  %ld\n", cmFreqSamp);
  printf("\n");

  printf("imov     %ld\n", imov);
  printf("plasRigid    %ld\n", plasRigid);
  printf("xcmPrint     %ld\n", xcmPrint);
  printf("ycmPrint     %ld\n", ycmPrint);
  printf("\n");
}

// ----------------------------------------------------------------------
// stateCheckSingle (Double, Triple etc) documents the number of times
// for every measurement cycle the system enters a particular microstate.
//
//
//
// ----------------------------------------------------------------------

void stateCheckSingle(double cmx)
{
  if (cmx < -xBoxMaxd2 || cmx > xBoxMaxd2)
  {
    state[0]++;
  }
  else 
  {
    state[1]++;
  }
}

void stateCheckDouble(double cmx2, double cmx3)
{
  bool cmx2E1 = cmx2<-xBoxMaxd2, cmx2E2 = cmx2> xBoxMaxd2, cmx2B = cmx2 > -xBoxMaxd2 && cmx2 < xBoxMaxd2;
  bool cmx3E1 = cmx3<-xBoxMaxd2, cmx3E2 = cmx3> xBoxMaxd2, cmx3B = cmx3 > -xBoxMaxd2 && cmx3 < xBoxMaxd2;

  if ((cmx2E1 && cmx2E1) || (cmx2E2 && cmx3E2))
  {
    state[0]++;
  }
  else if ((cmx2E1 && cmx3B) || (cmx3E1 && cmx2B))
  {
    state[1]++;
  }
  else if ((cmx2E2 && cmx3B) || (cmx3E2 && cmx2B))
  {
    state[1]++;
  }
  else if ((cmx3E1 && cmx2E2) || (cmx2E1 && cmx3E2))
  {
    state[2]++;
  }
  else if (cmx2B && cmx3B)
  {
    state[3]++;
  }
  else
  {
    printf("Missing a microstate.");
  }
}

void stateCheckTriple(double cmx2, double cmx3, double cmx4)
{ // "E" for inside of an ellipse 1 is -ve x or left, 2 is +ve x or right
  // "B" for inside the rectangle box area.
  bool cmx2E1 = cmx2<-xBoxMaxd2, cmx2E2 = cmx2> xBoxMaxd2, cmx2B = cmx2 > -xBoxMaxd2 && cmx2 < xBoxMaxd2;
  bool cmx3E1 = cmx3<-xBoxMaxd2, cmx3E2 = cmx3> xBoxMaxd2, cmx3B = cmx3 > -xBoxMaxd2 && cmx3 < xBoxMaxd2;
  bool cmx4E1 = cmx4<-xBoxMaxd2, cmx4E2 = cmx4> xBoxMaxd2, cmx4B = cmx4 > -xBoxMaxd2 && cmx4 < xBoxMaxd2;

  if ((cmx2E1 && cmx3E1 && cmx4E1) || (cmx2E2 && cmx3E2 && cmx4E2))
  {
    state[0]++;
  }
  else if ((cmx2B && cmx3E1 && cmx4E1 || cmx2B && cmx3E2 && cmx4E2) || (cmx3B && cmx2E1 && cmx4E1 || cmx3B && cmx2E2 && cmx4E2) || (cmx4B && cmx3E1 && cmx2E1 || cmx4B && cmx3E2 && cmx2E2))
  {
    state[1]++;
  }
  else if ((cmx2E1 && cmx3E2 && cmx4E2 || cmx2E2 && cmx3E1 && cmx4E1) || (cmx3E1 && cmx2E2 && cmx4E2 || cmx3E2 && cmx2E1 && cmx4E1) || (cmx4E1 && cmx3E2 && cmx2E2 || cmx4E2 && cmx3E1 && cmx2E1))
  {
    state[2]++;
  }
  else if ((cmx2E1 && cmx3B && cmx4B || cmx2E2 && cmx3B && cmx4B) || (cmx3E1 && cmx2B && cmx4B || cmx3E2 && cmx2B && cmx4B) || (cmx4E1 && cmx3B && cmx2B || cmx4E2 && cmx3B && cmx2B))
  {
    state[3]++;
  }
  else if ((cmx2E1 && cmx3B && cmx4E2 || cmx2E2 && cmx3B && cmx4E1) || (cmx3E1 && cmx2B && cmx4E2 || cmx3E2 && cmx2B && cmx4E1) || (cmx3E1 && cmx4B && cmx2E2 || cmx3E2 && cmx4B && cmx2E1))
  {
    state[4]++;
  }
  else if (cmx2B && cmx3B && cmx4B)
  {
    state[5]++;
  }
  else
  {
    printf("Missing a microstate.");
  }
}

/* squareEllipse and checkEllipse Code written by Zach Evans to create geometry of rectangle between two halves of an ellipse */

/*
 *
 * checkEllipse is the function that, if it is determined the check is required to take into account the ellipse, checks whether or not
 * the polymer is overlapping with the geometry.
 *
 */

int checkEllipse(double xPos, double yPos, double zPos)
{
  int reject = 1, accept = 0;
  double echeck;

  if (zPos > Hd2 || zPos < -Hd2) // If z position is greater than flat surface of container, reject move
  {
    return reject;
  }

  if (ecc >= 1.0)
  {
    if (xPos > xBoxMaxd2 || xPos < -xBoxMaxd2)
    {
      return (reject);
    }
  }

  if (xPos > xBoxMaxd2 && ecc < 1.0)
  { // If the polymer is outside of the leftmost semi-ellipse, reject
    if ((xPos - xBoxMaxd2) * (xPos - xBoxMaxd2) > amax2 * (1 - (yPos * yPos) / bmin2) && rectangleArea > 0.0)
    {
      return reject;
    }

    if (yPos * yPos > bmin2 * (1 - (xPos - xBoxMaxd2) * (xPos - xBoxMaxd2) / amax2) && rectangleArea > 0.0)
    {
      return reject;
    }

    if (rectangleArea > 0.0)
    {
      echeck = (((xPos - xBoxMaxd2) * (xPos - xBoxMaxd2)) / amax2) + ((yPos * yPos) / bmin2);
    }
    else
    {
      echeck = (((xPos) * (xPos)) / amax2) + ((yPos * yPos) / bmin2);
    }

    if (echeck > 1.0)
    {
      return (reject);
    }
  }

  else if (xPos < -xBoxMaxd2 && ecc < 1.0)
  { // Checking if outside of left elliptical end
    if ((xPos + xBoxMaxd2) * (xPos + xBoxMaxd2) > amax2 * (1 - (yPos * yPos) / bmin2) && rectangleArea > 0.0)
    {
      return reject;
    }

    if (yPos * yPos > bmin2 * (1 - (xPos + xBoxMaxd2) * (xPos + xBoxMaxd2) / amax2) && rectangleArea > 0.0)
    {
      return reject;
    }

    if (rectangleArea > 0.0)
    {
      echeck = (((xPos + xBoxMaxd2) * (xPos + xBoxMaxd2)) / amax2) + ((yPos * yPos) / bmin2);
    }
    else
    {
      echeck = (((xPos) * (xPos)) / amax2) + ((yPos * yPos) / bmin2);
    }

    if (echeck > 1.0)
    {
      return reject;
    }
  }

  return accept; // If inside container, accept move
}

/*
 * squareEllipse takes the position vectors of each monomer in the polymer and then determines:
 * A) Whether or not we are required to take into account the rectangular or the ellipse geometry of the system
 * B) Whether or not the move is accepted or rejected based on the move that was just made.
 */
int squareEllipse(double xPos, double yPos, double zPos)
{
  int reject = 1, accept = 0;

  if (zPos < -Hd2 || zPos > Hd2) // Check if outside of the flat surface of the container
  {
    return (reject);
  }

  if (xPos < -xBoxMaxd2 || xPos > xBoxMaxd2) // If monomer is outside of the rectangle, check if outside of ellipse
  {
    return checkEllipse(xPos, yPos, zPos);
  }

  if (yPos > yBoxMaxd2 || yPos < -yBoxMaxd2)
  {
    return (reject);
  }

  return accept;
}

// ----------------------------------------------------------------------
// This function determines whether the move just performed in the main
// function is permitted or not (checks overlap with a number of conditions)
// ----------------------------------------------------------------------
int check_accept(double rx[5000], double ry[5000], double rz[5000], long nseg)
{
  int accept, reject;
  long klow, khigh;

  accept = 0;
  reject = 1;

  if (ichain == 1)
  {
    // Checking if the T4 polymer overlaps with itself
    for (kk = 0; kk < nseg1 + nseg2 + nseg3 + nseg4; kk++)
    {
      // Check to see if iterative constant is greater than the size of all
      // polymers, if so, break inner loop and continue to next monomer in
      // checked polymer.
      if (kk > nseg1 && kk > nseg2 && kk > nseg3 && kk > nseg4)
      {
        break;
      }

      if (kk < nseg1)
      {
        if (squareEllipse(rx[kk], ry[kk], rz[kk]) == reject)
        {
          return (reject);
        }
        if (kk < k - 1 || kk > k + 1)
        {
          dx = rx[k] - rx[kk];
          dy = ry[k] - ry[kk];
          dz = rz[k] - rz[kk];
          dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 < 1.0)
          {
            return (reject);
          }
        }
      }

      // Checking if the plasmid overlaps with the T4 polymer
      if (kk < nseg2)
      {
        dx = rx[k] - r2x[kk];
        dy = ry[k] - r2y[kk];
        dz = rz[k] - r2z[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
        {
          return (reject);
        }
      }

      // Checking if third plasmid overlaps with T4 polymer
      if (kk < nseg3)
      {
        dx = rx[k] - r3x[kk];
        dy = ry[k] - r3y[kk];
        dz = rz[k] - r3z[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
        {
          return (reject);
        }
      }

      if (kk < nseg4)
      {
        dx = rx[k] - r4x[kk];
        dy = ry[k] - r4y[kk];
        dz = rz[k] - r4z[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
        {
          return (reject);
        }
      }
    }

    return (check_poly_energy(rx, ry, rz, nseg)); // apply rigidity
  }

  else
  {
    if (k == 0)
    {
      klow = nseg - 1;
      khigh = 1;
    }
    else if (k == nseg - 1)
    {
      klow = nseg - 2;
      khigh = 0;
    }
    else
    {
      klow = k - 1;
      khigh = k + 1;
    }

    for (kk = 0; kk < nseg1 + nseg2 + nseg3 + nseg4; kk++)
    {
      // Check to see if iterative constant is greater than the size of all
      // polymers, if so, break inner loop and continue to next monomer in
      // checked polymer.
      if (kk > nseg1 && kk > nseg2 && kk > nseg3 && kk > nseg4)
      {
        break;
      }

      // Check if nseg=2 plasmid escapes squareEllipse or overlaps
      if (kk < nseg)
      {
        if (squareEllipse(rx[kk], ry[kk], rz[kk]) == reject)
        {
          return (reject);
        }
        if (kk != k && kk != klow && kk != khigh)
        {
          dx = rx[k] - rx[kk];
          dy = ry[k] - ry[kk];
          dz = rz[k] - rz[kk];
          dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 < 1.0)
          {
            return (reject);
          }
        }
      }

      if (kk < nseg1 && ichain != 1)
      {
        // Check if polymer and plasmid overlap
        dx = rx[k] - r1x[kk];
        dy = ry[k] - r1y[kk];
        dz = rz[k] - r1z[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
        {
          return (reject);
        }
      }

      // Check if plasmids overlap
      if (kk < nseg2 && ichain != 2)
      {
        dx = rx[k] - r2x[kk];
        dy = ry[k] - r2y[kk];
        dz = rz[k] - r2z[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
        {
          return (reject);
        }
      }

      // Check if plasmids overlap
      if (kk < nseg3 && ichain != 3)
      {
        dx = rx[k] - r3x[kk];
        dy = ry[k] - r3y[kk];
        dz = rz[k] - r3z[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
        {
          return (reject);
        }
      }

      if (kk < nseg4 && ichain != 4)
      {
        dx = rx[k] - r4x[kk];
        dy = ry[k] - r4y[kk];
        dz = rz[k] - r4z[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
        {
          return (reject);
        }
      }
    }

    if (plasRigid == 1)
    {
      return (check_plasmid_energy(rx, ry, rz, nseg)); // apply rigidity
    }
    else
    {
      return (accept);
    }
  }
}

// ----------------------------------------------------------------------
// This function is called to determine whether or not a move is accepted
// based on considerations of polymer chain energy related to the
// move. See comments within the function for details.
// ----------------------------------------------------------------------
int check_poly_energy(double rx[5000], double ry[5000], double rz[5000], long nseg)
{
  int accept, reject; // will return either accept or reject at end of function

  accept = 0;
  reject = 1;

  // reset energy for all of, up to three angles being considered
  for (ind = 0; ind < 3; ind++)
  {
    energy_new[ind] = 0.0;
    energy_old[ind] = 0.0;
  }

  // reset total energies:
  E_new = 0.0;
  E_old = 0.0;

  // The following 5 blocks consider the possible scenarios of monomer
  // movement: 1.) first monomer, 2.) last monomer, 3.) second monomer,
  // 4.) second last monomer, and 5.) any other monomer
  //
  // shifting monomers at the ends requires only one change in angle
  // and thus energy. Monomers second from the end require two, and
  // all the rest require three.
  if (k == 0)
  {
    theta_new = calc_cosine(k, k + 1, k + 2, rx, ry, rz);
    theta_old = calc_cosine(-1, k + 1, k + 2, rx, ry, rz);
    energy_new[0] = kappa * (1.0 - theta_new);
    energy_old[0] = kappa * (1.0 - theta_old);
  }
  else if (k == nseg - 1)
  {
    theta_new = calc_cosine(k - 2, k - 1, k, rx, ry, rz);
    theta_old = calc_cosine(k - 2, k - 1, -1, rx, ry, rz);
    energy_new[0] = kappa * (1.0 - theta_new);
    energy_old[0] = kappa * (1.0 - theta_old);
  }
  else if (k == 1)
  {
    theta_new = calc_cosine(k - 1, k, k + 1, rx, ry, rz);
    theta_old = calc_cosine(k - 1, -1, k + 1, rx, ry, rz);
    energy_new[0] = kappa * (1.0 - theta_new);
    energy_old[0] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(k, k + 1, k + 2, rx, ry, rz);
    theta_old = calc_cosine(-1, k + 1, k + 2, rx, ry, rz);
    energy_new[1] = kappa * (1.0 - theta_new);
    energy_old[1] = kappa * (1.0 - theta_old);
  }
  else if (k == nseg - 2)
  {
    theta_new = calc_cosine(k - 2, k - 1, k, rx, ry, rz);
    theta_old = calc_cosine(k - 2, k - 1, -1, rx, ry, rz);
    energy_new[0] = kappa * (1.0 - theta_new);
    energy_old[0] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(k - 1, k, k + 1, rx, ry, rz);
    theta_old = calc_cosine(k - 1, -1, k + 1, rx, ry, rz);
    energy_new[1] = kappa * (1.0 - theta_new);
    energy_old[1] = kappa * (1.0 - theta_old);
  }
  else
  {
    theta_new = calc_cosine(k - 2, k - 1, k, rx, ry, rz);
    theta_old = calc_cosine(k - 2, k - 1, -1, rx, ry, rz);
    energy_new[0] = kappa * (1.0 - theta_new);
    energy_old[0] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(k - 1, k, k + 1, rx, ry, rz);
    theta_old = calc_cosine(k - 1, -1, k + 1, rx, ry, rz);
    energy_new[1] = kappa * (1.0 - theta_new);
    energy_old[1] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(k, k + 1, k + 2, rx, ry, rz);
    theta_old = calc_cosine(-1, k + 1, k + 2, rx, ry, rz);
    energy_new[2] = kappa * (1.0 - theta_new);
    energy_old[2] = kappa * (1.0 - theta_old);
  }

  // Get the total new and and old energies of the bond angles affected
  // by the move
  for (ind = 0; ind < 3; ind++)
  {
    E_new += energy_new[ind];
    E_old += energy_old[ind];
  }

  delta_E = E_new - E_old;

  // If change in energy is negative, accept the move. If positive, must do a
  // comparison of the change in energy (multiplied by Boltzmann factor and placed
  // in an exponential, with a random value [0,1] to determine if move is
  // accepted or rejected. (Higher delta_E leads to reduced chance of accepted move)
  if (delta_E <= 0)
    return (accept);
  else
  {
    if (ran3() <= exp(-1.0 * delta_E))
    {
      return (accept);
    }
    else
    {
      return (reject);
    }
  }
}

int check_plasmid_energy(double rx[5000], double ry[5000], double rz[5000], long nseg)
{
  int accept, reject; // will return either accept or reject at end of function

  accept = 0;
  reject = 1;

  // reset energy for all of, up to three angles being considered
  for (ind = 0; ind < 3; ind++)
  {
    energy_new[ind] = 0.0;
    energy_old[ind] = 0.0;
  }

  // reset total energies:
  E_new = 0.0;
  E_old = 0.0;

  // The following 3 blocks consider the possible scenarios of monomer
  // movement: 1.) first monomer ( k ==0 ), 2.) last monomer ( k == nseg* -1 )
  // 3.) the second monomer ( k==1 ) 4.) The second last monomer ( k == nseg* -2 )
  // shifting monomers at the ends requires a change in the indices. As k-1 for the k == 0 monomer
  // will give a value of -1, or the -1 indices of the polymer position array, which will provide incorrect
  // results.
  // Similar changes must be made for k values of 1, and nseg*-2

  if (k == 0)
  {
    theta_new = calc_cosine(nseg - 2, nseg - 1, k, rx, ry, rz);
    theta_old = calc_cosine(nseg - 2, nseg - 1, -1, rx, ry, rz);
    energy_new[0] = kappa * (1.0 - theta_new);
    energy_old[0] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(nseg - 1, k, k + 1, rx, ry, rz);
    theta_old = calc_cosine(nseg - 1, -1, k + 1, rx, ry, rz);
    energy_new[1] = kappa * (1.0 - theta_new);
    energy_old[1] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(k, k + 1, k + 2, rx, ry, rz);
    theta_old = calc_cosine(-1, k + 1, k + 2, rx, ry, rz);
    energy_new[2] = kappa * (1.0 - theta_new);
    energy_old[2] = kappa * (1.0 - theta_old);
  }
  else if (k == nseg - 1)
  {
    theta_new = calc_cosine(k - 2, k - 1, k, rx, ry, rz);
    theta_old = calc_cosine(k - 2, k - 1, -1, rx, ry, rz);
    energy_new[0] = kappa * (1.0 - theta_new);
    energy_old[0] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(k - 1, k, 0, rx, ry, rz);
    theta_old = calc_cosine(k - 1, -1, 0, rx, ry, rz);
    energy_new[1] = kappa * (1.0 - theta_new);
    energy_old[1] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(k, 0, 1, rx, ry, rz);
    theta_old = calc_cosine(-1, 0, 1, rx, ry, rz);
    energy_new[2] = kappa * (1.0 - theta_new);
    energy_old[2] = kappa * (1.0 - theta_old);
  }
  else if (k == 1)
  {
    theta_new = calc_cosine(nseg - 1, k - 1, k, rx, ry, rz);
    theta_old = calc_cosine(nseg - 1, k - 1, -1, rx, ry, rz);
    energy_new[0] = kappa * (1.0 - theta_new);
    energy_old[0] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(k - 1, k, k + 1, rx, ry, rz);
    theta_old = calc_cosine(k - 1, -1, k + 1, rx, ry, rz);
    energy_new[1] = kappa * (1.0 - theta_new);
    energy_old[1] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(k, k + 1, k + 2, rx, ry, rz);
    theta_old = calc_cosine(-1, k + 1, k + 2, rx, ry, rz);
    energy_new[2] = kappa * (1.0 - theta_new);
    energy_old[2] = kappa * (1.0 - theta_old);
  }
  else if (k == nseg - 2)
  {
    theta_new = calc_cosine(k - 2, k - 1, k, rx, ry, rz);
    theta_old = calc_cosine(k - 2, k - 1, -1, rx, ry, rz);
    energy_new[0] = kappa * (1.0 - theta_new);
    energy_old[0] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(k - 1, k, k + 1, rx, ry, rz);
    theta_old = calc_cosine(k - 1, -1, k + 1, rx, ry, rz);
    energy_new[1] = kappa * (1.0 - theta_new);
    energy_old[1] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(k, k + 1, 0, rx, ry, rz);
    theta_old = calc_cosine(-1, k + 1, 0, rx, ry, rz);
    energy_new[2] = kappa * (1.0 - theta_new);
    energy_old[2] = kappa * (1.0 - theta_old);
  }
  else
  {
    theta_new = calc_cosine(k - 2, k - 1, k, rx, ry, rz);
    theta_old = calc_cosine(k - 2, k - 1, -1, rx, ry, rz);
    energy_new[0] = kappa * (1.0 - theta_new);
    energy_old[0] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(k - 1, k, k + 1, rx, ry, rz);
    theta_old = calc_cosine(k - 1, -1, k + 1, rx, ry, rz);
    energy_new[1] = kappa * (1.0 - theta_new);
    energy_old[1] = kappa * (1.0 - theta_old);

    theta_new = calc_cosine(k, k + 1, k + 2, rx, ry, rz);
    theta_old = calc_cosine(-1, k + 1, k + 2, rx, ry, rz);
    energy_new[2] = kappa * (1.0 - theta_new);
    energy_old[2] = kappa * (1.0 - theta_old);
  }

  // Get the total new and and old energies of the bond angles affected
  // by the move
  for (ind = 0; ind < 3; ind++)
  {
    E_new += energy_new[ind];
    E_old += energy_old[ind];
  }

  delta_E = E_new - E_old;

  // If change in energy is negative, accept the move. If positive, must do a
  // comparison of the change in energy (multiplied by Boltzmann factor and placed
  // in an exponential, with a random value [0,1] to determine if move is
  // accepted or rejected. (Higher delta_E leads to reduced chance of accepted move)
  if (delta_E <= 0)
    return (accept);
  else
  {
    if (ran3() <= exp(-1.0 * delta_E))
    {
      return (accept);
    }
    else
    {
      return (reject);
    }
  }
}

/*
 * Generalized calc_cosine method in order to determine the angle between the links in the plasmids
 * and the linear polymer for the purposes of rigidity or bending resistance of the polymers.
 *
 * Written by Zach Evans 2023 June 30
 * zdjevans@protonmail.com
 */

// ----------------------------------------------------------------------
// Calculate the cosine of theta between two bonds (vectors) that connect
// three sequential monomers within the first polymer chain.
// ----------------------------------------------------------------------
double calc_cosine(int i1, int i2, int i3, double rx[5000], double ry[5000], double rz[5000])
{
  // determine which of the indices is negative (if any). The negative index
  // corresponds to the old position of the monomer. Set x, y, and z values of
  // each of the three monomers accordingly. Note y1 is an implicit variable
  // name in some stlib function, so use yone instead.

  if (i1 < 0)
  {
    x1 = xold;
    yone = yold, z1 = zold;
    x2 = rx[i2];
    y2 = ry[i2];
    z2 = rz[i2];
    x3 = rx[i3];
    y3 = ry[i3];
    z3 = rz[i3];
  }
  else if (i2 < 0)
  {
    x1 = rx[i1];
    yone = ry[i1];
    z1 = rz[i1];
    x2 = xold;
    y2 = yold;
    z2 = zold;
    x3 = rx[i3];
    y3 = ry[i3];
    z3 = rz[i3];
  }
  else if (i3 < 0)
  {
    x1 = rx[i1];
    yone = ry[i1];
    z1 = rz[i1];
    x2 = rx[i2];
    y2 = ry[i2];
    z2 = rz[i2];
    x3 = xold;
    y3 = yold;
    z3 = zold;
  }
  else
  {
    x1 = rx[i1];
    yone = ry[i1];
    z1 = rz[i1];
    x2 = rx[i2];
    y2 = ry[i2];
    z2 = rz[i2];
    x3 = rx[i3];
    y3 = ry[i3];
    z3 = rz[i3];
  }

  vax = x2 - x1;
  vay = y2 - yone;
  vaz = z2 - z1;
  vbx = x3 - x2;
  vby = y3 - y2;
  vbz = z3 - z2;

  va_sq = vax * vax + vay * vay + vaz * vaz;
  vb_sq = vbx * vbx + vby * vby + vbz * vbz;

  va_dot_vb = vax * vbx + vay * vby + vaz * vbz;

  return (va_dot_vb / (sqrt(va_sq * vb_sq)));
}

// ----------------------------------------------------------------------
//  Places polymer in original (unrealistic position); is called for each
//  unique window. Overlaps the polymers somewhere to match the current
//  window.
// ----------------------------------------------------------------------

void init_pos(void)
{

  double xadd, yadd, xmax, ymax, zplace, echeck;

  r1x[0] = -0.0;
  r1y[0] = -bmin + 4.0;
  r1z[0] = 1.0;
  xadd = 1.0;
  yadd = 1.0;
  zplace = 1.0;

  xmax = xBoxMaxd2;
  ymax = yBoxMaxd2;

  for (i = 1; i < nseg1; i++)
  {
    r1x[i] = r1x[i - 1] + xadd;
    r1y[i] = r1y[i - 1];
    r1z[i] = zplace;
    if (rectangleArea > 0.0)
    {
      echeck = (((r1x[i] - xBoxMaxd2) * (r1x[i] - xBoxMaxd2)) / amax2) + ((r1y[i] * r1y[i]) / bmin2); // Changed to be inside of the rectangle structure
      if (r1x[i] > xBoxMaxd2 || r1x[i] < -xBoxMaxd2)
      {
        if (echeck > 0.8)
        {
          if (r1y[i] + 1.0 > ymax || r1y[i] - 1.0 < -ymax)
          {
            zplace -= 1;
          }
          r1x[i] -= xadd;
          r1y[i] = r1y[i] + yadd;
          xadd *= -1.0;
        }
      }
      else
      {
        continue;
      }
    }

    else
    {
      echeck = ((r1x[i] * r1x[i]) / amax2) + ((r1y[i] * r1y[i]) / bmin2);

      if (echeck > 0.8)
      {
        if (r1y[i] + 1.0 > ymax || r1y[i] - 1.0 < -ymax)
        {
          zplace -= 1;
        }
        r1x[i] -= xadd;
        r1y[i] = r1y[i] + yadd;
        xadd *= -1.0;
      }
    }
  }

  for (i = 0; i < nseg2; i++)
  {
    r2z[i] = 2.0;

    if (i < nseg2 / 2)
    {
      r2x[i] = (double)-i;
      r2y[i] = 0.0;
    }
    if (i == nseg2 / 2 + 0.5 || i == nseg2 / 2)
    {
      r2x[i] = (double)-i;
      r2y[i] = -1.0;
    }
    if (i > nseg2 / 2 && i < nseg2 - 1)
    {
      r2x[i] = (double)+i - nseg2;
      r2y[i] = -2.0;
    }
    if (i == nseg2 - 1)
    {
      r2x[i] = (double)+i - nseg2;
      r2y[i] = -1.0;
    }
  }

  for (i = 0; i < nseg3; i++)
  {
    r3z[i] = 4.0;

    if (i < nseg3 / 2)
    {
      r3x[i] = (double)-i;
      r3y[i] = -0.0;
    }
    if (i == nseg3 / 2 + 0.5 || i == nseg3 / 2)
    {
      r3x[i] = (double)-i;
      r3y[i] = -1.0;
    }
    if (i > nseg3 / 2 && i < nseg3 - 1)
    {
      r3x[i] = (double)+i - nseg3;
      r3y[i] = -2.0;
    }
    if (i == nseg3 - 1)
    {
      r3x[i] = (double)+i - nseg3;
      r3y[i] = -1.0;
    }
  }

  for (i = 0; i < nseg4; i++)
  {
    r4z[i] = 6.0;

    if (i < nseg4 / 2)
    {
      r4x[i] = (double)-i;
      r4y[i] = 0.0;
    }
    if (i == nseg4 / 2 + 0.5 || i == nseg4 / 2)
    {
      r4x[i] = (double)-i;
      r4y[i] = -1.0;
    }
    if (i > nseg4 / 2 && i < nseg4 - 1)
    {
      r4x[i] = (double)+i - nseg4;
      r4y[i] = -2.0;
    }
    if (i == nseg4 - 1)
    {
      r4x[i] = (double)+i - nseg4;
      r4y[i] = -1.0;
    }
  }
}

// Initialize plasmids as circles rather than trapezoid
void init_pos_circular(void)
{
  double xadd, yadd, xmax, ymax, zplace;

  r1x[0] = -xBoxMaxd2;
  r1y[0] = -yBoxMaxd2 + 2.0;
  r1z[0] = 1.0;
  xadd = 1.0;
  yadd = 1.0;
  zplace = 1.0;

  for (i = 1; i < nseg1; i++)
  {
    r1x[i] = r1x[i - 1] + xadd;
    r1y[i] = r1y[i - 1];
    r1z[i] = zplace;
    xmax = xBoxMaxd2; // Changed to be inside of the rectangle structure
    // xmax = amax * sqrt(1.0 - pow(r1y[i] / bmin, 2.0)) - 3.0;

    if (r1x[i] > xmax || r1x[i] < -xmax)
    {
      if (r1y[i] + 2.0 > ymax || r1y[i] - 2.0 < -ymax)
      {
        zplace -= 1.0;
        yadd = -1.0 * yadd;
        // printf("%ld   %lf\n", i, yadd);
      }

      r1x[i] -= xadd;
      r1y[i] = r1y[i] + yadd;
      ymax = yBoxMaxd2;

      if (r1y[i] > ymax || r1y[i] < -ymax)
      {
        printf("Can't place polymer... exiting...\n");
        exit(0);
      }
      xadd *= -1.0;
    }
  }

  double theta_plasmid2 = 2.0 * PI / nseg2;
  double Rplasmid2 = 0.5 / tan(theta_plasmid2 / 2.0);

  double theta_plasmid3 = 2.0 * PI / nseg3;
  double Rplasmid3 = 0.5 / tan(theta_plasmid3 / 2.0);

  double theta_plasmid4 = 2.0 * PI / nseg4;
  double Rplasmid4 = 0.5 / tan(theta_plasmid4 / 2.0);

  for (i = 0; i < nseg2; i++)
  {
    r2z[i] = 2.0;
    r2x[i] = Rplasmid2 * cos(i * theta_plasmid2) - xBoxMaxd2 + Rplasmid2;
    r2y[i] = Rplasmid2 * sin(i * theta_plasmid2) - Rplasmid2 / 2.0;
  }

  for (i = 0; i < nseg3; i++)
  {
    r3z[i] = 4.0; // Initialized just above the first plasmid
    r3x[i] = Rplasmid3 * cos(i * theta_plasmid3) - xBoxMaxd2 + Rplasmid3;
    r3y[i] = Rplasmid3 * sin(i * theta_plasmid3) - Rplasmid3 / 2.0;
  }

  for (i = 0; i < nseg4; i++)
  {
    r4z[i] = 6.0; // Initialized just above the first plasmid
    r4x[i] = Rplasmid4 * cos(i * theta_plasmid4) - xBoxMaxd2 + Rplasmid4;
    r4y[i] = Rplasmid4 * sin(i * theta_plasmid4) - Rplasmid4 / 2.0;
  }
}

void write_data(void)
{
  FILE *fp;

  if ((fp = fopen("plas.dat", "w")) == NULL)
  {
    printf("Cannot open file: plas.dat\n");
    exit(0);
  }
  else
  {
    for (i = 0; i < ngridx; i++)
    {
      for (j = 0; j < ngridy; j++)
      {
        fprintf(fp, "%8.2lf  ", plas[i][j]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  if ((fp = fopen("areaProb.dat", "w")) == NULL)
  {
    printf("Cannot open file: areaProb.dat\n");
    exit(0);
  }
  else
  {
    int number = 6;

    for (j = 0; j < number; j++)
    {
      if (state[j] != 0)
      {
        fprintf(fp, "%ld\t%ld\n", j + 1, state[j]);
      }
    }
    fclose(fp);
  }

  if (nseg1 != 0)
  {
    if ((fp = fopen("prob1.dat", "w")) == NULL)
    {
      printf("Cannot open file: prob1.dat\n");
      exit(0);
    }
    else
    {
      for (i = 0; i < ngridx; i++)
      {
        for (j = 0; j < ngridy; j++)
        {
          fprintf(fp, "%8.2lf  ", prob1[i][j]);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }

  if (nseg1 != 0)
  {
    if ((fp = fopen("probmon.dat", "w")) == NULL)
    {
      printf("Cannot open file: probmon.dat\n");
      exit(0);
    }
    else
    {
      for (i = 0; i < ngridx; i++)
      {
        for (j = 0; j < ngridy; j++)
        {
          fprintf(fp, "%8.2lf  ", probmon[i][j]);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }

  if (nseg2 != 0 && nseg3 != 0)
  {
    if ((fp = fopen("x2x3cm.dat", "w")) == NULL)
    {
      printf("Cannot open file: x2x3cm.dat\n");
      exit(0);
    }
    else
    {
      for (j = neq / cmFreqSamp + 1; j < iter; j++)
      {
        fprintf(fp, "%8.2lf  ", plas12[j]);
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }

  if (nseg2 != 0 && nseg4 != 0)
  {
    if ((fp = fopen("x2x4cm.dat", "w")) == NULL)
    {
      printf("Cannot open file: x2x4cm.dat\n");
      exit(0);
    }
    else
    {
      for (j = neq / cmFreqSamp + 1; j < iter; j++)
      {
        fprintf(fp, "%8.2lf  ", plas13[j]);
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }

  if (nseg3 != 0 && nseg4 != 0)
  {
    if ((fp = fopen("x3x4cm.dat", "w")) == NULL)
    {
      printf("Cannot open file: x3x4cm.dat\n");
      exit(0);
    }
    else
    {
      for (j = neq / cmFreqSamp + 1; j < iter; j++)
      {
        fprintf(fp, "%8.2lf  ", plas23[j]);
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }

  if (nseg2 != 0)
  {
    if (xcmPrint == 1)
      if ((fp = fopen("x2cm.dat", "w")) == NULL)
      {
        printf("Cannot open file: x2cm.dat\n");
        exit(0);
      }
      else
      {
        for (j = neq / cmFreqSamp + 1; j < iter; j++)
        {
          fprintf(fp, "%8.2lf  ", x2cm[j]);
          fprintf(fp, "\n");
        }
        fclose(fp);
      }

    if (ycmPrint == 1)
      if ((fp = fopen("y2cm.dat", "w")) == NULL)
      {
        printf("Cannot open file: y2cm.dat\n");
        exit(0);
      }
      else
      {
        for (j = neq / cmFreqSamp + 1; j < iter; j++)
        {
          fprintf(fp, "%8.2lf  ", y2cm[j]);
          fprintf(fp, "\n");
        }
        fclose(fp);
      }
  }

  if (nseg3 != 0)
  {
    if (xcmPrint == 1)
      if ((fp = fopen("x3cm.dat", "w")) == NULL)
      {
        printf("Cannot open file: x3cm.dat\n");
        exit(0);
      }
      else
      {
        for (j = neq / cmFreqSamp + 1; j < iter; j++)
        {
          fprintf(fp, "%8.2lf  ", x3cm[j]);
          fprintf(fp, "\n");
        }
        fclose(fp);
      }

    if (ycmPrint == 1)
      if ((fp = fopen("y3cm.dat", "w")) == NULL)
      {
        printf("Cannot open file: y3cm.dat\n");
        exit(0);
      }
      else
      {
        for (j = neq / cmFreqSamp + 1; j < iter; j++)
        {
          fprintf(fp, "%8.2lf  ", y3cm[j]);
          fprintf(fp, "\n");
        }
        fclose(fp);
      }
  }

  if (nseg4 != 0)
  {
    if (xcmPrint == 1)
      if ((fp = fopen("x4cm.dat", "w")) == NULL)
      {
        printf("Cannot open file: x4cm.dat\n");
        exit(0);
      }
      else
      {
        for (j = neq / cmFreqSamp + 1; j < iter; j++)
        {

          fprintf(fp, "%8.2lf  ", x4cm[j]);
          fprintf(fp, "\n");
        }
        fclose(fp);
      }

    if (ycmPrint == 1)
      if ((fp = fopen("y4cm.dat", "w")) == NULL)
      {
        printf("Cannot open file: y4cm.dat\n");
        exit(0);
      }
      else
      {
        for (j = neq / cmFreqSamp + 1; j < iter; j++)
        {
          fprintf(fp, "%8.2lf  ", y4cm[j]);
          fprintf(fp, "\n");
        }
        fclose(fp);
      }
  }
}

double ran3()
{
  static int ma[60], mj, mk, mz, i, ii, k, inext, inextp, iff, mbig, mseed;

  if (iff == 0)
  {
    iff = 1;
    mbig = 1000000000;
    // mseed = 161803398;
    mseed = iseed; // iseed given in input file: simulation always proceeds same way?
    mz = 0;
    mj = mseed;
    mj = mj % mbig;
    ma[55] = mj;
    mk = 1;
    for (i = 1; i < 55; i++)
    {
      ii = (21 * i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < mz)
        mk = mk + mbig;
      mj = ma[ii];
    }
    for (k = 1; k <= 4; k++)
    {
      for (i = 1; i <= 55; i++)
      {
        ma[i] = ma[i] - ma[1 + ((i + 30) % 55)];
        if (ma[i] < mz)
          ma[i] = ma[i] + mbig;
      }
    }
    inext = 0;
    inextp = 31;
  }
  if (++inext == 56)
    inext = 1;
  if (++inextp == 56)
    inextp = 1;
  mj = ma[inext] - ma[inextp];
  if (mj < mz)
    mj = mj + mbig;
  ma[inext] = mj;
  return (double)mj / mbig;
}

void shift_move_plasmid(double rx[5000], double ry[5000], double rz[5000], long nseg)
{
  double delrx, delry, delrz;
  double xsold[5000];
  double ysold[5000];
  double zsold[5000];

  delrx = rshift_max * (2.0 * ran3() - 1.0);
  delry = rshift_max * (2.0 * ran3() - 1.0);
  delrz = rshift_max * (2.0 * ran3() - 1.0);

  for (i = 0; i < nseg; i++)
  {
    xsold[i] = rx[i];
    ysold[i] = ry[i];
    zsold[i] = rz[i];

    rx[i] += delrx;
    ry[i] += delry;
    rz[i] += delrz;
  }

  // printf("delrx = %lf, delry = %lf, delrz = %lf\n",delrx,delry,delrz);

  overlap = check_shift_chain(rx, ry, rz, nseg);

  nshift += 1;
  if (overlap == 0)
  {
    nacc += 1;
    nacc_shift += 1;
  }
  else if (overlap == 1)
  {
    for (i = 0; i < nseg; i++)
    {
      rx[i] = xsold[i];
      ry[i] = ysold[i];
      rz[i] = zsold[i];
    }
  }
}

int check_shift_chain(double rx[5000], double ry[5000], double rz[5000], long nseg)
{
  int accept = 0;
  int reject = 1;

  for (i = 0; i < nseg; i++)
  {

    if (squareEllipse(rx[i], ry[i], rz[i]) == reject)
    {
      return (reject);
    }

    for (kk = 0; kk < nseg1 + nseg2 + nseg3 + nseg4; kk++)
    {
      // Check to see if iterative constant is greater than the size of all
      // polymers, if so, break inner loop and continue to next monomer in
      // checked polymer.
      if (kk > nseg1 && kk > nseg2 && kk > nseg3 && kk > nseg4)
      {
        break;
      }

      if (kk < nseg1 && ichain != 1)
      {
        dx = rx[i] - r1x[kk];
        dy = ry[i] - r1y[kk];
        dz = rz[i] - r1z[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
        {
          return (reject);
        }
      }

      if (kk < nseg2 && ichain != 2)
      {
        dx = rx[i] - r2x[kk];
        dy = ry[i] - r2y[kk];
        dz = rz[i] - r2z[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
        {
          return (reject);
        }
      }

      if (kk < nseg3 && ichain != 3)
      {
        dx = rx[i] - r3x[kk];
        dy = ry[i] - r3y[kk];
        dz = rz[i] - r3z[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
        {
          return (reject);
        }
      }

      if (kk < nseg4 && ichain != 4)
      {
        dx = rx[i] - r4x[kk];
        dy = ry[i] - r4y[kk];
        dz = rz[i] - r4z[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
        {
          return (reject);
        }
      }
    }
  }
  return (accept);
}

void reptation_move_chain1()
{
  double rannum;
  rannum = ran3();
  if (rannum <= 0.5)
  {
    irep = 0;
  }
  else
  {
    irep = nseg1 - 1;
  }

  rannum = ran3();
  phi_prime = rannum * 2.0 * PI;

  rannum = ran3();
  if (kappa > -0.000001 && kappa < 0.000001)
  {
    costheta_prime = (2.0 * rannum) - 1.0;
  }
  else
  {
    costheta_prime = log((rannum * exp(kappa)) + ((1.0 - rannum) * exp(-1.0 * kappa))) / kappa;
  }

  //
  //      keep bond length = 1
  //
  //      rannum = ran3();
  //      dr_prime = pow((0.602*rannum)+0.729,(1.0/3.0));
  dr_prime = 1.0;

  xold = r1x[irep];
  yold = r1y[irep];
  zold = r1z[irep];

  calc_delta_xyz();

  if (irep == 0)
  {
    for (ind = 0; ind < nseg1 - 1; ind++)
    {
      r1x[ind] = r1x[ind + 1];
      r1y[ind] = r1y[ind + 1];
      r1z[ind] = r1z[ind + 1];
    }

    r1x[nseg1 - 1] = r1x[nseg1 - 2] + dx_fixed;
    r1y[nseg1 - 1] = r1y[nseg1 - 2] + dy_fixed;
    r1z[nseg1 - 1] = r1z[nseg1 - 2] + dz_fixed;

    overlap = check_accept_reptation(r1x, r1y, r1z, nseg1, nseg1 - 1);

    if (overlap == 0)
    {
      nacc_rep += 1;
    }
    else
    {
      for (ind = nseg1 - 1; ind > 0; ind--)
      {
        r1x[ind] = r1x[ind - 1];
        r1y[ind] = r1y[ind - 1];
        r1z[ind] = r1z[ind - 1];
      }

      r1x[0] = xold;
      r1y[0] = yold;
      r1z[0] = zold;
    }
  }
  else
  { // irep == nseg1-1
    for (ind = nseg1 - 1; ind > 0; ind--)
    {
      r1x[ind] = r1x[ind - 1];
      r1y[ind] = r1y[ind - 1];
      r1z[ind] = r1z[ind - 1];
    }

    r1x[0] = r1x[1] + dx_fixed;
    r1y[0] = r1y[1] + dy_fixed;
    r1z[0] = r1z[1] + dz_fixed;

    overlap = check_accept_reptation(r1x, r1y, r1z, nseg1, 0);

    if (overlap == 0)
    {
      nacc_rep += 1;
    }
    else
    {
      for (ind = 0; ind < nseg1 - 1; ind++)
      {
        r1x[ind] = r1x[ind + 1];
        r1y[ind] = r1y[ind + 1];
        r1z[ind] = r1z[ind + 1];
      }

      r1x[nseg1 - 1] = xold;
      r1y[nseg1 - 1] = yold;
      r1z[nseg1 - 1] = zold;
    }
  }
}

void calc_delta_xyz()
{
  dx_prime = dr_prime * sqrt(1.0 - (costheta_prime * costheta_prime)) * cos(phi_prime);
  dy_prime = dr_prime * sqrt(1.0 - (costheta_prime * costheta_prime)) * sin(phi_prime);
  dz_prime = dr_prime * costheta_prime;

  if (ichain == 1)
  {
    if (irep == 0)
    {
      ux = r1x[nseg1 - 1] - r1x[nseg1 - 2];
      uy = r1y[nseg1 - 1] - r1y[nseg1 - 2];
      uz = r1z[nseg1 - 1] - r1z[nseg1 - 2];
    }
    else
    {
      ux = r1x[0] - r1x[1];
      uy = r1y[0] - r1y[1];
      uz = r1z[0] - r1z[1];
    }
  }

  else if (ichain == 2)
  {
    if (irep == 0)
    {
      ux = r2x[nseg2 - 1] - r2x[nseg2 - 2];
      uy = r2y[nseg2 - 1] - r2y[nseg2 - 2];
      uz = r2z[nseg2 - 1] - r2z[nseg2 - 2];
    }
    else
    {
      ux = r2x[0] - r2x[1];
      uy = r2y[0] - r2y[1];
      uz = r2z[0] - r2z[1];
    }
  }

  else if (ichain == 3)
  {
    if (irep == 0)
    {
      ux = r3x[nseg3 - 1] - r3x[nseg3 - 2];
      uy = r3y[nseg3 - 1] - r3y[nseg3 - 2];
      uz = r3z[nseg3 - 1] - r3z[nseg3 - 2];
    }
    else
    {
      ux = r3x[0] - r3x[1];
      uy = r3y[0] - r3y[1];
      uz = r3z[0] - r3z[1];
    }
  }

  else if (ichain == 4)
  {
    if (irep == 0)
    {
      ux = r4x[nseg4 - 1] - r4x[nseg4 - 2];
      uy = r4y[nseg4 - 1] - r4y[nseg4 - 2];
      uz = r4z[nseg4 - 1] - r4z[nseg4 - 2];
    }
    else
    {
      ux = r4x[0] - r4x[1];
      uy = r4y[0] - r4y[1];
      uz = r4z[0] - r4z[1];
    }
  }

  u = sqrt(ux * ux + uy * uy + uz * uz);
  uxy = sqrt(ux * ux + uy * uy);
  cos_beta = uz / u;
  sin_beta = sqrt(1.0 - (cos_beta * cos_beta));

  if (ux > -0.000001 && ux < 0.000001 && uy > -0.000001 && uy < 0.000001)
  {
    cos_alpha = 1.0;
    sin_alpha = 0.0;
  }
  else
  {
    // cos_alpha = (ux/uxy);
    // sin_alpha = sqrt(1.0 -(cos_alpha*cos_alpha));

    if (uy >= 0.0)
    {
      cos_alpha = (ux / uxy);
      sin_alpha = sqrt(1.0 - (cos_alpha * cos_alpha));
    }
    else
    {
      alpha = acos(fabs(ux) / uxy);
      if (ux < 0)
        alpha += PI;
      else
        alpha = (2.0 * PI) - alpha;
      cos_alpha = cos(alpha);
      sin_alpha = sin(alpha);
    }
  }

  // inverted matrix
  dx_fixed = (cos_beta * cos_alpha * dx_prime) - (sin_alpha * dy_prime) + (sin_beta * cos_alpha * dz_prime);
  dy_fixed = (cos_beta * sin_alpha * dx_prime) + (cos_alpha * dy_prime) + (sin_beta * sin_alpha * dz_prime);
  dz_fixed = (cos_beta * dz_prime) - (sin_beta * dx_prime);
}

int check_accept_reptation(double rx[5000], double ry[5000], double rz[5000], long nseg, long krep)
{
  int accept, reject; // will return either accept or reject at end of function

  accept = 0;
  reject = 1;

  if (ichain == 1)
  {
    for (kk = 0; kk < nseg1 + nseg2 + nseg3 + nseg4; kk++)
    {
      // Check to see if iterative constant is greater than the size of all
      // polymers, if so, break inner loop and continue to next monomer in
      // checked polymer.
      if (kk > nseg1 && kk > nseg2 && kk > nseg3 && kk > nseg4)
      {
        break;
      }

      if (squareEllipse(rx[kk], ry[kk], rz[kk]) == reject && kk < nseg)
      {
        return (reject);
      }

      if ((kk < krep - 1 || kk > krep + 1) && kk < nseg)
      {
        dz = rz[krep] - rz[kk];
        if (fabs(dz) < 1.0)
        {
          dx = rx[krep] - rx[kk];
          dy = ry[krep] - ry[kk];
          dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 < 1.0)
          {
            return (reject); // if overlap with other monomer within chain, reject
          }
        }
      }

      if (kk < nseg2)
      {
        dz = rz[krep] - r2z[kk];
        if (fabs(dz) < 1.0)
        {
          dx = rx[krep] - r2x[kk];
          dy = ry[krep] - r2y[kk];
          dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 < 1.0)
            return (reject); // if overlap with monomer in other chain, reject
        }
      }

      if (kk < nseg3)
      {
        dz = rz[krep] - r3z[kk];
        if (fabs(dz) < 1.0)
        {
          dx = rx[krep] - r3x[kk];
          dy = ry[krep] - r3y[kk];
          dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 < 1.0)
          {
            return (reject); // if overlap with monomer in other chain, reject
          }
        }
      }

      if (kk < nseg4)
      {
        dz = rz[krep] - r4z[kk];
        if (fabs(dz) < 1.0)
        {
          dx = rx[krep] - r4x[kk];
          dy = ry[krep] - r4y[kk];
          dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 < 1.0)
          {
            return (reject); // if overlap with monomer in other chain, reject
          }
        }
      }
    }
  }
  else
  { // plasmids

    for (kk = 0; kk < nseg1 + nseg2 + nseg3 + nseg4; kk++)
    {
      // Check to see if iterative constant is greater than the size of all
      // polymers, if so, break inner loop and continue to next monomer in
      // checked polymer.
      if (kk > nseg1 && kk > nseg2 && kk > nseg3 && kk > nseg4)
      {
        break;
      }

      if (kk < nseg)
      {
        if (squareEllipse(rx[kk], ry[kk], rz[kk]) == reject)
        {
          return (reject);
        }
      }

      if (kk < nseg1)
      {
        dz = rz[krep] - r1z[kk];
        if (fabs(dz) < 1.0)
        {
          dx = rx[krep] - r1x[kk];
          dy = ry[krep] - r1y[kk];
          dr2 = dx * dx + dy * dy + dz * dz;

          if (dr2 < 1.0)
          {
            return (reject);
          }
        }
      }

      if (kk < nseg2)
      {
        if ((kk < krep - 1 || kk > krep + 1) && ichain == 2)
        {
          dz = rz[krep] - r2z[kk];
          if (fabs(dz) < 1.0)
          {
            dx = rx[krep] - r2x[kk];
            dy = ry[krep] - r2y[kk];
            dr2 = dx * dx + dy * dy + dz * dz;

            if (dr2 < 1.0)
            {
              return (reject);
            }
          }
        }
        else if (ichain != 2)
        {
          dz = rz[krep] - r2z[kk];
          if (fabs(dz) < 1.0)
          {
            dx = rx[krep] - r2x[kk];
            dy = ry[krep] - r2y[kk];
            dr2 = dx * dx + dy * dy + dz * dz;

            if (dr2 < 1.0)
            {
              return (reject);
            }
          }
        }
      }

      if (kk < nseg3)
      {
        if ((kk < krep - 1 || kk > krep + 1) && ichain == 3)
        {
          dz = rz[krep] - r3z[kk];
          if (fabs(dz) < 1.0)
          {
            dx = rx[krep] - r3x[kk];
            dy = ry[krep] - r3y[kk];
            dr2 = dx * dx + dy * dy + dz * dz;

            if (dr2 < 1.0)
            {
              return (reject);
            }
          }
        }
        else if (ichain != 3)
        {
          dz = rz[krep] - r3z[kk];
          if (fabs(dz) < 1.0)
          {
            dx = rx[krep] - r3x[kk];
            dy = ry[krep] - r3y[kk];
            dr2 = dx * dx + dy * dy + dz * dz;

            if (dr2 < 1.0)
            {
              return (reject);
            }
          }
        }
      }

      if (kk < nseg4)
      {
        if ((kk < krep - 1 || kk > krep + 1) && ichain == 4)
        {
          dz = rz[krep] - r4z[kk];
          if (fabs(dz) < 1.0)
          {
            dx = rx[krep] - r4x[kk];
            dy = ry[krep] - r4y[kk];
            dr2 = dx * dx + dy * dy + dz * dz;

            if (dr2 < 1.0)
            {
              return (reject);
            }
          }
        }
        else if (ichain != 4)
        {
          dz = rz[krep] - r4z[kk];
          if (fabs(dz) < 1.0)
          {
            dx = rx[krep] - r4x[kk];
            dy = ry[krep] - r4y[kk];
            dr2 = dx * dx + dy * dy + dz * dz;

            if (dr2 < 1.0)
            {
              return (reject);
            }
          }
        }
      }
    }
  }

  return (accept);
}

void crank_move_polymer(double rx[5000], double ry[5000], double rz[5000])
{
  double delrx, delry, delrz, Rx, Ry, Rz, Rmag, rdotRn, Rnx, Rny, Rnz;
  double ux, uy, uz, vx, vy, vz, vmag, wx, wy, wz, wmag;
  double cosphi, sinphi, delphi;

  delrx = rx[k] - rx[k - 1];
  delry = ry[k] - ry[k - 1];
  delrz = rz[k] - rz[k - 1];

  Rx = rx[k + 1] - rx[k - 1];
  Ry = ry[k + 1] - ry[k - 1];
  Rz = rz[k + 1] - rz[k - 1];
  Rmag = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);

  Rnx = Rx / Rmag;
  Rny = Ry / Rmag;
  Rnz = Rz / Rmag;

  rdotRn = delrx * Rnx + delry * Rny + delrz * Rnz;
  ux = rdotRn * Rnx;
  uy = rdotRn * Rny;
  uz = rdotRn * Rnz;

  vx = delrx - ux;
  vy = delry - uy;
  vz = delrz - uz;
  vmag = sqrt(vx * vx + vy * vy + vz * vz);
  // if (vmag < 0.00000001) printf("vmag = %lf\n",vmag);

  wx = uy * vz - uz * vy;
  wy = uz * vx - ux * vz;
  wz = ux * vy - uy * vx;
  wmag = sqrt(wx * wx + wy * wy + wz * wz);

  if (wmag > 0.00000001)
  {
    delphi = (2.0 * ran3() - 1.0) * delphi_max;
    cosphi = cos(delphi);
    sinphi = sin(delphi);

    delrx = ux + cosphi * vx + sinphi * vmag * wx / wmag;
    delry = uy + cosphi * vy + sinphi * vmag * wy / wmag;
    delrz = uz + cosphi * vz + sinphi * vmag * wz / wmag;

    rx[k] = rx[k - 1] + delrx;
    ry[k] = ry[k - 1] + delry;
    rz[k] = rz[k - 1] + delrz;
  }
  else
  { // bonds are parallel
    rx[k] = xold;
    ry[k] = yold;
    rz[k] = zold;
  }
}

/*
 * Generalized crank_move method for the plasmid to reduce the amount of code needed to run the program.
 * Uses the crank shaft and checks to see if the ring structure of the plasmid is maintained.
 *
 * Written by Zach Evans 2023 June 30
 * zdjevans@protonmail.com
 */

void crank_move_plasmid(double rx[5000], double ry[5000], double rz[5000], long nseg)
{

  double drx, dry, drz, Rx, Ry, Rz, Rmag, rdotRn, Rnx, Rny, Rnz;
  double ux, uy, uz, vx, vy, vz, vmag, wx, wy, wz, wmag;
  double cosphi, sinphi, delphi;

  if (k > 0 && k < nseg - 1)
  {
    drx = rx[k] - rx[k - 1];
    dry = ry[k] - ry[k - 1];
    drz = rz[k] - rz[k - 1];

    Rx = rx[k + 1] - rx[k - 1];
    Ry = ry[k + 1] - ry[k - 1];
    Rz = rz[k + 1] - rz[k - 1];
    Rmag = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);

    Rnx = Rx / Rmag;
    Rny = Ry / Rmag;
    Rnz = Rz / Rmag;

    rdotRn = drx * Rnx + dry * Rny + drz * Rnz;
    ux = rdotRn * Rnx;
    uy = rdotRn * Rny;
    uz = rdotRn * Rnz;

    vx = drx - ux;
    vy = dry - uy;
    vz = drz - uz;
    vmag = sqrt(vx * vx + vy * vy + vz * vz);
    // if (vmag < 0.00000001) printf("vmag = %lf\n", vmag);

    wx = uy * vz - uz * vy;
    wy = uz * vx - ux * vz;
    wz = ux * vy - uy * vx;
    wmag = sqrt(wx * wx + wy * wy + wz * wz);

    delphi = (2.0 * ran3() - 1.0) * delphi_max;
    cosphi = cos(delphi);
    sinphi = sin(delphi);

    drx = ux + cosphi * vx + sinphi * vmag * wx / wmag;
    dry = uy + cosphi * vy + sinphi * vmag * wy / wmag;
    drz = uz + cosphi * vz + sinphi * vmag * wz / wmag;

    rx[k] = rx[k - 1] + drx;
    ry[k] = ry[k - 1] + dry;
    rz[k] = rz[k - 1] + drz;
  }
  else if (k == 0)
  {

    drx = rx[k] - rx[nseg - 1];
    dry = ry[k] - ry[nseg - 1];
    drz = rz[k] - rz[nseg - 1];

    Rx = rx[k + 1] - rx[nseg - 1];
    Ry = ry[k + 1] - ry[nseg - 1];
    Rz = rz[k + 1] - rz[nseg - 1];
    Rmag = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);

    Rnx = Rx / Rmag;
    Rny = Ry / Rmag;
    Rnz = Rz / Rmag;

    rdotRn = drx * Rnx + dry * Rny + drz * Rnz;
    ux = rdotRn * Rnx;
    uy = rdotRn * Rny;
    uz = rdotRn * Rnz;

    vx = drx - ux;
    vy = dry - uy;
    vz = drz - uz;
    vmag = sqrt(vx * vx + vy * vy + vz * vz);
    // if (vmag < 0.00000001) printf("vmag = %lf\n", vmag);

    wx = uy * vz - uz * vy;
    wy = uz * vx - ux * vz;
    wz = ux * vy - uy * vx;
    wmag = sqrt(wx * wx + wy * wy + wz * wz);

    delphi = (2.0 * ran3() - 1.0) * delphi_max;
    cosphi = cos(delphi);
    sinphi = sin(delphi);

    drx = ux + cosphi * vx + sinphi * vmag * wx / wmag;
    dry = uy + cosphi * vy + sinphi * vmag * wy / wmag;
    drz = uz + cosphi * vz + sinphi * vmag * wz / wmag;

    rx[k] = rx[nseg - 1] + drx;
    ry[k] = ry[nseg - 1] + dry;
    rz[k] = rz[nseg - 1] + drz;
  }
  else if (k == nseg4 - 1)
  {

    drx = rx[k] - rx[k - 1];
    dry = ry[k] - ry[k - 1];
    drz = rz[k] - rz[k - 1];

    Rx = rx[0] - rx[k - 1];
    Ry = ry[0] - ry[k - 1];
    Rz = rz[0] - rz[k - 1];
    Rmag = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);

    Rnx = Rx / Rmag;
    Rny = Ry / Rmag;
    Rnz = Rz / Rmag;

    rdotRn = drx * Rnx + dry * Rny + drz * Rnz;
    ux = rdotRn * Rnx;
    uy = rdotRn * Rny;
    uz = rdotRn * Rnz;

    vx = drx - ux;
    vy = dry - uy;
    vz = drz - uz;
    vmag = sqrt(vx * vx + vy * vy + vz * vz);
    // if (vmag < 0.00000001) printf("vmag = %lf\n", vmag);

    wx = uy * vz - uz * vy;
    wy = uz * vx - ux * vz;
    wz = ux * vy - uy * vx;
    wmag = sqrt(wx * wx + wy * wy + wz * wz);

    delphi = (2.0 * ran3() - 1.0) * delphi_max;
    cosphi = cos(delphi);
    sinphi = sin(delphi);

    drx = ux + cosphi * vx + sinphi * vmag * wx / wmag;
    dry = uy + cosphi * vy + sinphi * vmag * wy / wmag;
    drz = uz + cosphi * vz + sinphi * vmag * wz / wmag;

    rx[k] = rx[k - 1] + drx;
    ry[k] = ry[k - 1] + dry;
    rz[k] = rz[k - 1] + drz;
  }
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double **m;

  /* allocate pointers to rows */
  m = (double **)malloc((size_t)((nrow + NR_END) * sizeof(double *)));
  if (!m)
  {
    printf("allocation failure 1 in matrix()... exiting...");
    exit(1);
  }
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (double *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
  if (!m[nrl])
  {
    printf("allocation failure 2 in matrix()... exiting...");
    exit(1);
  }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}