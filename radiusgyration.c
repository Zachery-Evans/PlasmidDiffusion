/*
This code was written by Dr. James Polson, it has been edited and modified by Zach Evans during the summer of 2023

This program was quickly written to calculate the radius of gyration of a linear polymer, as such, there are some redundant if statements 
and confusing variable names that exist in a different context. 

Input parameters are to be changes in a file labeled "mc.inp", the parameters are as follows:

1) Length of the polymer 

2) Bending rigidity of the polymer (Physically relevant scales from 0 -> 10)

3) Number of Monte Carlo Cycles
4) Number of equilibrium cycles
5) Maximum angle that can be used for a crank shaft movement (default is PI)
6) Seed of the random number generator

7) Frequency of sampling the position of the polymer e.g. 1 snapshot every x moves
8) Frequency of sampling the center mass of the polymer (and Radius of Gyration)

9) 1 for writing an xyz file, 0 for omitting the xyz file. 
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.141592653589793
#define NR_END 1
#define FREE_ARG char *

// function heade_
double ran3(void);
void init_pos(void);
void init_pos_circular(void);
void write_log(void);
void write_data(void);
void input(void);
int check_accept(double[], double[], double[], long);

int check_poly_energy(double[], double[], double[], long);
double calc_cosine(int, int, int, double[], double[], double[]);

double **dmatrix(long, long, long, long);
void free_dmatrix(double **, long, long, long, long);

void reptation_move_chain1(void);

void shift_move_chain(void);

int check_accept_reptation(double[], double[], double[], long, long);
void calc_delta_xyz(void);

void crank_move_polymer(double[], double[], double[]);

long nseg1, nseg2, nseg3, nseg4, nbin, i, j, k, ii, ncyc, overlap, nacc, kk, itest, iseed;
long neq, nbintot, ibin, ichain, nsamp, nacc_shift, nshift;
long imov, kmaxtest, cmFreqSamp, freq_samp, ncmt, ngridx, ngridy;

double L, H, Ld2, Hd2, xt, yt, zt, dx, dy, dz, re, dr2, drxy2, dr2min, dr2max;
double xBoxMax, yBoxMax;
double kappa, xold, yold, zold, delphi_max;
double z1min, z1max, z2min, z2max, zcm1, zcm2, z1bcm, z2bcm;


FILE *fpmov;

double r1x[5000], r1y[5000], r1z[5000];

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
double Rgsum, Rgsq[10000];

// ------------------------------------------------------------------------
// main function
// ------------------------------------------------------------------------
int main()
{
  long imon, indx, indy;
  double Rgsq, Rgsq_avg;
  double xcm1, ycm1, zcm1;
  clock_t start, end;

  input();

  rep_prob = 0.95;

  nsamp = 0;

  write_log();

  // *************************************************************

  if (imov == 1)
  { // Don't include this in cluster
    fpmov = fopen("chain.xyz", "w");
    start = clock();
  }

  imon = 0;

  init_pos(); // function call

  // Don't include if statement below in cluster
  if (imov == 1)
  {
    if (ii % freq_samp == 0 && ii > neq)
    {
      fprintf(fpmov, "%ld\n", nseg1);
      fprintf(fpmov, "Polymer:  %ld\n", ii);

      for (i = 0; i < nseg1; i++)
      {
        fprintf(fpmov, "N    %lf  %lf  %lf\n", r1x[i], r1y[i], r1z[i]);
      }
    }
  }

  nacc = 0; // will hold number of accepted moves
  nacc_rep = 0;
  nacc_shift = 0;
  nshift = 0;
  nrep = 0;
  ichain = 1;
  for (ii = 0; ii < ncyc; ii++)
  {
    // if (ii % 100 == 0) printf("ii = %ld\n",ii);
    for (j = 0; j < nseg1 + nseg2 + nseg3 + nseg4; j++)
    {
      k = (nseg1)*ran3();

      kmaxtest = nseg1 - 2;

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
        }
      }
      else
      {
        nrep++;
        if (ichain == 1)
        {
          reptation_move_chain1();
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
        zcm1 += r1z[i];
      }
      xcm1 /= nseg1;
      ycm1 /= nseg1;
      zcm1 /= nseg1;

      nsamp += 1;
    }

    for (i = 0; i < nseg1; i++)
    {
      Rgsq += (r1x[i] - xcm1) * (r1x[i] - xcm1) + (r1y[i] - ycm1) * (r1y[i] - ycm1) + (r1z[i] - zcm1) * (r1z[i] - zcm1);
    }

    Rgsq = Rgsq / nseg1;

    if (ii % freq_samp == 0 && ii > neq)
    {
      Rgsq_avg += Rgsq;
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
      }
    }
  }

  printf("Acc. ratio = %lf\n", 1.0 * nacc / ((ncyc * (nseg1)) - nrep));
  printf("Number of reptation attempts = %ld\n", nrep);
  printf("Reptation Acc. ratio = %lf\n", 1.0 * nacc_rep / nrep);

  FILE *fp;

  if ((fp = fopen("rgsq.dat", "w")) == NULL)
  {
    printf("Cannot open file: rgsq.dat\n");
    exit(0);
  }
  else
  {
    for (i = 0; i < ngridx; i++)
    {
      for (j = 0; j < ngridy; j++)
      {
        fprintf(fp, "%lf\n", Rgsq_avg / (ncyc - neq));
      }
    }
    fclose(fp);
  }

  printf("Radius of Gyration = %lf\n", Rgsq_avg / (ncyc - neq));

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

    fscanf(fp, "%lf%*s", &kappa);

    fscanf(fp, "\n%ld%*s", &ncyc);
    fscanf(fp, "%ld%*s", &neq);
    fscanf(fp, "%lf%*s", &delphi_max);
    fscanf(fp, "%ld%*s", &iseed);

    fscanf(fp, "\n%ld%*s", &freq_samp);
    fscanf(fp, "\n%ld%*s", &cmFreqSamp);
    fscanf(fp, "\n%ld%*s", &freq_samp);

    fscanf(fp, "\n%ld%*s", &imov);
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
  printf("kappa    %lf\n", kappa);

  printf("\n");

  printf("ncyc     %ld\n", ncyc);
  printf("neq      %ld\n", neq);
  printf("delphi_max   %lf\n", delphi_max);
  printf("iseed    %ld\n", iseed);
  printf("\n");

  printf("freq_samp  %ld\n", freq_samp);
  printf("cmFreqSamp  %ld\n", cmFreqSamp);
  printf("\n");

  printf("imov     %ld\n", imov);
  printf("\n");
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
    }
    return (check_poly_energy(rx, ry, rz, nseg)); // apply rigidity
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
  double xadd, yadd, xmax, ymax, zplace;

  r1x[0] = - (double) nseg1 / 2.0 ;
  r1y[0] = 0.0;
  r1z[0] = 1.0;
  xadd = 1.0;
  yadd = 1.0;
  zplace = 1.0;

  for (i = 1; i < nseg1; i++)
  {
    r1x[i] = r1x[i - 1] + xadd;
    r1y[i] = r1y[i - 1];
    r1z[i] = zplace;
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
