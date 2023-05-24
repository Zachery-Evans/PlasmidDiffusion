#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define PI 3.141592653589793
#define NR_END 1
#define FREE_ARG char *

// function headers
double ran3(void);
void init_pos(void);
void write_log(void);
void write_data(void);
void input(void);
int check_accept(void);
int check_shift_chain2(void);

int check_energy(void);
double calc_cosine_chain1(int, int, int);
double calc_cosine_chain2(int, int, int);

double **dmatrix(long, long, long, long);
void free_dmatrix(double **, long, long, long, long);

void reptation_move_chain1(void);
void reptation_move_chain2(void);
void reptation_move_chain3(void);
void shift_move_chain(void);
void shift_move_chain2(void);
void shift_move_chain3(void);
int check_accept_reptation(long);
void calc_delta_xyz(void);
void crank_move_chain1(void);
void crank_move_chain2(void);
void crank_move_chain3(void);

long nseg1, nseg2, nseg3, nbin, i, j, k, ii, ncyc, overlap, nacc, kk, itest, iseed;
long neq, nbintot, ibin, ichain, nsamp, nacc_shift, nshift;
long imov, kmaxtest, freq_samp, freq_mon, freq_mov, ncmt, ngridx, ngridy;

double L, H, Ld2, Hd2, rmax, xt, yt, zt, dx, dy, dz, re, dr2, drxy2, dr2min, dr2max;
double qmin, qmax, re2av, re2, drmin, drmax, gridspace, gridspacex_real, gridspacey_real;
double amax, bmin, amax2, bmin2, ecc, Area, rectangleArea, rectangleXYRatio, xBoxMax, yBoxMax, rshift_max;
double xBoxMaxd2, yBoxMaxd2;
double kappa, xold, yold, zold, delphi_max;
double z1min, z1max, z2min, z2max, zcm1, zcm2, z1bcm, z2bcm;
double **prob1, **prob2, **probmon;

// FILE *fpmov;

double r1x[5000];
double r1y[5000];
double r1z[5000];
double r2x[5000];
double r2y[5000];
double r2z[5000];
double r3x[5000];
double r3y[5000];
double r3z[5000];

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

// ------------------------------------------------------------------------
// main function
// ------------------------------------------------------------------------
int main()
{
  long imon, indx, indy;
  double xcm1, ycm1, xcm2, ycm2;

  FILE *xp1, *yp1, *xp2, *yp2, *xp3, *yp3;

  input();

  rep_prob = 0.95;
  // defining the area of the ellipse as the total area subtracted by the area
  // of the rectangle now between two halves of the ellipse.

  // Requires that the y direction of the box is the same height as the semi-minor
  // axis of the ellipse.

  amax = sqrt((Area - rectangleArea) / PI) * pow(1.0 - ecc * ecc, -0.25);
  bmin = sqrt((Area - rectangleArea) / PI) * pow(1.0 - ecc * ecc, +0.25);
  amax2 = amax * amax;
  bmin2 = bmin * bmin;
  yBoxMaxd2 = bmin; // Width of the rectangle section equivalent to the semi-minor axis
  xBoxMax = rectangleArea / yBoxMaxd2 * 2.0;
  xBoxMaxd2 = xBoxMax / 2.0;

  // printf("Length of the box: %lf\n", xBoxMax);
  // printf("1/2 Length of the box: %lf\n", xBoxMaxd2);
  // printf("Semi-major axis: %lf\n", amax);
  // printf("Semi-minor axis: %lf\n", bmin);
  // printf("Height of box: %lf\n", yBoxMax);
  Hd2 = H / 2.0;

  ngridx = 2.0 * (amax + xBoxMaxd2) / gridspace + 0.00001;
  ngridy = 2.0 * bmin / gridspace + 0.00001;
  gridspacex_real = 2.0 * (amax + xBoxMaxd2) / ngridx;
  gridspacey_real = 2.0 * bmin / ngridy;
  // printf("ngridx = %ld, ngridy = %ld, gridspacex_real = %lf, gridspacey_real = %lf\n",
  // ngridx, ngridy, gridspacex_real, gridspacey_real);
  prob1 = dmatrix(0, ngridx - 1, 0, ngridy - 1);
  prob2 = dmatrix(0, ngridx - 1, 0, ngridy - 1);
  probmon = dmatrix(0, ngridx - 1, 0, ngridy - 1);

  for (i = 0; i < ngridx; i++)
  {
    for (j = 0; j < ngridy; j++)
    {
      prob1[i][j] = 0.0;
      prob2[i][j] = 0.0;
      probmon[i][j] = 0.0;
    }
  }

  nsamp = 0;

  dr2min = drmin * drmin;
  dr2max = drmax * drmax;

  write_log();

  // *************************************************************

  // if (imov == 1) // Don't include this in cluster
  // fpmov = fopen("chain.xyz", "w");

  imon = 0;

  init_pos(); // function call
              /* Don't include this in cluster
                if (imov == 1)
                {
                  if (ii % freq_mov == 0 && ii > -1)
                  {
                    fprintf(fpmov, "%ld\n", nseg1 + nseg2);
                    fprintf(fpmov, "polymer:  %ld\n", ii);
            
                    for (i = 0; i < nseg1; i++)
                    {
                      fprintf(fpmov, "N    %lf  %lf  %lf\n", r1x[i], r1y[i], r1z[i]);
                    }
                    for (i = 0; i < nseg2; i++)
                    {
                      fprintf(fpmov, "O    %lf  %lf  %lf\n", r2x[i], r2y[i], r2z[i]);
                    }
                  }
                }
              */
  nacc = 0;   // will hold number of accepted moves
  nacc_rep = 0;
  nacc_shift = 0;
  nshift = 0;
  nrep = 0;

  xp1 = fopen("xcm1.dat", "w");
  xp2 = fopen("xcm2.dat", "w");
  yp1 = fopen("ycm1.dat", "w");
  yp2 = fopen("ycm2.dat", "w");

  for (ii = 0; ii < ncyc; ii++)
  {
    // if (ii % 100 == 0) printf("ii = %ld\n",ii);

    for (j = 0; j < nseg1 + nseg2; j++)
    {

      k = (nseg1 + nseg2) * ran3();

      if (k < nseg1)
      {
        ichain = 1;
        kmaxtest = nseg1 - 2;
      }
      else
      {
        ichain = 2;
        k -= nseg1;
        kmaxtest = nseg2 - 2;
      }

      if (ran3() >= rep_prob && (k >= 2 || k < kmaxtest))
      {

        if (ichain == 1)
        {
          xold = r1x[k];
          yold = r1y[k];
          zold = r1z[k];
          crank_move_chain1();
        }
        else
        {
          xold = r2x[k];
          yold = r2y[k];
          zold = r2z[k];
          crank_move_chain2();
        }

        overlap = check_accept();

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
          else
          {
            r2x[k] = xold;
            r2y[k] = yold;
            r2z[k] = zold;
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
        else
        {
          shift_move_chain2();
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

      indx = (xcm1 + amax + xBoxMaxd2) / gridspacex_real;
      indy = (ycm1 + bmin) / gridspacey_real;
      /*
       * printf statements after the following conditions are met:
       * if (indx >= ngridx || indy >= ngridy)
       * Warning for invalid data
       */
      if (indx >= ngridx || indy >= ngridy)
        printf("1:  indx = %ld/%ld, indy = %ld/%ld\n", indx, ngridx, indy, ngridy);
      if (indx >= 0 && indx < ngridx && indy >= 0 && indy < ngridy)
        prob1[indx][indy] += 1.0;

      for (i = 0; i < nseg1; i++)
      {
        indx = (r1x[i] + amax + xBoxMaxd2) / gridspacex_real;
        indy = (r1y[i] + bmin) / gridspacey_real;
        if (indx >= 0 && indx < ngridx && indy >= 0 && indy < ngridy)
          probmon[indx][indy] += 1.0;
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

      indx = (xcm2 + amax + xBoxMaxd2) / gridspacex_real;
      indy = (ycm2 + bmin) / gridspacey_real;

      if (indx >= ngridx || indy >= ngridy)
      {
        printf("2:  indx = %ld/%ld, indy = %ld/%ld\n", indx, ngridx, indy, ngridy);
        printf("    xcm2 = %lf, 2*amax = %lf,  ycm2 = %lf, 2*bmin = %lf\n",
               xcm2 + amax / 2.0, 2 * amax, ycm2 + bmin / 2.0, 2 * bmin);
      }
      if (indx >= 0 && indx < ngridx && indy >= 0 && indy < ngridy)
        prob2[indx][indy] += 1.0;

      nsamp += 1;

      fprintf(xp1, "%lf\n", xcm1);
      fprintf(xp2, "%lf\n", xcm2);
      fprintf(yp1, "%lf\n", ycm1);
      fprintf(yp2, "%lf\n", ycm2);
    }
    /*
        if (imov == 1)
        {
          if (ii % freq_mov == 0 && ii > -1)
          {
            fprintf(fpmov, "%ld\n", nseg1 + nseg2);
            fprintf(fpmov, "polymer:  %ld\n", ii);

            for (i = 0; i < nseg1; i++)
            {
              fprintf(fpmov, "N    %lf  %lf  %lf\n", r1x[i], r1y[i], r1z[i]);
            }
            for (i = 0; i < nseg2; i++)
            {
              fprintf(fpmov, "O    %lf  %lf  %lf\n", r2x[i], r2y[i], r2z[i]);
            }
          }
        }
      */
  }

  fclose(xp1);
  fclose(xp2);
  fclose(yp1);
  fclose(yp2);

  printf("Acc. ratio = %lf\n", 1.0 * nacc / ((ncyc * (nseg1 + nseg2)) - nrep));
  printf("Number of reptation attempts = %ld\n", nrep);
  printf("Reptation Acc. ratio = %lf\n", 1.0 * nacc_rep / nrep);
  printf("Shift 2 Acc. ratio = %ld / %ld = %lf\n", nacc_shift, nshift, 1.0 * nacc_shift / nshift);

  write_data();

  /*
    if (imov == 1)
      fclose(fpmov);
      */
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
    fscanf(fp, "%ls%*s", &nseg3);
    fscanf(fp, "%lf%*s", &Area);
    fscanf(fp, "%lf%*s", &rectangleArea);
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
    fscanf(fp, "\n%ld%*s", &freq_mon);
    fscanf(fp, "\n%ld%*s", &freq_mov);

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
  printf("nseg2    %ld\n", nseg2);
  printf("nseg2    %ld\n", nseg3);
  printf("Area     %lf\n", Area);
  printf("Rectangle %lf\n", rectangleArea);
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
  printf("freq_mon   %ld\n", freq_mon);
  printf("freq_mov   %ld\n", freq_mov);
  printf("\n");

  printf("imov     %ld\n", imov);
  printf("\n");
}

/* squareEllipse and checkEllipse Code written by Zach Evans in attempt to create geometry of rectangle between two halves of an ellipse */

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

  if (xPos > xBoxMaxd2)
  { // If the polymer is outside of the leftmost semi-ellipse, reject
    if ((xPos - xBoxMaxd2) * (xPos - xBoxMaxd2) > amax2 * (1 - (yPos * yPos) / bmin2))
    {
      return reject;
    }
    if (yPos * yPos > bmin2 * (1 - (xPos - xBoxMaxd2) * (xPos - xBoxMaxd2) / amax2))
    {
      return reject;
    }

    echeck = (((xPos - xBoxMaxd2) * (xPos - xBoxMaxd2)) / amax2) + ((yPos * yPos) / bmin2);

    if (echeck > 1.0)
    {
      return (reject);
    }
  }

  else if (xPos < -xBoxMaxd2)
  { // Checking if outside of left elliptical end
    if ((xPos + xBoxMaxd2) * (xPos + xBoxMaxd2) > amax2 * (1 - (yPos * yPos) / bmin2))
    {
      return reject;
    }
    if (yPos * yPos > bmin2 * (1 - (xPos + xBoxMaxd2) * (xPos + xBoxMaxd2) / amax2))
    {
      return reject;
    }

    echeck = ((xPos + xBoxMaxd2) * (xPos + xBoxMaxd2) / amax2) + (yPos * yPos / bmin2);

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
  else if (yPos > yBoxMaxd2 || yPos < -yBoxMaxd2)
  {
    return (reject);
  }

  return accept;
}

// ----------------------------------------------------------------------
// This function determines whether the move just performed in the main
// function is permitted or not (checks overlap with a number of conditions)
// ----------------------------------------------------------------------
int check_accept(void)
{
  int accept, reject;
  long klow, khigh;

  accept = 0;
  reject = 1;

  if (ichain == 1)
  {
    /* Checked in squareEllipse
    if (r1z[k] < -Hd2 || r1z[k] > Hd2)
    {
      return (reject);
    }
    */

    /* This eccentricity check does not apply in the region of the square.
    echeck = r1x[k] * r1x[k] / amax2 + r1y[k] * r1y[k] / bmin2;

    if (echeck > 1.0)
    {
      return (reject);
    }
    */

    /*
        for (kk=k-1;kk<=k+1;kk+=2) {
          if (kk >= 0 && kk < nseg1) {
            dx = r1x[k] - r1x[kk];
            dy = r1y[k] - r1y[kk];
            dz = r1z[k] - r1z[kk];
            dr2 = dx*dx + dy*dy + dz*dz;
            if (dr2 < dr2min || dr2 > dr2max) return(reject);
          }
        }
    */

    // Checking if the T4 polymer overlaps with itself
    for (kk = 0; kk < nseg1; kk++)
    {
      if (squareEllipse(r1x[kk], r1y[kk], r1z[kk]) == reject)
      {
        return (reject);
      }
      if (kk < k - 1 || kk > k + 1)
      {
        dx = r1x[k] - r1x[kk];
        dy = r1y[k] - r1y[kk];
        dz = r1z[k] - r1z[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
          return (reject);
      }
    }
    // Checking if the plasmid overlaps with the T4 polymer
    for (kk = 0; kk < nseg2; kk++)
    {
      dx = r1x[k] - r2x[kk];
      dy = r1y[k] - r2y[kk];
      dz = r1z[k] - r2z[kk];
      dr2 = dx * dx + dy * dy + dz * dz;
      if (dr2 < 1.0)
        return (reject);
    }
  }
  else
  {
    /* squareEllipse checks this
    if (r2z[k] < -Hd2 || r2z[k] > Hd2)
    {
      return (reject);
    }
    */

    /*
    echeck = r2x[k] * r2x[k] / amax2 + r2y[k] * r2y[k] / bmin2;
    if (echeck > 1.0)
      return (reject);


      for (kk=k-1;kk<=k+1;kk+=2) {
        if (kk >= 0 && kk < nseg2) {
          dx = r2x[k] - r2x[kk];
          dy = r2y[k] - r2y[kk];
          dz = r2z[k] - r2z[kk];
          dr2 = dx*dx + dy*dy + dz*dz;
          if (dr2 < dr2min || dr2 > dr2max) return(reject);
        }
      }
    */

    if (k == 0)
    {
      klow = nseg2 - 1;
      khigh = 1;
    }
    else if (k == nseg2 - 1)
    {
      klow = nseg2 - 2;
      khigh = 0;
    }
    else
    {
      klow = k - 1;
      khigh = k + 1;
    }

    for (kk = 0; kk < nseg2; kk++)
    {
      if (squareEllipse(r2x[kk], r2y[kk], r2z[kk]) == reject)
      {
        return (reject);
      }
      if (kk != k && kk != klow && kk != khigh)
      {
        dx = r2x[k] - r2x[kk];
        dy = r2y[k] - r2y[kk];
        dz = r2z[k] - r2z[kk];
        dr2 = dx * dx + dy * dy + dz * dz;

        if (dr2 < 1.0)
        {
          return (reject);
        }
      }
    }

    for (kk = 0; kk < nseg1; kk++)
    {
      dx = r2x[k] - r1x[kk];
      dy = r2y[k] - r1y[kk];
      dz = r2z[k] - r1z[kk];
      dr2 = dx * dx + dy * dy + dz * dz;
      if (dr2 < 1.0)
        return (reject);
    }
  }

  if (ichain == 1)
  {
    return (check_energy());
  }
  else
  {
    return (accept);
  }
}

// ----------------------------------------------------------------------
// This function is called to determine whether or not a move is accepted
// based on considerations of polymer chain energy related to the
// move. See comments within the function for details.
// ----------------------------------------------------------------------
int check_energy(void)
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

  // if considering a movement in the first chain
  if (ichain == 1)
  {
    // The following 5 blocks consider the possible scenarios of monomer
    // movement: 1.) first monomer, 2.) last monomer, 3.) second monomer,
    // 4.) second last monomer, and 5.) any other monomer
    //
    // shifting monomers at the ends requires only one change in angle
    // and thus energy. Monomers second from the end require two, and
    // all the rest require three.
    if (k == 0)
    {
      theta_new = calc_cosine_chain1(k, k + 1, k + 2);
      theta_old = calc_cosine_chain1(-1, k + 1, k + 2);
      energy_new[0] = kappa * (1.0 - theta_new);
      energy_old[0] = kappa * (1.0 - theta_old);
    }
    else if (k == nseg1 - 1)
    {
      theta_new = calc_cosine_chain1(k - 2, k - 1, k);
      theta_old = calc_cosine_chain1(k - 2, k - 1, -1);
      energy_new[0] = kappa * (1.0 - theta_new);
      energy_old[0] = kappa * (1.0 - theta_old);
    }
    else if (k == 1)
    {
      theta_new = calc_cosine_chain1(k - 1, k, k + 1);
      theta_old = calc_cosine_chain1(k - 1, -1, k + 1);
      energy_new[0] = kappa * (1.0 - theta_new);
      energy_old[0] = kappa * (1.0 - theta_old);

      theta_new = calc_cosine_chain1(k, k + 1, k + 2);
      theta_old = calc_cosine_chain1(-1, k + 1, k + 2);
      energy_new[1] = kappa * (1.0 - theta_new);
      energy_old[1] = kappa * (1.0 - theta_old);
    }
    else if (k == nseg1 - 2)
    {
      theta_new = calc_cosine_chain1(k - 2, k - 1, k);
      theta_old = calc_cosine_chain1(k - 2, k - 1, -1);
      energy_new[0] = kappa * (1.0 - theta_new);
      energy_old[0] = kappa * (1.0 - theta_old);

      theta_new = calc_cosine_chain1(k - 1, k, k + 1);
      theta_old = calc_cosine_chain1(k - 1, -1, k + 1);
      energy_new[1] = kappa * (1.0 - theta_new);
      energy_old[1] = kappa * (1.0 - theta_old);
    }
    else
    {
      theta_new = calc_cosine_chain1(k - 2, k - 1, k);
      theta_old = calc_cosine_chain1(k - 2, k - 1, -1);
      energy_new[0] = kappa * (1.0 - theta_new);
      energy_old[0] = kappa * (1.0 - theta_old);

      theta_new = calc_cosine_chain1(k - 1, k, k + 1);
      theta_old = calc_cosine_chain1(k - 1, -1, k + 1);
      energy_new[1] = kappa * (1.0 - theta_new);
      energy_old[1] = kappa * (1.0 - theta_old);

      theta_new = calc_cosine_chain1(k, k + 1, k + 2);
      theta_old = calc_cosine_chain1(-1, k + 1, k + 2);
      energy_new[2] = kappa * (1.0 - theta_new);
      energy_old[2] = kappa * (1.0 - theta_old);
    }
  }
  // if considering a movement in the second chain:
  else
  {
    if (k == 0)
    {
      theta_new = calc_cosine_chain2(k, k + 1, k + 2);
      theta_old = calc_cosine_chain2(-1, k + 1, k + 2);
      energy_new[0] = kappa * (1.0 - theta_new);
      energy_old[0] = kappa * (1.0 - theta_old);
    }
    else if (k == nseg2 - 1)
    {
      theta_new = calc_cosine_chain2(k - 2, k - 1, k);
      theta_old = calc_cosine_chain2(k - 2, k - 1, -1);
      energy_new[0] = kappa * (1.0 - theta_new);
      energy_old[0] = kappa * (1.0 - theta_old);
    }
    else if (k == 1)
    {
      theta_new = calc_cosine_chain2(k - 1, k, k + 1);
      theta_old = calc_cosine_chain2(k - 1, -1, k + 1);
      energy_new[0] = kappa * (1.0 - theta_new);
      energy_old[0] = kappa * (1.0 - theta_old);

      theta_new = calc_cosine_chain2(k, k + 1, k + 2);
      theta_old = calc_cosine_chain2(-1, k + 1, k + 2);
      energy_new[1] = kappa * (1.0 - theta_new);
      energy_old[1] = kappa * (1.0 - theta_old);
    }
    else if (k == nseg2 - 2)
    {
      theta_new = calc_cosine_chain2(k - 2, k - 1, k);
      theta_old = calc_cosine_chain2(k - 2, k - 1, -1);
      energy_new[0] = kappa * (1.0 - theta_new);
      energy_old[0] = kappa * (1.0 - theta_old);

      theta_new = calc_cosine_chain2(k - 1, k, k + 1);
      theta_old = calc_cosine_chain2(k - 1, -1, k + 1);
      energy_new[1] = kappa * (1.0 - theta_new);
      energy_old[1] = kappa * (1.0 - theta_old);
    }
    else
    {
      theta_new = calc_cosine_chain2(k - 2, k - 1, k);
      theta_old = calc_cosine_chain2(k - 2, k - 1, -1);
      energy_new[0] = kappa * (1.0 - theta_new);
      energy_old[0] = kappa * (1.0 - theta_old);

      theta_new = calc_cosine_chain2(k - 1, k, k + 1);
      theta_old = calc_cosine_chain2(k - 1, -1, k + 1);
      energy_new[1] = kappa * (1.0 - theta_new);
      energy_old[1] = kappa * (1.0 - theta_old);

      theta_new = calc_cosine_chain2(k, k + 1, k + 2);
      theta_old = calc_cosine_chain2(-1, k + 1, k + 2);
      energy_new[2] = kappa * (1.0 - theta_new);
      energy_old[2] = kappa * (1.0 - theta_old);
    }
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

// ----------------------------------------------------------------------
// Calculate the cosine of theta between two bonds (vectors) that connect
// three sequential monomers within the first polymer chain.
// ----------------------------------------------------------------------
double calc_cosine_chain1(int i1, int i2, int i3)
{
  // determine which of the indices is negative (if any). The negative index
  // corresponds to the old position of the monomer. Set x, y, and z values of
  // each of the three monomers accordingly. Note y1 is an implicit variable
  // name in some stlib function, so use yone instead.
  if (i1 < 0)
  {
    x1 = xold;
    yone = yold, z1 = zold;
    x2 = r1x[i2];
    y2 = r1y[i2];
    z2 = r1z[i2];
    x3 = r1x[i3];
    y3 = r1y[i3];
    z3 = r1z[i3];
  }
  else if (i2 < 0)
  {
    x1 = r1x[i1];
    yone = r1y[i1];
    z1 = r1z[i1];
    x2 = xold;
    y2 = yold;
    z2 = zold;
    x3 = r1x[i3];
    y3 = r1y[i3];
    z3 = r1z[i3];
  }
  else if (i3 < 0)
  {
    x1 = r1x[i1];
    yone = r1y[i1];
    z1 = r1z[i1];
    x2 = r1x[i2];
    y2 = r1y[i2];
    z2 = r1z[i2];
    x3 = xold;
    y3 = yold;
    z3 = zold;
  }
  else
  {
    x1 = r1x[i1];
    yone = r1y[i1];
    z1 = r1z[i1];
    x2 = r1x[i2];
    y2 = r1y[i2];
    z2 = r1z[i2];
    x3 = r1x[i3];
    y3 = r1y[i3];
    z3 = r1z[i3];
  }

  // va is vector A, connects the first two monomers being considered.
  // vb is vector B, connects the second two monomers being considered.
  // Use rules of dot products to calculate the cosine of the angle between
  // vectors A and B, and return this value.

  vax = x2 - x1;
  vay = y2 - yone;
  vaz = z2 - z1;
  vbx = x3 - x2;
  vby = y3 - y2;
  vbz = z3 - z2;

  va_sq = vax * vax + vay * vay + vaz * vaz;
  vb_sq = vbx * vbx + vby * vby + vbz * vbz;

  va_dot_vb = vax * vbx + vay * vby + vaz * vbz;

  return (va_dot_vb / (sqrt(va_sq) * sqrt(vb_sq)));
}

// ----------------------------------------------------------------------
// Calculate the cosine of theta between two bonds (vectors) that connect
// three sequential monomers within the second polymer chain.
// ----------------------------------------------------------------------
double calc_cosine_chain2(int i1, int i2, int i3)
{
  if (i1 < 0)
  {
    x1 = xold;
    yone = yold, z1 = zold;
    x2 = r2x[i2];
    y2 = r2y[i2];
    z2 = r2z[i2];
    x3 = r2x[i3];
    y3 = r2y[i3];
    z3 = r2z[i3];
  }
  else if (i2 < 0)
  {
    x1 = r2x[i1];
    yone = r2y[i1];
    z1 = r2z[i1];
    x2 = xold;
    y2 = yold;
    z2 = zold;
    x3 = r2x[i3];
    y3 = r2y[i3];
    z3 = r2z[i3];
  }
  else if (i3 < 0)
  {
    x1 = r2x[i1];
    yone = r2y[i1];
    z1 = r2z[i1];
    x2 = r2x[i2];
    y2 = r2y[i2];
    z2 = r2z[i2];
    x3 = xold;
    y3 = yold;
    z3 = zold;
  }
  else
  {
    x1 = r2x[i1];
    yone = r2y[i1];
    z1 = r2z[i1];
    x2 = r2x[i2];
    y2 = r2y[i2];
    z2 = r2z[i2];
    x3 = r2x[i3];
    y3 = r2y[i3];
    z3 = r2z[i3];
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

  return (va_dot_vb / (sqrt(va_sq) * sqrt(vb_sq)));
}

// ----------------------------------------------------------------------
//  Places polymer in original (unrealistic position); is called for each
//  unique window. Overlaps the polymers somewhere to match the current
//  window.
// ----------------------------------------------------------------------
void init_pos(void)
{

  double xadd, yadd, bond_length, xmax, ymax;

  bond_length = 1.0;

  r1x[0] = 0.0;
  r1y[0] = -bmin / 2.0 + 1.0;
  r1z[0] = 1.0;
  xadd = 1.0;
  yadd = 1.02;

  for (i = 1; i < nseg1; i++)
  {
    r1x[i] = r1x[i - 1] + xadd;
    r1y[i] = r1y[i - 1];
    xmax = xBoxMaxd2; // Changed to be inside of the rectangle structure
    // xmax = amax * sqrt(1.0 - pow(r1y[i] / bmin, 2.0)) - 3.0;
    if (r1x[i] > xmax || r1x[i] < -xmax)
    {
      r1x[i] -= xadd;
      r1y[i] = r1y[i] + yadd;
      ymax = yBoxMaxd2;
      // ymax = bmin * sqrt(1.0 - pow(r1x[i] / amax, 2.0));
      if (r1y[i] > ymax || r1y[i] < -ymax)
      {
        printf("Can't place polymer... exiting...\n");
        exit(0);
      }
      xadd *= -1.0;
    }
    r1z[i] = 1.0;
  }

  double theta_plasmid = 2.0 * PI / nseg2;
  double Rplasmid = 0.5 / tan(theta_plasmid / 2.0);
  for (i = 0; i < nseg2; i++)
  {
    r2z[i] = -1.0;
    r2x[i] = Rplasmid * cos(i * theta_plasmid);
    r2y[i] = Rplasmid * sin(i * theta_plasmid);
  }
}

void write_data(void)
{
  FILE *fp;

  if ((fp = fopen("prob1.dat", "w")) == NULL)
  {
    printf("Cannot open file: prob1.dat\n");
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

  if ((fp = fopen("prob2.dat", "w")) == NULL)
  {
    printf("Cannot open file: prob2.dat\n");
  }
  else
  {
    for (i = 0; i < ngridx; i++)
    {
      for (j = 0; j < ngridy; j++)
        fprintf(fp, "%8.2lf  ", prob2[i][j]);
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  if ((fp = fopen("probmon.dat", "w")) == NULL)
  {
    printf("Cannot open file: probmon.dat\n");
  }
  else
  {
    for (i = 0; i < ngridx; i++)
    {
      for (j = 0; j < ngridy; j++)
        fprintf(fp, "%8.2lf  ", probmon[i][j]);
      fprintf(fp, "\n");
    }
    fclose(fp);
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

void shift_move_chain2()
{
  double delrx, delry, delrz;
  double xsold[5000];
  double ysold[5000];
  double zsold[5000];

  for (i = 0; i < nseg2; i++)
  {
    xsold[i] = r2x[i];
    ysold[i] = r2y[i];
    zsold[i] = r2z[i];
  }

  delrx = rshift_max * (2.0 * ran3() - 1.0);
  delry = rshift_max * (2.0 * ran3() - 1.0);
  delrz = rshift_max * (2.0 * ran3() - 1.0);

  for (i = 0; i < nseg2; i++)
  {
    r2x[i] += delrx;
    r2y[i] += delry;
    r2z[i] += delrz;
  }

  // printf("delrx = %lf, delry = %lf, delrz = %lf\n",delrx,delry,delrz);

  overlap = check_shift_chain2();

  nshift += 1;
  if (overlap == 0)
  {
    nacc += 1;
    nacc_shift += 1;
  }
  else if (overlap == 1)
  {
    for (i = 0; i < nseg2; i++)
    {
      r2x[i] = xsold[i];
      r2y[i] = ysold[i];
      r2z[i] = zsold[i];
    }
  }
}

int check_shift_chain2()
{
  int accept = 0;
  int reject = 1;
  double echeck;

  for (i = 0; i < nseg2; i++)
  {
    /* r2z check done inside of squareEllipse

    if (r2z[i] < -Hd2 || r2z[i] > Hd2)
    {
      return (reject);
    }

    //This echeck no longer relevant for this shape.

    echeck = r2x[i] * r2x[i] / amax2 + r2y[i] * r2y[i] / bmin2;
    if (echeck > 1.0)
      return (reject);
    */

    for (kk = 0; kk < nseg1; kk++)
    {
      if (squareEllipse(r2x[kk], r2y[kk], r2z[kk]) == reject)
      {
        return (reject);
      }

      dx = r2x[i] - r1x[kk];
      dy = r2y[i] - r1y[kk];
      dz = r2z[i] - r1z[kk];
      dr2 = dx * dx + dy * dy + dz * dz;
      if (dr2 < 1.0)
        return (reject);
    }
  }

  return (accept);
}

void reptation_move_chain1()
{
  double rannum;
  rannum = ran3();
  if (rannum <= 0.5)
    irep = 0;
  else
    irep = nseg1 - 1;

  rannum = ran3();
  phi_prime = rannum * 2.0 * PI;

  rannum = ran3();
  if (kappa > -0.000001 && kappa < 0.000001)
    costheta_prime = (2.0 * rannum) - 1.0;
  else
    costheta_prime = log((rannum * exp(kappa)) + ((1.0 - rannum) * exp(-1.0 * kappa))) / kappa;

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

    overlap = check_accept_reptation(nseg1 - 1);

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

    overlap = check_accept_reptation(0);

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

void reptation_move_chain2()
{
  double rannum;
  rannum = ran3();
  if (rannum <= 0.5)
  {
    irep = 0;
  }
  else
  {
    irep = nseg2 - 1;
  }
  rannum = ran3();
  phi_prime = rannum * 2 * PI;

  rannum = ran3();
  if (kappa > -0.000001 && kappa < 0.000001)
    costheta_prime = (2.0 * rannum) - 1.0;
  else
    costheta_prime = log((rannum * exp(kappa)) + ((1.0 - rannum) * exp(-1.0 * kappa))) / kappa;

  dr_prime = pow((0.602 * rannum) + 0.729, (1.0 / 3.0));

  //
  //      keep bond length = 1
  //
  //      rannum = ran3();
  //      dr_prime = pow((0.602*rannum)+0.729,(1.0/3.0));
  dr_prime = 1.0;

  xold = r2x[irep];
  yold = r2y[irep];
  zold = r2z[irep];

  calc_delta_xyz();

  if (irep == 0)
  {
    for (ind = 0; ind < nseg2 - 1; ind++)
    {
      r2x[ind] = r2x[ind + 1];
      r2y[ind] = r2y[ind + 1];
      r2z[ind] = r2z[ind + 1];
    }

    r2x[nseg2 - 1] = r2x[nseg2 - 2] + dx_fixed;
    r2y[nseg2 - 1] = r2y[nseg2 - 2] + dy_fixed;
    r2z[nseg2 - 1] = r2z[nseg2 - 2] + dz_fixed;

    overlap = check_accept_reptation(nseg2 - 1);

    if (overlap == 0)
    {
      nacc_rep += 1;
    }
    else
    {
      for (ind = nseg2 - 1; ind > 0; ind--)
      {
        r2x[ind] = r2x[ind - 1];
        r2y[ind] = r2y[ind - 1];
        r2z[ind] = r2z[ind - 1];
      }

      r2x[0] = xold;
      r2y[0] = yold;
      r2z[0] = zold;
    }
  }
  else
  { // irep == nseg2-1
    for (ind = nseg2 - 1; ind > 0; ind--)
    {
      r2x[ind] = r2x[ind - 1];
      r2y[ind] = r2y[ind - 1];
      r2z[ind] = r2z[ind - 1];
    }

    r2x[0] = r2x[1] + dx_fixed;
    r2y[0] = r2y[1] + dy_fixed;
    r2z[0] = r2z[1] + dz_fixed;

    overlap = check_accept_reptation(0);

    if (overlap == 0)
    {
      nacc_rep += 1;
    }
    else
    {
      for (ind = 0; ind < nseg2 - 1; ind++)
      {
        r2x[ind] = r2x[ind + 1];
        r2y[ind] = r2y[ind + 1];
        r2z[ind] = r2z[ind + 1];
      }

      r2x[nseg2 - 1] = xold;
      r2y[nseg2 - 1] = yold;
      r2z[nseg2 - 1] = zold;
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
  else
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

  /*beta = atan(ux / uz);
  alpha = atan(uy / ux);*/

  /*if (uz > -0.000001 && uz < 0.000001) {
    beta = PI / 2.0;
  }
  else if (ux > -0.000001 && ux < 0.000001) {
    if (uz >= 0.0)
      beta = 0.0;
    else
      beta = PI;
  }
  else if (uz >= 0.0) {
    beta = atan(fabs(ux)/ uz);
  }
  else {
    beta = atan(fabs(uz) / fabs(ux)) + (PI / 2.0);
  }*/

  /*if (ux > -0.000001 && ux < 0.000001) {
    if (uy >= 0.0)
      alpha = PI / 2.0;
    else
      alpha = 3.0 * PI / 2.0;
  }
  else if (uy > -0.000001 && uy < 0.000001) {
    if (ux >= 0.0)
      alpha = 0.0;
    else
      alpha = PI;
  }
  else if (ux >= 0.0 && uy >= 0.0)
    alpha = atan(uy/ux);
  else if (ux < 0.0 && uy >= 0.0)
    alpha = atan (fabs(ux)/uy) + (PI / 2.0);
  else if (ux < 0.0 && uy < 0.0)
    alpha = atan(uy/ux) + PI;
  else
    alpha = atan(ux/fabs(uy)) + (3.0 * PI / 2.0);
  */

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

  /*cos_beta = cos(beta);
  cos_alpha = cos(alpha);
  sin_beta = sin(beta);
  sin_alpha = sin(alpha);*/

  // Incorret order of rotations
  /*dx_fixed = (cos_beta*cos_alpha*dx_prime) + (cos_beta*sin_alpha*dy_prime) - (sin_beta*dz_prime);
  dy_fixed = (cos_alpha*dy_prime) - (sin_alpha*dx_prime);
  dz_fixed = (sin_beta*cos_alpha*dx_prime) + (sin_beta*sin_alpha*dy_prime) + (cos_beta*dz_prime);*/

  // Inverted matrix with negative angles
  /*dx_fixed = (cos_beta*cos_alpha*dx_prime) + (sin_alpha*dy_prime) - (sin_beta*cos_alpha*dz_prime);
  dy_fixed = (cos_alpha*dy_prime) - (cos_beta*sin_alpha*dx_prime) + (sin_beta*sin_alpha*dz_prime);
  dz_fixed = (cos_beta*dz_prime) + (sin_beta*dx_prime);*/

  // inverted matrix
  dx_fixed = (cos_beta * cos_alpha * dx_prime) - (sin_alpha * dy_prime) + (sin_beta * cos_alpha * dz_prime);
  dy_fixed = (cos_beta * sin_alpha * dx_prime) + (cos_alpha * dy_prime) + (sin_beta * sin_alpha * dz_prime);
  dz_fixed = (cos_beta * dz_prime) - (sin_beta * dx_prime);
}

int check_accept_reptation(long krep)
{
  int accept, reject; // will return either accept or reject at end of function
  double echeck;

  accept = 0;
  reject = 1;

  if (ichain == 1)
  {
    /*
    if (r1z[krep] < -Hd2 || r1z[krep] > Hd2)
      return (reject);
    echeck = r1x[krep] * r1x[krep] / amax2 + r1y[krep] * r1y[krep] / bmin2;
    if (echeck > 1.0)
      return (reject);
    */

    for (kk = 0; kk < nseg1; kk++)
    {

      if (squareEllipse(r1x[kk], r1y[kk], r1z[kk]) == reject)
      {
        return (reject);
      }

      if (kk < krep - 1 || kk > krep + 1)
      {
        dz = r1z[krep] - r1z[kk];
        if (fabs(dz) < 1.0)
        {
          dx = r1x[krep] - r1x[kk];
          dy = r1y[krep] - r1y[kk];
          dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 < 1.0)
            return (reject); // if overlap with other monomer within chain, reject
        }
      }
    }

    for (kk = 0; kk < nseg2; kk++)
    {
      dz = r1z[krep] - r2z[kk];
      if (fabs(dz) < 1.0)
      {
        dx = r1x[krep] - r2x[kk];
        dy = r1y[krep] - r2y[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
          return (reject); // if overlap with monomer in other chain, reject
      }
    }
  }
  else
  { // ichain == 2

    /*
    if (r2z[krep] < -Hd2 || r2z[krep] > Hd2)
      return (reject);
    echeck = r2x[krep] * r2x[krep] / amax2 + r2y[krep] * r2y[krep] / bmin2;
    if (echeck > 1.0)
      return (reject);
    */

    for (kk = 0; kk < nseg2; kk++)
    {
      if (squareEllipse(r2x[kk], r2y[kk], r2z[kk]) == reject)
      {
        return (reject);
      }
      if (kk < krep - 1 || kk > krep + 1)
      {
        dz = r2z[krep] - r2z[kk];
        if (fabs(dz) < 1.0)
        {
          dx = r2x[krep] - r2x[kk];
          dy = r2y[krep] - r2y[kk];
          dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 < 1.0)
            return (reject);
        }
      }
    }

    for (kk = 0; kk < nseg1; kk++)
    {
      dz = r2z[krep] - r1z[kk];
      if (fabs(dz) < 1.0)
      {
        dx = r2x[krep] - r1x[kk];
        dy = r2y[krep] - r1y[kk];
        dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 < 1.0)
          return (reject);
      }
    }
  }

  return (accept);
}

void crank_move_chain1()
{

  double rx, ry, rz, Rx, Ry, Rz, Rmag, rdotRn, Rnx, Rny, Rnz;
  double ux, uy, uz, vx, vy, vz, vmag, wx, wy, wz, wmag;
  double cosphi, sinphi, delphi;

  rx = r1x[k] - r1x[k - 1];
  ry = r1y[k] - r1y[k - 1];
  rz = r1z[k] - r1z[k - 1];

  Rx = r1x[k + 1] - r1x[k - 1];
  Ry = r1y[k + 1] - r1y[k - 1];
  Rz = r1z[k + 1] - r1z[k - 1];
  Rmag = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);

  Rnx = Rx / Rmag;
  Rny = Ry / Rmag;
  Rnz = Rz / Rmag;

  rdotRn = rx * Rnx + ry * Rny + rz * Rnz;
  ux = rdotRn * Rnx;
  uy = rdotRn * Rny;
  uz = rdotRn * Rnz;

  vx = rx - ux;
  vy = ry - uy;
  vz = rz - uz;
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

    rx = ux + cosphi * vx + sinphi * vmag * wx / wmag;
    ry = uy + cosphi * vy + sinphi * vmag * wy / wmag;
    rz = uz + cosphi * vz + sinphi * vmag * wz / wmag;

    r1x[k] = r1x[k - 1] + rx;
    r1y[k] = r1y[k - 1] + ry;
    r1z[k] = r1z[k - 1] + rz;
  }
  else
  { // bonds are parallel

    r1x[k] = xold;
    r1y[k] = yold;
    r1z[k] = zold;
  }
}

void crank_move_chain2()
{

  double rx, ry, rz, Rx, Ry, Rz, Rmag, rdotRn, Rnx, Rny, Rnz;
  double ux, uy, uz, vx, vy, vz, vmag, wx, wy, wz, wmag;
  double cosphi, sinphi, delphi;

  if (k > 0 && k < nseg2 - 1)
  {
    rx = r2x[k] - r2x[k - 1];
    ry = r2y[k] - r2y[k - 1];
    rz = r2z[k] - r2z[k - 1];

    Rx = r2x[k + 1] - r2x[k - 1];
    Ry = r2y[k + 1] - r2y[k - 1];
    Rz = r2z[k + 1] - r2z[k - 1];
    Rmag = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);

    Rnx = Rx / Rmag;
    Rny = Ry / Rmag;
    Rnz = Rz / Rmag;

    rdotRn = rx * Rnx + ry * Rny + rz * Rnz;
    ux = rdotRn * Rnx;
    uy = rdotRn * Rny;
    uz = rdotRn * Rnz;

    vx = rx - ux;
    vy = ry - uy;
    vz = rz - uz;
    vmag = sqrt(vx * vx + vy * vy + vz * vz);
    // if (vmag < 0.00000001) printf("vmag = %lf\n", vmag);

    wx = uy * vz - uz * vy;
    wy = uz * vx - ux * vz;
    wz = ux * vy - uy * vx;
    wmag = sqrt(wx * wx + wy * wy + wz * wz);

    delphi = (2.0 * ran3() - 1.0) * delphi_max;
    cosphi = cos(delphi);
    sinphi = sin(delphi);

    rx = ux + cosphi * vx + sinphi * vmag * wx / wmag;
    ry = uy + cosphi * vy + sinphi * vmag * wy / wmag;
    rz = uz + cosphi * vz + sinphi * vmag * wz / wmag;

    r2x[k] = r2x[k - 1] + rx;
    r2y[k] = r2y[k - 1] + ry;
    r2z[k] = r2z[k - 1] + rz;
  }
  else if (k == 0)
  {

    rx = r2x[k] - r2x[nseg2 - 1];
    ry = r2y[k] - r2y[nseg2 - 1];
    rz = r2z[k] - r2z[nseg2 - 1];

    Rx = r2x[k + 1] - r2x[nseg2 - 1];
    Ry = r2y[k + 1] - r2y[nseg2 - 1];
    Rz = r2z[k + 1] - r2z[nseg2 - 1];
    Rmag = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);

    Rnx = Rx / Rmag;
    Rny = Ry / Rmag;
    Rnz = Rz / Rmag;

    rdotRn = rx * Rnx + ry * Rny + rz * Rnz;
    ux = rdotRn * Rnx;
    uy = rdotRn * Rny;
    uz = rdotRn * Rnz;

    vx = rx - ux;
    vy = ry - uy;
    vz = rz - uz;
    vmag = sqrt(vx * vx + vy * vy + vz * vz);
    // if (vmag < 0.00000001) printf("vmag = %lf\n", vmag);

    wx = uy * vz - uz * vy;
    wy = uz * vx - ux * vz;
    wz = ux * vy - uy * vx;
    wmag = sqrt(wx * wx + wy * wy + wz * wz);

    delphi = (2.0 * ran3() - 1.0) * delphi_max;
    cosphi = cos(delphi);
    sinphi = sin(delphi);

    rx = ux + cosphi * vx + sinphi * vmag * wx / wmag;
    ry = uy + cosphi * vy + sinphi * vmag * wy / wmag;
    rz = uz + cosphi * vz + sinphi * vmag * wz / wmag;

    r2x[k] = r2x[nseg2 - 1] + rx;
    r2y[k] = r2y[nseg2 - 1] + ry;
    r2z[k] = r2z[nseg2 - 1] + rz;
  }
  else if (k == nseg2 - 1)
  {

    rx = r2x[k] - r2x[k - 1];
    ry = r2y[k] - r2y[k - 1];
    rz = r2z[k] - r2z[k - 1];

    Rx = r2x[0] - r2x[k - 1];
    Ry = r2y[0] - r2y[k - 1];
    Rz = r2z[0] - r2z[k - 1];
    Rmag = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);

    Rnx = Rx / Rmag;
    Rny = Ry / Rmag;
    Rnz = Rz / Rmag;

    rdotRn = rx * Rnx + ry * Rny + rz * Rnz;
    ux = rdotRn * Rnx;
    uy = rdotRn * Rny;
    uz = rdotRn * Rnz;

    vx = rx - ux;
    vy = ry - uy;
    vz = rz - uz;
    vmag = sqrt(vx * vx + vy * vy + vz * vz);
    // if (vmag < 0.00000001) printf("vmag = %lf\n", vmag);

    wx = uy * vz - uz * vy;
    wy = uz * vx - ux * vz;
    wz = ux * vy - uy * vx;
    wmag = sqrt(wx * wx + wy * wy + wz * wz);

    delphi = (2.0 * ran3() - 1.0) * delphi_max;
    cosphi = cos(delphi);
    sinphi = sin(delphi);

    rx = ux + cosphi * vx + sinphi * vmag * wx / wmag;
    ry = uy + cosphi * vy + sinphi * vmag * wy / wmag;
    rz = uz + cosphi * vz + sinphi * vmag * wz / wmag;

    r2x[k] = r2x[k - 1] + rx;
    r2y[k] = r2y[k - 1] + ry;
    r2z[k] = r2z[k - 1] + rz;
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
