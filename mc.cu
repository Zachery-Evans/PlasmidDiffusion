#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265359
#define MAX_MONOMERS 5000
const int N = 100;

double ran3(void);
__global__ void crankmethod();
__global__ void crossProduct(double vector1[], double vector2[], double result[]);
__global__ void position_check(double x_pos[], double y_pos[], double z_pos[]);
__global__ void overlap(double x_pos[], double y_pos[], double z_pos[]);
__global__ void vectorMagnitude(double xvar, double yvar, double zvar);
__global__ void dotProduct(double xvar_a, double yvar_a, double zvar_a, double xvar_b, double yvar_b, double zvar_b);

int k;
long SEEDF = 18394783;
double xold[MAX_MONOMERS], yold[MAX_MONOMERS], zold[MAX_MONOMERS], x[MAX_MONOMERS], y[MAX_MONOMERS], z[MAX_MONOMERS], r[MAX_MONOMERS];
double x[MAX_MONOMERS], y[MAX_MONOMERS], z[MAX_MONOMERS];
__global__ double xdev[MAX_MONOMERS], ydev[MAX_MONOMERS], zdev[MAX_MONOMERS];
__global__ double vecMag, dotProd;
__global__ long posCheckTrue;
__global__ long overlapCheckTrue;

int main(void)
{
    crankmethod();

    return 0;
}

void crankmethod(void)
{
    int icyc, ncyc = 5e5, eq_cyc = 1e4, i, j;
    double rt2 = sqrt(2) / 2, delrVec[3], uVec[3], vVec[3], wVec[3], uHat[3], vHat[3], wHat[3], vprimeVec[3], delrprimeVec[3]; // x = [0] y = [1], z = [2]
    double Rgsq_avg, Rgsq, delmax = 0.15, dsrgRun = 0, acptRun = 0, acptRatio = 0, sigma = 1.0, phi, delx, dely, delz;
    double dr_pos, vVec_mag, uHat_dot_delrVec, delrVec_mag, uVec_mag, delrVec_dot_uHat, dx_pos, dy_pos, dz_pos, dist_tot, xcm, ycm, zcm;

    FILE *gp; // Monomer Position Data
    if ((gp = fopen("mc.dat", "w")) == NULL)
    {
        printf("Cannot open file: mc.dat");
        exit(0);
    }

    else
    {
        for (i = 0; i < N; i++)
        { // Initial conditions of Polymer is centered on z axis
            y[i] = 0.0;
            if (i % 2 != 0)
            {
                x[i] = rt2;
                z[i] = i * rt2;
            }

            else
            {
                x[i] = 0.0;
                z[i] = i * rt2;
            }
        }

        for (icyc = 0; icyc < ncyc; icyc++) // For some predetermined number of cycles
        {
            for (j = 0; j < N; j++) // For every particle on chain
            {
                k = N * ran3(); // Select **random particle**. Since k is int, will change any 0 < #.0 < 1.0 -> #.0

                phi = PI * (2.0 * ran3() - 1.0); // Random angle from 0 -> pi using similar method as (2.0 * ran3() - 1.0)

                if (k == 0) // DO NOT TOUCH !!!
                {
                    delrVec[0] = x[k + 1] - x[k];
                    delrVec[1] = y[k + 1] - y[k];
                    delrVec[2] = z[k + 1] - z[k];

                    uVec[0] = x[k + 1] - x[k + 2];
                    uVec[1] = y[k + 1] - y[k + 2];
                    uVec[2] = z[k + 1] - z[k + 2];
                }

                else if (k == N) // DO NOT TOUCH !!!
                {
                    delrVec[0] = x[k] - x[k - 1];
                    delrVec[1] = y[k] - y[k - 1];
                    delrVec[2] = z[k] - z[k - 1];

                    uVec[0] = x[k - 1] - x[k - 2];
                    uVec[1] = y[k - 1] - y[k - 2];
                    uVec[2] = z[k - 1] - z[k - 2];
                }

                else // DO NOT TOUCH !!!
                {
                    delrVec[0] = x[k - 1] - x[k]; // delrVec is vector from k-1 to k
                    delrVec[1] = y[k - 1] - y[k];
                    delrVec[2] = z[k - 1] - z[k];

                    uVec[0] = x[k + 1] - x[k - 1]; // uVec is vector from k+1 to k-1
                    uVec[1] = y[k + 1] - y[k - 1];
                    uVec[2] = z[k + 1] - z[k - 1];
                }

                vectorMagnitude<<<1,1>>>(delrVec[0], delrVec[1], delrVec[2]);
                delrVec_mag = vecMag;
                // printf("%s %lf\n", "delrVec_mag:", delrVec_mag);

                vectorMagnitude<<<1,1>>>(uVec[0], uVec[1], uVec[2]);
                uVec_mag = vecMag;

                uHat[0] = uVec[0] / uVec_mag;
                uHat[1] = uVec[1] / uVec_mag;
                uHat[2] = uVec[2] / uVec_mag;

                dotProduct<<<1,1>>>(uHat[0], uHat[1], uHat[2], delrVec[0], delrVec[1], delrVec[2]);
                uHat_dot_delrVec = dotProd;

                // printf("%lf\n", uHat_dot_delrVec);

                vVec[0] = delrVec[0] - uHat[0] * uHat_dot_delrVec;
                vVec[1] = delrVec[1] - uHat[1] * uHat_dot_delrVec;
                vVec[2] = delrVec[2] - uHat[2] * uHat_dot_delrVec;

                vectorMagnitude<<<1,1>>>(vVec[0], vVec[1], vVec[2]);
                vVec_mag = vecMag;
                // printf("%s %f\n", "vVec_mag:", vVec_mag); // Is zero sometimes somehow

                vHat[0] = vVec[0] / vVec_mag;
                vHat[1] = vVec[1] / vVec_mag;
                vHat[2] = vVec[2] / vVec_mag;

                crossProduct(uHat, vHat, wHat); // Does a cross product of uVec and vVec, sets wVec parameters to x,y,z values

                wVec[0] = wHat[0] * vVec_mag;
                wVec[1] = wHat[1] * vVec_mag;
                wVec[2] = wHat[2] * vVec_mag;

                vprimeVec[0] = vVec[0] * cos(phi) + wVec[0] * sin(phi);
                vprimeVec[1] = vVec[1] * cos(phi) + wVec[1] * sin(phi);
                vprimeVec[2] = vVec[2] * cos(phi) + wVec[2] * sin(phi);

                // printf("%lf %lf %lf\n", vprimeVec[0], vprimeVec[1], vprimeVec[2]);

                dotProduct<<<1,1>>>(delrVec[0], delrVec[1], delrVec[2], uHat[0], uHat[1], uHat[2]);
                delrVec_dot_uHat = dotProd;

                delrprimeVec[0] = delrVec_dot_uHat * uHat[0] + vprimeVec[0];
                delrprimeVec[1] = delrVec_dot_uHat * uHat[1] + vprimeVec[1];
                delrprimeVec[2] = delrVec_dot_uHat * uHat[2] + vprimeVec[2];

                delx = delrVec[0] - delrprimeVec[0]; // Random number from 0 -> 1 multiplied by maximum movement for each Cartesian coordinate
                dely = delrVec[1] - delrprimeVec[1];
                delz = delrVec[2] - delrprimeVec[2];

                xold[k] = x[k]; // Save current position as old position
                yold[k] = y[k];
                zold[k] = z[k];

                x[k] += delx; // Change position of current particle
                y[k] += dely;
                z[k] += delz;

                position_check<<<1,1>>>(x, y, z);
                if (posCheckTrue)
                { // If positon check returns 1, -> accept
                    acptRun += 1;
                    // printf("%d\n", k);
                }

                else
                { // If position check returns 0 -> disregard
                    x[k] = xold[k];
                    y[k] = yold[k];
                    z[k] = zold[k];
                    dsrgRun += 1;
                    // printf("%d\n", k);
                }
            }

            Rgsq = 0.0;
            xcm = 0.0;
            ycm = 0.0;
            zcm = 0.0;

            for (i = 0; i < N; i++)
            {
                xcm += x[i];
                ycm += y[i];
                zcm += z[i];
            }

            xcm /= N;
            ycm /= N;
            zcm /= N;

            for (i = 0; i < N; i++)
            {
                Rgsq += (x[i] - xcm) * (x[i] - xcm) + (y[i] - ycm) * (y[i] - ycm) + (z[i] - zcm) * (z[i] - zcm);
            }

            Rgsq = Rgsq / N;

            if (icyc % 100 == 0)
            {
                fprintf(gp, "%d %lf\n", icyc, Rgsq);
            }

            if (icyc > eq_cyc)
            {
                Rgsq_avg += Rgsq;
            }
        }

        acptRatio = acptRun / (acptRun + dsrgRun);
        Rgsq_avg /= (ncyc - eq_cyc);

        printf("%s %lf\n", "\nAcceptance Ratio:", acptRatio);

        printf("\n%s %lf %s %lf\n", "Disregarded Number:", dsrgRun, "Accepted Number:", acptRun);

        printf("\n%s %lf\n", "Average Radius of Gyration:", Rgsq_avg);

        fclose(gp);
    }
}

__global__ 
void vectorMagnitude(double xvar, double yvar, double zvar)
{
    vecMag = sqrt(xvar * xvar + yvar * yvar + zvar * zvar);
}

__global__
void dotProduct(double xvar_a, double yvar_a, double zvar_a, double xvar_b, double yvar_b, double zvar_b)
{
    dotProd = xvar_a * xvar_b + yvar_a * yvar_b + zvar_a * zvar_b;
}

__global__
void crossProduct(double vector1[3], double vector2[3], double result[3])
{
    result[0] = (vector1[1] * vector2[2] - vector1[2] * vector2[1]);
    result[1] = -(vector1[0] * vector2[2] - vector1[2] * vector2[0]);
    result[2] = (vector1[0] * vector2[1] - vector1[1] * vector2[0]);
}

__global__
void position_check(double x_pos[], double y_pos[], double z_pos[])
{
    double dx_pos, dy_pos, dz_pos, dist_tot;

    // check bonded particles distance

    if (k != 0)
    {
        dx_pos = x_pos[k] - x_pos[k - 1];
        dy_pos = y_pos[k] - y_pos[k - 1];
        dz_pos = z_pos[k] - z_pos[k - 1];

        dist_tot = dx_pos * dx_pos + dy_pos * dy_pos + dz_pos * dz_pos;

        if (dist_tot < 0.99 || 1.01 < dist_tot) // 0.9*0.9 since dist_tot is not in sqrt() for speed
        {
            posCheckTrue = 0;
        }
    }

    if (k != N - 1)
    {
        dx_pos = x_pos[k + 1] - x_pos[k];
        dy_pos = y_pos[k + 1] - y_pos[k];
        dz_pos = z_pos[k + 1] - z_pos[k];

        dist_tot = dx_pos * dx_pos + dy_pos * dy_pos + dz_pos * dz_pos;

        // if dist w/n acceptable values
        if (dist_tot < 0.99 || 1.01 < dist_tot)
        {
            posCheckTrue = 0;
        }
    }

    // check non-bonded particles distance
    for (int i = 0; i < N - 1; i++) // NOTE: "N-1" instead of "N"
    {
        if (i != k)
        {
            dx_pos = x_pos[k] - x_pos[i]; // Measure x y z distance from particle a to particle b
            dy_pos = y_pos[k] - y_pos[i];
            dz_pos = z_pos[k] - z_pos[i];

            dist_tot = dx_pos * dx_pos + dy_pos * dy_pos + dz_pos * dz_pos; // Calculate magnitude of dist squared

            if (dist_tot < 1.0)
            {
                posCheckTrue = 0;
            }
        }
    }
    posCheckTrue = 1;
}

// overlap checks if the particle overlaps with the one that came before it.
__global__
void overlap(double x_pos[], double y_pos[], double z_pos[])
{
    double dx_pos, dy_pos, dz_pos, dist_tot;

    // check non-bonded particles distance
    for (int i = 0; i < k - 1; i++) // NOTE: "k-1" instead of "k"
    {
        dx_pos = x_pos[k] - x_pos[i]; // Measure x y z distance from particle a to particle b
        dy_pos = y_pos[k] - y_pos[i];
        dz_pos = z_pos[k] - z_pos[i];

        dist_tot = dx_pos * dx_pos + dy_pos * dy_pos + dz_pos * dz_pos; // Calculate magnitude of dist squared

        if (dist_tot < 1.0)
        {
            overlapCheckTrue = 0;
        }
    }
    overlapCheckTrue = 1;
}

double ran3()
{
    long iseed = SEEDF;
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
