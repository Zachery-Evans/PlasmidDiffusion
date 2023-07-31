#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <curand.h>
#include <curand_kernel.h>

#define PI 3.14159265359
#define MAX_MONOMERS 5000
#define MAX_INT 4294967295

const int N = 500;

long icyc, ncyc = 1e7, eq_cyc = 1e5;

void crankmethod(void);
void init_pos(void);
double ran3(curandGenerator_t);
__global__ void overlap(double[], double[], double[]);
__global__ void vectorMagnitude(double, double, double);
__global__ void crossProduct(double[], double[], double[]);
__global__ void position_check(double[], double[], double[], long);
__global__ void dotProduct(double, double, double, double, double, double);

__device__ int k;
__device__ double *rt2dev;
__device__ double *vecMag, *dotProd;
__device__ long posCheckTrue, overlapCheckTrue;
double x[MAX_MONOMERS], y[MAX_MONOMERS], z[MAX_MONOMERS];
double xOld[MAX_MONOMERS], yOld[MAX_MONOMERS], zOld[MAX_MONOMERS];
__device__ double xdev[MAX_MONOMERS], ydev[MAX_MONOMERS], zdev[MAX_MONOMERS];
__device__ double xdevOld[MAX_MONOMERS], ydevOld[MAX_MONOMERS], zdevOld[MAX_MONOMERS];
double *delrVec, *uVec, *vVec, *wVec, *uHat, *vHat, *wHat, *vprimeVec, *delrprimeVec; // x = [0] y = [1], z = [2]
double Rgsq_avg, Rgsq, dsrgRun = 0, acptRun = 0, acptRatio = 0, phi, delx, dely, delz;
double *vVec_mag, *delrVec_mag, *uHat_dot_delrVec, *uVec_mag, *delrVec_dot_uHat, dist_tot, xcm, ycm, zcm;

double rand_dev;
curandGenerator_t rgen;

int main(void)
{
    curandCreateGeneratorHost(&rgen, CURAND_RNG_PSEUDO_MTGP32);
    curandSetPseudoRandomGeneratorSeed(rgen, 138679324);
    rand_dev = (long)malloc(sizeof(long));

    for (int i = 0; i < N; i++)
    {
        printf("%lf \n", ran3(rgen));
    }

    for (icyc = 0; icyc < ncyc; icyc++) // For some predetermined number of cycles
    {
        for (int j = 0; j < N; j++) // For every particle on chain
        {

            k = N * rand_dev; // Select **random particle**. Since k is int, will change any 0 < #.0 < 1.0 -> #.0

            phi = PI * (2.0 * rand_dev - 1.0); // Random angle from 0 -> pi using similar method as (2.0 * ran3() - 1.0)

            if (k == 0) // DO NOT TOUCH !!!
            {
                delrVec[0] = xdev[k + 1] - xdev[k];
                delrVec[1] = ydev[k + 1] - ydev[k];
                delrVec[2] = zdev[k + 1] - zdev[k];

                uVec[0] = xdev[k + 1] - xdev[k + 2];
                uVec[1] = ydev[k + 1] - ydev[k + 2];
                uVec[2] = zdev[k + 1] - zdev[k + 2];
            }

            else if (k == N) // DO NOT TOUCH !!!
            {
                delrVec[0] = xdev[k] - xdev[k - 1];
                delrVec[1] = ydev[k] - ydev[k - 1];
                delrVec[2] = zdev[k] - zdev[k - 1];

                uVec[0] = xdev[k - 1] - xdev[k - 2];
                uVec[1] = ydev[k - 1] - ydev[k - 2];
                uVec[2] = zdev[k - 1] - zdev[k - 2];
            }

            else // DO NOT TOUCH !!!
            {
                delrVec[0] = xdev[k - 1] - xdev[k]; // delrVec is vector from k-1 to k
                delrVec[1] = ydev[k - 1] - ydev[k];
                delrVec[2] = zdev[k - 1] - zdev[k];

                uVec[0] = xdev[k + 1] - xdev[k - 1]; // uVec is vector from k+1 to k-1
                uVec[1] = ydev[k + 1] - ydev[k - 1];
                uVec[2] = zdev[k + 1] - zdev[k - 1];
            }

            vectorMagnitude<<<1, 1>>>(delrVec[0], delrVec[1], delrVec[2]);
            *delrVec_mag = *vecMag;
            // printf("%s %lf\n", "delrVec_mag:", delrVec_mag);

            vectorMagnitude<<<1, 1>>>(uVec[0], uVec[1], uVec[2]);
            *uVec_mag = *vecMag;

            uHat[0] = uVec[0] / *uVec_mag;
            uHat[1] = uVec[1] / *uVec_mag;
            uHat[2] = uVec[2] / *uVec_mag;

            dotProduct<<<1, 1>>>(uHat[0], uHat[1], uHat[2], delrVec[0], delrVec[1], delrVec[2]);
            *uHat_dot_delrVec = *dotProd;

            // printf("%lf\n", uHat_dot_delrVec);

            vVec[0] = delrVec[0] - uHat[0] * *uHat_dot_delrVec;
            vVec[1] = delrVec[1] - uHat[1] * *uHat_dot_delrVec;
            vVec[2] = delrVec[2] - uHat[2] * *uHat_dot_delrVec;

            vectorMagnitude<<<1, 1>>>(vVec[0], vVec[1], vVec[2]);
            // printf("%s %f\n", "vVec_mag:", vVec_mag); // Is zero sometimes somehow

            vHat[0] = vVec[0] / *vVec_mag;
            vHat[1] = vVec[1] / *vVec_mag;
            vHat[2] = vVec[2] / *vVec_mag;

            crossProduct<<<1, 1>>>(uHat, vHat, wHat); // Does a cross product of uVec and vVec, sets wVec parameters to x,y,z values

            wVec[0] = wHat[0] * *vVec_mag;
            wVec[1] = wHat[1] * *vVec_mag;
            wVec[2] = wHat[2] * *vVec_mag;

            vprimeVec[0] = vVec[0] * cos(phi) + wVec[0] * sin(phi);
            vprimeVec[1] = vVec[1] * cos(phi) + wVec[1] * sin(phi);
            vprimeVec[2] = vVec[2] * cos(phi) + wVec[2] * sin(phi);

            // printf("%lf %lf %lf\n", vprimeVec[0], vprimeVec[1], vprimeVec[2]);

            dotProduct<<<1, 1>>>(delrVec[0], delrVec[1], delrVec[2], uHat[0], uHat[1], uHat[2]);
            *delrVec_dot_uHat = *dotProd;

            delrprimeVec[0] = *delrVec_dot_uHat * uHat[0] + vprimeVec[0];
            delrprimeVec[1] = *delrVec_dot_uHat * uHat[1] + vprimeVec[1];
            delrprimeVec[2] = *delrVec_dot_uHat * uHat[2] + vprimeVec[2];

            delx = delrVec[0] - delrprimeVec[0]; // Random number from 0 -> 1 multiplied by maximum movement for each Cartesian coordinate
            dely = delrVec[1] - delrprimeVec[1];
            delz = delrVec[2] - delrprimeVec[2];

            xdevOld[k] = xdev[k]; // Save current position as old position
            ydevOld[k] = ydev[k];
            zdevOld[k] = zdev[k];

            xdev[k] += delx; // Change position of current particle
            ydev[k] += dely;
            zdev[k] += delz;

            position_check<<<1, 1>>>(xdev, ydev, zdev, k);
            
            if (posCheckTrue)
            { // If positon check returns 1, -> accept
                acptRun += 1;
                // printf("%d\n", k);
            }

            else
            { // If position check returns 0 -> disregard
                xdev[k] = xdevOld[k];
                ydev[k] = ydevOld[k];
                zdev[k] = zdevOld[k];
                dsrgRun += 1;
                // printf("%d\n", k);
            }

            // Reset position check and overlap check to null to preserve status.
            posCheckTrue = NULL;
            overlapCheckTrue = NULL;
        }

        Rgsq = 0.0;
        xcm = 0.0;
        ycm = 0.0;
        zcm = 0.0;

        for (int i = 0; i < N; i++)
        {
            xcm += xdev[i];
            ycm += ydev[i];
            zcm += zdev[i];
        }

        xcm /= N;
        ycm /= N;
        zcm /= N;

        for (int i = 0; i < N; i++)
        {
            Rgsq += (xdev[i] - xcm) * (xdev[i] - xcm) + (ydev[i] - ycm) * (ydev[i] - ycm) + (zdev[i] - zcm) * (zdev[i] - zcm);
        }

        Rgsq = Rgsq / N;

        if (icyc > eq_cyc)
        {
            Rgsq_avg += Rgsq;
        }
    }
}

double ran3(curandGenerator_t randomGenerator)
{
    curandGenerateUniformDouble(rgen, &rand_dev, 1);
    return rand_dev;
}

void init_pos()
{
    for (int i = 0; i < N; i++)
    { // Initial conditions of Polymer is centered on z axis
        y[i] = 0.0;
        if (i % 2 != 0)
        {
            x[i] = *rt2dev;
            z[i] = i * *rt2dev;
        }
        else
        {
            x[i] = 0.0;
            z[i] = i * *rt2dev;
        }
    }
}

__global__ void vectorMagnitude(double xvar, double yvar, double zvar)
{
    *vecMag = sqrt(xvar * xvar + yvar * yvar + zvar * zvar);
}

__global__ void dotProduct(double xvar_a, double yvar_a, double zvar_a, double xvar_b, double yvar_b, double zvar_b)
{
    *dotProd = xvar_a * xvar_b + yvar_a * yvar_b + zvar_a * zvar_b;
}

__global__ void crossProduct(double *vector1[3], double *vector2[3], double *result[3])
{
    *result[0] = (*vector1[1] * *vector2[2] - *vector1[2] * *vector2[1]);
    *result[1] = -(*vector1[0] * *vector2[2] - *vector1[2] * *vector2[0]);
    *result[2] = (*vector1[0] * *vector2[1] - *vector1[1] * *vector2[0]);
}

__global__ void position_check(double x_pos[], double y_pos[], double z_pos[], long k)
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
    if (posCheckTrue == NULL)
    {
        posCheckTrue = 1;
    }
    else
    {
        posCheckTrue = 0;
    }
}

// overlap checks if the particle overlaps with the one that came before it.
__global__ void overlap(double x_pos[], double y_pos[], double z_pos[])
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

    if (overlapCheckTrue == NULL)
    {
        overlapCheckTrue = 1;
    }
    else
    {
        overlapCheckTrue = 0;
    }
}
