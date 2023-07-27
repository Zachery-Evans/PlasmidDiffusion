#include "operators.h"
#include "constants.h"
#include <math.h>
#include <stdio.h>

double vectorMagnitude(double xvar, double yvar, double zvar)
{
    return sqrt(xvar * xvar + yvar * yvar + zvar * zvar);
}

double dotProduct(double xvar_a, double yvar_a, double zvar_a, double xvar_b, double yvar_b, double zvar_b)
{
    return xvar_a * xvar_b + yvar_a * yvar_b + zvar_a * zvar_b;
}

void crossProduct(double vector1[3], double vector2[3], double result[3])
{
    result[0] = (vector1[1] * vector2[2] - vector1[2] * vector2[1]);
    result[1] = -(vector1[0] * vector2[2] - vector1[2] * vector2[0]);
    result[2] = (vector1[0] * vector2[1] - vector1[1] * vector2[0]);
}

int position_check(double x_pos[], double y_pos[], double z_pos[])
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
            return 0;
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
            return 0;
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
                return 0;
            }
        }
    }
    return 1;
}

// overlap checks if the particle overlaps with the one that came before it.
int overlap(double x_pos[], double y_pos[], double z_pos[])
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
            return 0;
        }
    }
    return 1;
}
