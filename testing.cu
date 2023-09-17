#include <stdio.h>
#include <curand_kernel.h> // Include the CUDA Random Number Generation library
#include <curand.h>

// Device function to generate a random number using curand
__device__ double randomNumber(curandState *state)
{
	// Generate a random double in the range [0, 1)
	return curand_uniform_double(state);
}

// CUDA kernel that uses the random number generator
__global__ void randomKernel(double *output, int n)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid < n)
	{
		// Initialize the random number generator state for each thread
		curandState state;
		curand_init(clock64(), tid, 0, &state);

		// Generate a random number for each thread and store it in the output array
		output[tid] = randomNumber(&state);
	}
}

__device__ void randomKernelDevice(double *output, int n)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid < n)
	{
		// Initialize the random number generator state for each thread
		curandState state;
		curand_init(clock64(), tid, 0, &state);

		// Generate a random number for each thread and store it in the output array
		output[tid] = randomNumber(&state);
	}
}

void init_pos_dev(double *xout, double *yout, double *zout, long N)
{
	for (int i = 0; i < N; i++)
	{
		yout[i] = 0.0;

		if (i % 2 != 0)
		{
			xout[i] = 0.70710678; // This is sqrt(2), not sure how to put this into device variable yet.
			zout[i] = i * 0.70710678;
		}
		else
		{
			xout[i] = 0.0;
			zout[i] = i * 0.70710678;
		}
	}
}

void crank_dev(double *xout, double *yout, double *zout, long k, long N)
{
	// k is the random number, N is the number of monomers in the polymer chain
	double xold = 0.0, yold = 0.0, zold = 0.0, delrVec[3], uVec[3], uHat[3], uVec_mag = 0.0, delrVec_mag = 0.0;
	double vVec[3], vprimeVec[3], vVec_mag, vHat[3], wVec[3], wVec_mag[3], wHat[3], delrprimeVec[3];
	double delx, dely, delz, uHat_dot_delrVec, delrVec_dot_uHat, phi;
	// randomKernel<<<1,1>>>(phi, 1);
	phi = 3.141592653;

	if (k == 0) // DO NOT TOUCH !!!
	{
		delrVec[0] = xout[k + 1] - xout[k];
		delrVec[1] = yout[k + 1] - yout[k];
		delrVec[2] = zout[k + 1] - zout[k];

		uVec[0] = xout[k + 1] - xout[k + 2];
		uVec[1] = yout[k + 1] - yout[k + 2];
		uVec[2] = zout[k + 1] - zout[k + 2];
	}

	else if (k == N) // DO NOT TOUCH !!!
	{
		delrVec[0] = xout[k] - xout[k - 1];
		delrVec[1] = yout[k] - yout[k - 1];
		delrVec[2] = zout[k] - zout[k - 1];

		uVec[0] = xout[k - 1] - xout[k - 2];
		uVec[1] = yout[k - 1] - yout[k - 2];
		uVec[2] = zout[k - 1] - zout[k - 2];
	}

	else // DO NOT TOUCH !!!
	{
		delrVec[0] = xout[k - 1] - xout[k]; // delrVec is vector from k-1 to k
		delrVec[1] = yout[k - 1] - yout[k];
		delrVec[2] = zout[k - 1] - zout[k];

		uVec[0] = xout[k + 1] - xout[k - 1]; // uVec is vector from k+1 to k-1
		uVec[1] = yout[k + 1] - yout[k - 1];
		uVec[2] = zout[k + 1] - zout[k - 1];
	}

	delrVec_mag = delrVec[0] * delrVec[0] + delrVec[1] * delrVec[1] + delrVec[2] * delrVec[2];
	// printf("%s %lf\n", "delrVec_mag:", delrVec_mag);

	uVec_mag = uVec[0] * uVec[0] + uVec[1] * uVec[1] + uVec[2] * uVec[2];

	uHat[0] = uVec[0] / uVec_mag;
	uHat[1] = uVec[1] / uVec_mag;
	uHat[2] = uVec[2] / uVec_mag;

	uHat_dot_delrVec = uHat[0] * delrVec[0] + uHat[1] * delrVec[1] + uHat[2] * delrVec[2];

	// printf("%lf\n", uHat_dot_delrVec);

	vVec[0] = delrVec[0] - uHat[0] * uHat_dot_delrVec;
	vVec[1] = delrVec[1] - uHat[1] * uHat_dot_delrVec;
	vVec[2] = delrVec[2] - uHat[2] * uHat_dot_delrVec;

	vVec_mag = vVec[0] * vVec[0] + vVec[1] * vVec[1] + vVec[2] * vVec[2];
	// printf("%s %f\n", "vVec_mag:", vVec_mag); // Is zero sometimes somehow

	vHat[0] = vVec[0] / vVec_mag;
	vHat[1] = vVec[1] / vVec_mag;
	vHat[2] = vVec[2] / vVec_mag;

	wHat[0] = uHat[1] * vHat[2] - uHat[2] * vHat[1]; // Does a cross product of uVec and vVec, sets wVec parameters to x,y,z values
	wHat[1] = uHat[0] * vHat[2] - uHat[2] * vHat[0];
	wHat[2] = uHat[0] * vHat[1] - uHat[1] * vHat[0];

	wVec[0] = wHat[0] * vVec_mag;
	wVec[1] = wHat[1] * vVec_mag;
	wVec[2] = wHat[2] * vVec_mag;

	vprimeVec[0] = vVec[0] * cos(phi) + wVec[0] * sin(phi);
	vprimeVec[1] = vVec[1] * cos(phi) + wVec[1] * sin(phi);
	vprimeVec[2] = vVec[2] * cos(phi) + wVec[2] * sin(phi);

	// printf("%lf %lf %lf\n", vprimeVec[0], vprimeVec[1], vprimeVec[2]);

	delrVec_dot_uHat = delrVec[0] * uHat[0] + delrVec[1] * uHat[1] + delrVec[2] * uHat[2];

	delrprimeVec[0] = delrVec_dot_uHat * uHat[0] + vprimeVec[0];
	delrprimeVec[1] = delrVec_dot_uHat * uHat[1] + vprimeVec[1];
	delrprimeVec[2] = delrVec_dot_uHat * uHat[2] + vprimeVec[2];

	delx = delrVec[0] - delrprimeVec[0]; // Random number from 0 -> 1 multiplied by maximum movement for each Cartesian coordinate
	dely = delrVec[1] - delrprimeVec[1];
	delz = delrVec[2] - delrprimeVec[2];

	xold = xout[k]; // Save current position as old position
	yold = yout[k];
	zold = zout[k];

	xout[k] += delx; // Change position of current particle
	yout[k] += dely;
	zout[k] += delz;
}

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

int check_accept(double rx, double ry, double rz, long nseg)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	int accept, reject;
	long klow, khigh;
	double dx, dy, dz, dr, dr2;

	accept = 0;
	reject = 1;

	if (tid < n)
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

	return accept; // apply rigidity
}

// overlap checks if the particle overlaps with the one that came before it.
__global__ void overlap(double x_pos[], double y_pos[], double z_pos[], long k, long N)
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
			overlapCheck = 0;
		}
	}

	if (overlapCheck == NULL)
	{
		overlapCheck = 1;
	}
	else
	{
		overlapCheck = 0;
	}
}

int posCheck, overlapCheck;

int main()
{
	long n = 10;				  // Number of elements to generate
	double *ran3Device = nullptr; // Device array to store random numbers
	double *ran3 = nullptr;		  // Host array to store random numbers
	double *xhost = nullptr, *yhost = nullptr, *zhost = nullptr;
	double *x = nullptr, *y = nullptr, *z = nullptr;
	long threadsPerBlock, blocksPerGrid, ncyc = 10, freq_samp = 1;
	long mon;

	int imov = 1;

	FILE *fp;

	xhost = (double *)malloc(n * sizeof(double));
	if (xhost == nullptr)
	{
		printf("Memory allocation error on the host xhost.\n");
		return 1;
	}

	yhost = (double *)malloc(n * sizeof(double));
	if (yhost == nullptr)
	{
		printf("Memory allocation error on the host yhost.\n");
		return 1;
	}

	zhost = (double *)malloc(n * sizeof(double));
	if (zhost == nullptr)
	{
		printf("Memory allocation error on the host zhost.\n");
		return 1;
	}

	// Allocate memory on the host
	ran3 = (double *)malloc(n * sizeof(double));
	if (ran3 == nullptr)
	{
		printf("Memory allocation error on the host ran3.\n");
		return 1;
	}

	// Allocate memory on the device
	cudaError_t cudaStatus = cudaMalloc((void **)&ran3Device, n * sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
		free(ran3);
		return 1;
	}

	cudaStatus = cudaMalloc((void **)&x, n * sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc x failed: %s\n", cudaGetErrorString(cudaStatus));
		free(ran3);
		return 1;
	}

	cudaStatus = cudaMalloc((void **)&y, n * sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc y failed: %s\n", cudaGetErrorString(cudaStatus));
		free(ran3);
		return 1;
	}

	cudaStatus = cudaMalloc((void **)&z, n * sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc z failed: %s\n", cudaGetErrorString(cudaStatus));
		free(ran3);
		return 1;
	}

	if (imov == 1)
	{
		fp = fopen("chain.xyz", "w");
		fprintf(fp, "%ld\n", n);
	}

	init_pos_dev(xhost, yhost, zhost, n);

	// Begin MC Cycles
	for (long ii = 0; ii < ncyc; ii++)
	{
		// Launch the kernel
		threadsPerBlock = 1;
		blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;
		randomKernel<<<blocksPerGrid, threadsPerBlock>>>(ran3Device, n);

		// Check for kernel launch errors
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			printf("Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			cudaFree(ran3Device);
			free(ran3);
			return 1;
		}

		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			printf("Kernal launch failed: %s\n", cudaGetErrorString(cudaStatus));
			cudaFree(x);
			cudaFree(y);
			cudaFree(z);
			free(xhost);
			free(yhost);
			free(zhost);
			return 1;
		}

		// Copy the results back to the host
		cudaStatus = cudaMemcpy(ran3, ran3Device, n * sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			printf("cudaMemcpy ran3 failed: %s\n", cudaGetErrorString(cudaStatus));
			cudaFree(ran3Device);
			free(ran3);
			return 1;
		}

		for (int jj = 0; jj < n; jj++)
		{
			mon = n * ran3[jj];
			printf("%ld\n", mon);
			// printf("%ld\n", jj);
			crank_dev(xhost, yhost, zhost, mon, n);

			// Print the generated random numbers on the host

			if (imov == 1 && ii % freq_samp == 0)
			{
				fprintf(fp, "Polymer: %ld\n", ii);
				for (int i = 0; i < n; ++i)
				{
					fprintf(fp, "%lf  %lf  %lf\n", xhost[i], yhost[i], zhost[i]);
				}
			}
		}
	}

	if (imov == 1)
	{
		fclose(fp);
	}

	// Free device and host memory
	cudaFree(ran3Device);
	free(ran3);

	return 0;
}
