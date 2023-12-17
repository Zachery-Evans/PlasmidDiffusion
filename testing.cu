#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h> // Include the cuda commands. Require these next few libraries
#include <curand_kernel.h> // Include the CUDA Random Number Generation library
#include <curand.h>

#define PI 3.14159265

void input(void);
void init_pos(double *, double *, double *, long);

__global__ void crank_dev(double *, double *, double *, long, long);
__global__ void randomKernel(double *, int);
__global__ void reptation_dev(double *, double *, double *, double, long, long, int, int *);
__global__ void reptation_dev_undo(double *, double *, double *, long, int);
__global__ void reptation_dev_shift(double *, double *, double *, long, int);
__global__ void check_accept(double *, double *, double *, long, long, int *);
__device__ void calc_delta_xyz(double *, double *, double *, long, int, double *, double *, double *, double, double);

__global__ void check_accept_reptation(double *, double *, double *, long, long, int *);

long nseg1, nseg2, nseg3, nseg4, i, j, k, ii, ncyc, nacc, kk, iseed;
long neq, ichain, nsamp, nacc_shift, nshift, xcmPrint, ycmPrint;
long imov, plasRigid, kmaxtest, freq_samp, cmFreqSamp, ngridx, ngridy;

long *monDev = nullptr;

double L, H, Ld2, Hd2, rmax, dx, dy, dz, dr2;
double drmin, drmax, gridspace, gridspacex_real, gridspacey_real;
double amax, bmin, amax2, bmin2, ecc, Area, rectangleArea, xBoxMax, yBoxMax, rshift_max;
double xBoxMaxd2, yBoxMaxd2;
double kappa, xold, yold, zold, delphi_max;
double **prob1, **plas, **probmon;

int host_check;

__device__ double *dev_Ld2 = nullptr, *dev_Hd2 = nullptr, *dev_rmax = nullptr, *dev_dx = nullptr, *dev_dy = nullptr, *dev_dz = nullptr, *dev_dr2 = nullptr;
__device__ double *dev_amax = nullptr, *dev_bmin = nullptr, *dev_amax2 = nullptr, *dev_bmin2 = nullptr, *dev_ecc = nullptr, *dev_Area = nullptr;
__device__ double *dev_rectangleArea = nullptr, *dev_xBoxMaxd2 = nullptr, *dev_yBoxMaxd2 = nullptr, *dev_zBoxMaxd2 = nullptr;
long *dev_N = nullptr;
int *dev_validity_check = nullptr;

int main()
{
	input();

	int accept = 0;
	long n = nseg1; // Number of elements to generate
	long threadsPerBlock, blocksPerGrid, mon;
	double *ran3Device = nullptr; // Device array to store random numbers
	double *ran3 = nullptr;		  // Host array to store random numbers
	double *xhost = nullptr, *yhost = nullptr, *zhost = nullptr;
	double *x = nullptr, *y = nullptr, *z = nullptr;
	double *dev_kappa = nullptr, rep_prob = 0.95;
	int *validity_check = nullptr, irep;

	cudaError_t cudaStatus;

	int imov = 1;

	FILE *fp;

	cudaStatus = cudaMalloc(&dev_N, sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc dev_N failed: %s\n", cudaGetErrorString(cudaStatus));
		return 1;
	}

	cudaStatus = cudaMalloc(&dev_kappa, sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc dev_kappa failed: %s\n", cudaGetErrorString(cudaStatus));
		return 1;
	}

	cudaStatus = cudaMalloc(&dev_validity_check, sizeof(int));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc dev_validity_check failed: %s\n", cudaGetErrorString(cudaStatus));
		return 1;
	}

	cudaStatus = cudaMalloc(&dev_rectangleArea, sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc dev_rectArea failed: %s\n", cudaGetErrorString(cudaStatus));
		return 1;
	}

	cudaStatus = cudaMalloc(&dev_Ld2, sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc dev_Ld2 failed: %s\n", cudaGetErrorString(cudaStatus));
		return 1;
	}

	cudaStatus = cudaMalloc(&dev_Hd2, sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc dev_Hd2 failed: %s\n", cudaGetErrorString(cudaStatus));
		return 1;
	}

	cudaStatus = cudaMalloc(&dev_xBoxMaxd2, sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc dev_xBoxMaxd2 failed: %s\n", cudaGetErrorString(cudaStatus));
		return 1;
	}

	cudaStatus = cudaMalloc(&dev_yBoxMaxd2, sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc dev_yBoxMaxd2 failed: %s\n", cudaGetErrorString(cudaStatus));
		return 1;
	}

		cudaStatus = cudaMalloc(&dev_zBoxMaxd2, sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc dev_zBoxMaxd2 failed: %s\n", cudaGetErrorString(cudaStatus));
		return 1;
	}

	bmin2 = bmin * bmin;

	if (ecc < 1.0)
	{
		amax2 = bmin2 / (1 - ecc * ecc);
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

	//  printf("%lf \t %lf\n", amax, bmin);
	//  printf("Length of the box: %lf\n", xBoxMax);
	//  printf("1/2 Length of the box: %lf\n", xBoxMaxd2);
	//  printf("Semi-major axis: %lf\n", amax);
	//  printf("Semi-minor axis: %lf\n", bmin);
	//  printf("Height of box: %lf\n", yBoxMax);
	Hd2 = H / 2.0;

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

	validity_check = (int *)malloc(sizeof(int));

	// Allocate memory on the host
	ran3 = (double *)malloc(n * sizeof(double));
	if (ran3 == nullptr)
	{
		printf("Memory allocation error on the host ran3.\n");
		return 1;
	}

	// Allocate memory on the device
	cudaStatus = cudaMalloc(&ran3Device, n * sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
		cudaFree(ran3Device);
		return 1;
	}

	cudaStatus = cudaMalloc(&x, n * sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc x failed: %s\n", cudaGetErrorString(cudaStatus));
		free(x);
		return 1;
	}

	cudaStatus = cudaMalloc(&y, n * sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc y failed: %s\n", cudaGetErrorString(cudaStatus));
		free(y);
		return 1;
	}

	cudaStatus = cudaMalloc(&z, n * sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc z failed: %s\n", cudaGetErrorString(cudaStatus));
		free(z);
		return 1;
	}
	/*
		cudaStatus = cudaMalloc(&monDev, sizeof(long));
		if (cudaStatus != cudaSuccess)
		{
			printf("cudaMalloc monDev failed: %s\n", cudaGetErrorString(cudaStatus));
			return 1;
		}
	*/
	init_pos(xhost, yhost, zhost, n);

	if (imov == 1)
	{
		fp = fopen("chain.xyz", "w");
		fprintf(fp, "%ld\n", n);
		fprintf(fp, "Polymer: %ld\n", 0);
		for (int i = 0; i < n; ++i)
		{
			fprintf(fp, "N    %lf  %lf  %lf\n", xhost[i], yhost[i], zhost[i]);
		}
	}

	threadsPerBlock = n;
	blocksPerGrid = 1;

	// Begin MC Cycles
	double *random_number = nullptr, *dev_random_number = nullptr;
	random_number = (double *)malloc(2 * sizeof(double));
	cudaStatus = cudaMalloc(&dev_random_number, 2 * sizeof(double));

	for (long ii = 0; ii < ncyc; ii++)
	{
		// Launch the kernel
		randomKernel<<<blocksPerGrid, threadsPerBlock>>>(ran3Device, n);

		// Copy the results back to the host
		cudaStatus = cudaMemcpy(ran3, ran3Device, n * sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			printf("cudaMemcpy ran3 failed: %s\n", cudaGetErrorString(cudaStatus));
			return 1;
		}

		// Check for kernel launch errors
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			printf("Line 200: Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return 1;
		}

		cudaMemcpy(x, xhost, n * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(y, yhost, n * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(z, zhost, n * sizeof(double), cudaMemcpyHostToDevice);

		for (int jj = 0; jj < n; jj++)
		{
			*validity_check = 0; // reset validity check for each MC cycle, since if > 0 then failure

			randomKernel<<<1, 2>>>(dev_random_number, 2);
			cudaMemcpy(random_number, dev_random_number, 2 * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(dev_validity_check, validity_check, sizeof(int), cudaMemcpyHostToDevice);

			mon = (long)n * ran3[jj]; // Choose random number out of total number of random numbers

			// printf("%ld\n", mon);
			//  printf("%ld\n", jj);

			kmaxtest = n - 2;

			if (random_number[0] >= rep_prob && (mon >= 2 || mon < kmaxtest))
			{
				crank_dev<<<1, 1>>>(x, y, z, n, mon);
				check_accept<<<blocksPerGrid, threadsPerBlock>>>(x, y, z, n, mon, dev_validity_check);
			}
			else
			{
				if (random_number[1] <= 0.5)
				{
					irep = 0;
					reptation_dev<<<1, 1>>>(x, y, z, kappa, nseg1, mon, irep, dev_validity_check);
				}
				else
				{
					irep = n - 1;
					reptation_dev<<<1, 1>>>(x, y, z, kappa, nseg1, mon, irep, dev_validity_check);
				}
			}

			cudaMemcpy(validity_check, dev_validity_check, sizeof(int), cudaMemcpyDeviceToHost);

			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess)
			{
				printf("dev_validity check failed memcpy: %s\n", cudaGetErrorString(cudaStatus));
				return 1;
			}

			if (*validity_check > accept)
			{
				cudaMemcpy(x, xhost, n * sizeof(double), cudaMemcpyHostToDevice); // Reset to configuration before MC cycle
				cudaMemcpy(y, yhost, n * sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(z, zhost, n * sizeof(double), cudaMemcpyHostToDevice);
			}
		}

		cudaMemcpy(xhost, x, n * sizeof(double), cudaMemcpyDeviceToHost); // Copy to new configuration
		cudaMemcpy(yhost, y, n * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(zhost, z, n * sizeof(double), cudaMemcpyDeviceToHost);

		if (imov == 1 && ii % freq_samp == 0)
		{
			fprintf(fp, "%ld\n", n);
			fprintf(fp, "Polymer: %ld\n", ii);
			for (int i = 0; i < n; ++i)
			{
				fprintf(fp, "N    %lf  %lf  %lf\n", xhost[i], yhost[i], zhost[i]);
			}
		}
	}

	if (imov == 1)
	{
		fclose(fp);
	}

	return 0;
}

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

void init_pos(double *xout, double *yout, double *zout, long N)
{
	double rt2 = sqrt(2) / 2;

	for (int i = 0; i < N; i++)
	{
		yout[i] = 0.000000000;

		if (i % 2 != 0)
		{
			xout[i] = rt2;
			zout[i] = i * rt2;
		}
		else
		{
			xout[i] = 0.0;
			zout[i] = i * rt2;
		}
	}
}

__global__ void crank_dev(double *xout, double *yout, double *zout, long n, long mon)
{
	double delphi_max = 3.141592653;
	double *randomNumber = nullptr;
	double xold, yold, zold;
	cudaMalloc(&randomNumber, sizeof(double));
	randomKernelDevice(randomNumber, 1);
	long tid = blockDim.x * blockIdx.x + threadIdx.x;

	if (tid < n)
	{
		double rx, ry, rz, Rx, Ry, Rz, Rmag, rdotRn, Rnx, Rny, Rnz;
		double ux, uy, uz, vx, vy, vz, vmag, wx, wy, wz, wmag;
		double cosphi, sinphi, delphi;

		xold = xout[tid];
		yold = yout[tid];
		zold = zout[tid];

		rx = xout[mon] - xout[mon - 1];
		ry = yout[mon] - yout[mon - 1];
		rz = zout[mon] - zout[mon - 1];

		Rx = xout[mon + 1] - xout[mon - 1];
		Ry = yout[mon + 1] - yout[mon - 1];
		Rz = zout[mon + 1] - zout[mon - 1];

		Rmag = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);

		Rnx = Rx / Rmag;
		Rny = Ry / Rmag;
		Rnz = Rz / Rmag;

		rx = rx;
		ry = ry;
		rz = rz;

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
			delphi = (randomNumber[0]) * delphi_max;
			cosphi = cos(delphi);
			sinphi = sin(delphi);

			rx = ux + cosphi * vx + sinphi * vmag * wx / wmag;
			ry = uy + cosphi * vy + sinphi * vmag * wy / wmag;
			rz = uz + cosphi * vz + sinphi * vmag * wz / wmag;

			xout[mon] = xout[mon - 1] + rx;
			yout[mon] = yout[mon - 1] + ry;
			zout[mon] = zout[mon - 1] + rz;
		}
		else
		{ // bonds are parallel
			xout[mon] = xold;
			yout[mon] = yold;
			zout[mon] = zold;
		}
		__syncthreads();
	}
	cudaFree(randomNumber);
}

__global__ void reptation_dev(double *r1x, double *r1y, double *r1z, double rept_kappa, long nseg, long ind, int irep_dev, int *check)
{
	double *rannum = nullptr, xold, yold, zold;
	double phi_prime, costheta_prime, *dx_fixed = nullptr, *dy_fixed = nullptr, *dz_fixed = nullptr;
	long tid = blockDim.x * blockIdx.x + threadIdx.x;

	if (tid < nseg)
	{
		cudaMalloc(&rannum, sizeof(double));
		cudaMalloc(&dx_fixed, sizeof(double));
		cudaMalloc(&dy_fixed, sizeof(double));
		cudaMalloc(&dz_fixed, sizeof(double));
		randomKernelDevice(rannum, 1);
		phi_prime = rannum[0] * 2.0 * 3.141592653;
		randomKernelDevice(rannum, 1);

		if (rept_kappa > -0.000001 && rept_kappa < 0.000001)
		{
			costheta_prime = (2.0 * rannum[0]) - 1.0;
		}
		else
		{
			costheta_prime = log((rannum[0] * exp(rept_kappa)) + ((1.0 - rannum[0]) * exp(-1.0 * rept_kappa))) / rept_kappa;
		}
		cudaFree(rannum);

		xold = r1x[irep_dev];
		yold = r1y[irep_dev];
		zold = r1z[irep_dev];

		// printf("%ld\n", irep_dev);
		calc_delta_xyz(r1x, r1y, r1z, nseg, irep_dev, dx_fixed, dy_fixed, dz_fixed, costheta_prime, phi_prime);

		if (irep_dev == 0)
		{

			for (ind = 0; ind < nseg - 1; ind++)
			{
				r1x[ind] = r1x[ind + 1];
				r1y[ind] = r1y[ind + 1];
				r1z[ind] = r1z[ind + 1];
			}

			// reptation_dev_shift<<<1, nseg>>>(r1x, r1y, r1z, nseg, irep_dev);
			//   printf("here\n");

			r1x[nseg - 1] = r1x[nseg - 2] + *dx_fixed;
			r1y[nseg - 1] = r1y[nseg - 2] + *dy_fixed;
			r1z[nseg - 1] = r1z[nseg - 2] + *dz_fixed;

			check_accept_reptation<<<1, nseg>>>(r1x, r1y, r1z, nseg, nseg - 1, check);

			if (*check != 0)
			{
				reptation_dev_undo<<<1, nseg>>>(r1x, r1y, r1z, nseg, irep_dev);

				r1x[0] = xold;
				r1y[0] = yold;
				r1z[0] = zold;
			}
		}
		else
		{ // irep == nseg1-1

			for (ind = nseg - 1; ind > 0; ind--)
			{
				r1x[ind] = r1x[ind - 1];
				r1y[ind] = r1y[ind - 1];
				r1z[ind] = r1z[ind - 1];
			}

			// reptation_dev_shift<<<1, nseg>>>(r1x, r1y, r1z, nseg, irep_dev);
			// printf("here\n");

			r1x[0] = r1x[1] + *dx_fixed;
			r1y[0] = r1y[1] + *dy_fixed;
			r1z[0] = r1z[1] + *dz_fixed;

			check_accept_reptation<<<1, nseg>>>(r1x, r1y, r1z, nseg, 0, check);

			if (*check != 0)
			{
				reptation_dev_undo<<<1, nseg>>>(r1x, r1y, r1z, nseg, irep_dev);

				r1x[nseg - 1] = xold;
				r1y[nseg - 1] = yold;
				r1z[nseg - 1] = zold;
			}
		}
	}
	cudaFree(dx_fixed);
	cudaFree(dy_fixed);
	cudaFree(dz_fixed);
	__syncthreads();
}

__global__ void reptation_dev_shift(double *r1x, double *r1y, double *r1z, long nseg, int irep)
{
	long tid = blockDim.x * blockIdx.x + threadIdx.x;

	if (tid < nseg)
	{
		if (irep == 0 && tid < nseg - 1)
		{
			r1x[tid] = r1x[tid + 1];
			r1y[tid] = r1y[tid + 1];
			r1z[tid] = r1z[tid + 1];
		}
		else if (irep != 0 && tid > 0)
		{
			r1x[tid] = r1x[tid - 1];
			r1y[tid] = r1y[tid - 1];
			r1z[tid] = r1z[tid - 1];
		}
	}
	__syncthreads();
}

__global__ void reptation_dev_undo(double *r1x, double *r1y, double *r1z, long nseg, int irep)
{
	long tid = blockDim.x * blockIdx.x + threadIdx.x;

	if (tid < nseg)
	{
		if (irep == 0 && tid < nseg - 1)
		{
			r1x[tid] = r1x[tid - 1];
			r1y[tid] = r1y[tid - 1];
			r1z[tid] = r1z[tid - 1];
		}
		else if (irep != 0 && tid > 0)
		{
			r1x[tid] = r1x[tid + 1];
			r1y[tid] = r1y[tid + 1];
			r1z[tid] = r1z[tid + 1];
		}
	}
	__syncthreads();
}

__device__ void calc_delta_xyz(double *r1x, double *r1y, double *r1z, long n_xyz, int irep, double *dx_fixed, double *dy_fixed, double *dz_fixed, double costheta_prime, double phi_prime)
{
	double dx_prime, dy_prime, dz_prime, dr_prime = 1.0, ux, uy, uz;
	double u, uxy, cos_beta, sin_beta, cos_alpha, sin_alpha, alpha;

	dx_prime = dr_prime * sqrt(1.0 - (costheta_prime * costheta_prime)) * cos(phi_prime);
	dy_prime = dr_prime * sqrt(1.0 - (costheta_prime * costheta_prime)) * sin(phi_prime);
	dz_prime = dr_prime * costheta_prime;

	if (irep == 0)
	{
		ux = r1x[n_xyz - 1] - r1x[n_xyz - 2];
		uy = r1y[n_xyz - 1] - r1y[n_xyz - 2];
		uz = r1z[n_xyz - 1] - r1z[n_xyz - 2];
	}
	else
	{
		ux = r1x[0] - r1x[1];
		uy = r1y[0] - r1y[1];
		uz = r1z[0] - r1z[1];
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
				alpha += 3.141592653;
			else
				alpha = (2.0 * 3.141592653) - alpha;
			cos_alpha = cos(alpha);
			sin_alpha = sin(alpha);
		}
	}

	// inverted matrix
	*dx_fixed = (cos_beta * cos_alpha * dx_prime) - (sin_alpha * dy_prime) + (sin_beta * cos_alpha * dz_prime);
	*dy_fixed = (cos_beta * sin_alpha * dx_prime) + (cos_alpha * dy_prime) + (sin_beta * sin_alpha * dz_prime);
	*dz_fixed = (cos_beta * dz_prime) - (sin_beta * dx_prime);
}

__global__ void check_accept_reptation(double *rx, double *ry, double *rz, long nseg, long krep, int *check)
{

	int reject = 1; // will return either accept or reject at end of function
	double dx, dy, dz, dr2;
	long tid = blockDim.x * blockIdx.x + threadIdx.x;

	if (tid < nseg)
	{
		// Check to see if iterative constant is greater than the size of all
		// polymers, if so, break inner loop and continue to next monomer in
		// checked polymer.
		/*
		if (squareEllipse(rx[kk], ry[kk], rz[kk]) == reject && kk < nseg)
		{
			return (reject);
		}
		*/

		if ((tid < krep - 1 || tid > krep + 1) && tid < nseg)
		{
			dz = rz[krep] - rz[tid];
			if (fabs(dz) < 1.0)
			{
				dx = rx[krep] - rx[tid];
				dy = ry[krep] - ry[tid];
				dr2 = dx * dx + dy * dy + dz * dz;
				if (dr2 < 1.0)
				{
					*check += reject; // if overlap with other monomer within chain, reject
				}
			}
		}
	}
}

/*
 *
 * checkEllipse is the function that, if it is determined the check is required to take into account the ellipse, checks whether or not
 * the polymer is overlapping with the geometry.
 *
 */
__device__ int checkEllipse(double xPos, double yPos, double zPos)
{
	int reject = 1, accept = 0;
	double echeck;

	if (zPos > *dev_Hd2 || zPos < -*dev_Hd2) // If z position is greater than flat surface of container, reject move
	{
		return reject;
	}

	if (*dev_ecc >= 1.0)
	{
		if (xPos > *dev_xBoxMaxd2 || xPos < -*dev_xBoxMaxd2)
		{
			return (reject);
		}
	}

	if (xPos > *dev_xBoxMaxd2 && *dev_ecc < 1.0)
	{ // If the polymer is outside of the leftmost semi-ellipse, reject
		if ((xPos - *dev_xBoxMaxd2) * (xPos - *dev_xBoxMaxd2) > *dev_amax2 * (1 - (yPos * yPos) / *dev_bmin2) && *dev_rectangleArea > 0.0)
		{
			return reject;
		}

		if (yPos * yPos > *dev_bmin2 * (1 - (xPos - *dev_xBoxMaxd2) * (xPos - *dev_xBoxMaxd2) / *dev_amax2) && *dev_rectangleArea > 0.0)
		{
			return reject;
		}

		if (*dev_rectangleArea > 0.0)
		{
			echeck = (((xPos - *dev_xBoxMaxd2) * (xPos - *dev_xBoxMaxd2)) / *dev_amax2) + ((yPos * yPos) / *dev_bmin2);
		}
		else
		{
			echeck = (((xPos) * (xPos)) / *dev_amax2) + ((yPos * yPos) / *dev_bmin2);
		}

		if (echeck > 1.0)
		{
			return (reject);
		}
	}

	else if (xPos < -*dev_xBoxMaxd2 && *dev_ecc < 1.0)
	{ // Checking if outside of left elliptical end
		if ((xPos + *dev_xBoxMaxd2) * (xPos + *dev_xBoxMaxd2) > *dev_amax2 * (1 - (yPos * yPos) / *dev_bmin2) && *dev_rectangleArea > 0.0)
		{
			return reject;
		}

		if (yPos * yPos > *dev_bmin2 * (1 - (xPos + *dev_xBoxMaxd2) * (xPos + *dev_xBoxMaxd2) / *dev_amax2) && *dev_rectangleArea > 0.0)
		{
			return reject;
		}

		if (*dev_rectangleArea > 0.0)
		{
			echeck = (((xPos + *dev_xBoxMaxd2) * (xPos + *dev_xBoxMaxd2)) / *dev_amax2) + ((yPos * yPos) / *dev_bmin2);
		}
		else
		{
			echeck = (((xPos) * (xPos)) / *dev_amax2) + ((yPos * yPos) / *dev_bmin2);
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
__device__ int squareEllipse(double xPos, double yPos, double zPos)
{
	int reject = 1, accept = 0;

	if (zPos < -*dev_Hd2 || zPos > *dev_Hd2) // Check if outside of the flat surface of the container
	{
		return (reject);
	}

	if (xPos < -*dev_xBoxMaxd2 || xPos > *dev_xBoxMaxd2) // If monomer is outside of the rectangle, check if outside of ellipse
	{
		return checkEllipse(xPos, yPos, zPos);
	}

	if (yPos > *dev_yBoxMaxd2 || yPos < -*dev_yBoxMaxd2)
	{
		return (reject);
	}

	return accept;
}

__global__ void check_accept(double *rx, double *ry, double *rz, long N, long monomer, int *check)
{
	long tid = blockIdx.x * blockDim.x + threadIdx.x;

	int reject = 1;
	double dx, dy, dz, dr2;

	if (tid < N)
	{
		/* Not restricting to a geometry for now.
		if (squareEllipse(rx[tid], ry[tid], rz[tid]) == reject)
		{
			*check = (reject);
		}
		*/

		if (tid < monomer - 1 || tid > monomer + 1)
		{
			// double testx = rx[monomer] - rx[tid];
			// double testy = ry[monomer] - ry[tid];
			// double testz = rz[monomer] - rz[tid];
			// double testdiff = testx * testx + testy * testy + testz * testz;
			// printf("Mon: %ld\tDistsq: %lf\tx: %lf\t y: %lf\t z: %lf\n", monomer, testdiff, testx, testy, testz);

			dx = rx[monomer] - rx[tid];
			dy = ry[monomer] - ry[tid];
			dz = rz[monomer] - rz[tid];
			dr2 = sqrt(dx * dx + dy * dy + dz * dz);

			if (dr2 < 1.0)
			{
				*check += reject;
			}
		}
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
// Copies the relevant input variables to the device
// ----------------------------------------------------------------------
/*
__host__ void input_dev(void)
{
	// dev_Ld2, dev_Hd2, dev_rmax, dev_dx, dev_dy, dev_dz, dev_dr2;
	// dev_amax, dev_bmin, dev_amax2, dev_bmin2, dev_ecc, dev_Area;
	// dev_rectangleArea, dev_xBoxMaxd2, dev_yBoxMaxd2;
	cudaError_t cudaStatus;
	cudaMemcpy(&rectangleArea, dev_rectangleArea, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&Ld2, dev_Ld2, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&Hd2, dev_Hd2, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&xBoxMaxd2, dev_xBoxMaxd2, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&yBoxMaxd2, dev_yBoxMaxd2, sizeof(double), cudaMemcpyDeviceToHost);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		printf("Input_dev Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}
*/