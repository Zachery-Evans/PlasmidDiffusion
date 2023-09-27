#include <stdio.h>
#include <curand_kernel.h> // Include the CUDA Random Number Generation library
#include <curand.h>

#define PI 3.14159265

void input(void);
void init_pos(double *, double *, double *, long);
void crank_dev(double *, double *, double *, long, long);

__global__ void check_accept(double *, double *, double *, long, volatile bool *);
__global__ void randomKernel(double *, int);

long nseg1, nseg2, nseg3, nseg4, i, j, k, ii, ncyc, nacc, kk, iseed;
long neq, ichain, nsamp, nacc_shift, nshift, xcmPrint, ycmPrint;
long imov, plasRigid, kmaxtest, freq_samp, cmFreqSamp, ngridx, ngridy;

__device__ long *monDev = nullptr;

double L, H, Ld2, Hd2, rmax, dx, dy, dz, dr2;
double drmin, drmax, gridspace, gridspacex_real, gridspacey_real;
double amax, bmin, amax2, bmin2, ecc, Area, rectangleArea, xBoxMax, yBoxMax, rshift_max;
double xBoxMaxd2, yBoxMaxd2;
double kappa, xold, yold, zold, delphi_max;
double **prob1, **plas, **probmon;

int main()
{
	long n = 10; // Number of elements to generate
	long threadsPerBlock, blocksPerGrid, freq_samp = 1, *mon = nullptr;
	double *ran3Device = nullptr; // Device array to store random numbers
	double *ran3 = nullptr;		  // Host array to store random numbers
	double *xhost = nullptr, *yhost = nullptr, *zhost = nullptr;
	double *x = nullptr, *y = nullptr, *z = nullptr;

	bool *validityCheckHost = nullptr;

	int imov = 1;

	FILE *fp;

	input();

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

	validityCheckHost = (bool *)malloc(n * sizeof(bool));
	if (validityCheckHost == nullptr)
	{
		printf("Memory allocation error on the host validityCheckHost.\n");
		return 1;
	}

	mon = (long *)malloc(sizeof(long));
	if (mon == nullptr)
	{
		printf("Memory allocation error on the host validityCheckHost.\n");
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
		free(x);
		return 1;
	}

	cudaStatus = cudaMalloc((void **)&y, n * sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc y failed: %s\n", cudaGetErrorString(cudaStatus));
		free(y);
		return 1;
	}

	cudaStatus = cudaMalloc((void **)&z, n * sizeof(double));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc z failed: %s\n", cudaGetErrorString(cudaStatus));
		free(z);
		return 1;
	}

	cudaStatus = cudaMalloc((void **)&validityCheck, n * sizeof(bool));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc validityCheck failed: %s\n", cudaGetErrorString(cudaStatus));
		free(validityCheck);
		return 1;
	}

	cudaStatus = cudaMalloc((void **)&monDev, sizeof(long));
	if (cudaStatus != cudaSuccess)
	{
		printf("cudaMalloc monDev failed: %s\n", cudaGetErrorString(cudaStatus));
		cudaFree(monDev);
		return 1;
	}

	if (imov == 1)
	{
		fp = fopen("chain.xyz", "w");
		fprintf(fp, "%ld\n", n);
	}

	init_pos(xhost, yhost, zhost, n);

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
			volatile bool *validityCheck = false;
			*mon = (long)n * ran3[jj];
			// printf("%ld\n", mon);
			//  printf("%ld\n", jj);
			cudaMemcpy(monDev, mon, n * sizeof(double), cudaMemcpyHostToDevice);

			crank_dev(xhost, yhost, zhost, *monDev, n);

			cudaMemcpy(x, xhost, n * sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(y, yhost, n * sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(z, xhost, n * sizeof(double), cudaMemcpyHostToDevice);

			check_accept<<<threadsPerBlock, blocksPerGrid>>>(x, y, z, n, validityCheck);

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
 *
 * checkEllipse is the function that, if it is determined the check is required to take into account the ellipse, checks whether or not
 * the polymer is overlapping with the geometry.
 *
 */
__device__ int checkEllipse(double xPos, double yPos, double zPos)
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
__device__ int squareEllipse(double xPos, double yPos, double zPos)
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

__global__ void check_accept(double *rx, double *ry, double *rz, long n, volatile bool *check)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	int accept, reject;
	long klow, khigh;
	double dx, dy, dz, dr, dr2;

	accept = 0;
	reject = 1;

	if (tid < n)
	{
		/* Not restricting to a geometry for now.
		if (squareEllipse(rx[kk], ry[kk], rz[kk]) == reject)
		{
			return (reject);
		}
		*/
		if (tid < *monDev - 1 || tid > *monDev + 1)
		{
			dx = rx[*monDev] - rx[tid];
			dy = ry[*monDev] - ry[tid];
			dz = rz[*monDev] - rz[tid];
			dr2 = dx * dx + dy * dy + dz * dz;
			if (dr2 < 1.0)
			{
				*check = false;
			}
		}
	}

	*check = true;
}

// overlap checks if the particle overlaps with the one that came before it.
__global__ void overlap(double *x_pos, double *y_pos, double *z_pos, long k, long N, volatile bool *check)
{
	double dx_pos, dy_pos, dz_pos, dist_tot;
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (tid < N)
	{
		// check non-bonded particles distance

		if (k != tid)
		{
			dx_pos = x_pos[k] - x_pos[tid]; // Measure x y z distance from particle a to particle b
			dy_pos = y_pos[k] - y_pos[tid];
			dz_pos = z_pos[k] - z_pos[tid];

			dist_tot = dx_pos * dx_pos + dy_pos * dy_pos + dz_pos * dz_pos; // Calculate magnitude of dist squared

			if (dist_tot < 1.0)
			{
				check[tid] = false;
			}
		}

		check[tid] = true;
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