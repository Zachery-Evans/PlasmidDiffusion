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
	if (tid < n) {
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
	if (tid < n) {
		// Initialize the random number generator state for each thread
		curandState state;
		curand_init(clock64(), tid, 0, &state);

		// Generate a random number for each thread and store it in the output array
		output[tid] = randomNumber(&state);
		cudaDeviceSynchronize();
	}
}


__global__ void init_pos_dev(double *xout, double *yout, double *zout, long N)
{
	for(int i = 0; i < N; i++)
	{	
		yout[i] = 0.0;

		if(i % 2 != 0)
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

__global__ void crank_dev(double *x, double *y, double *z, long k, long N) {

	double xold=0.0, yold=0.0, zold=0.0, delrVec[3], uVec[3], uHat[3], uVec_mag=0.0,  delrVec_mag=0.0;
	double vVec[3], vprimeVec[3], vVec_mag, vHat[3], wVec[3], wVec_mag[3], wHat[3], delrprimeVec[3];
	double delx, dely, delz, uHat_dot_delrVec, delrVec_dot_uHat, *phi;

	__syncthreads();

	randomKernel<<<1,1>>>(phi, 1);
	*phi = *phi * 3.141592653;

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

	delrVec_mag = delrVec[0]*delrVec[0] + delrVec[1]*delrVec[1] + delrVec[2]*delrVec[2];
	// printf("%s %lf\n", "delrVec_mag:", delrVec_mag);

	uVec_mag = uVec[0] * uVec[0] + uVec[1] * uVec[1] + uVec[2] * uVec[2];

	uHat[0] = uVec[0] / uVec_mag;
	uHat[1] = uVec[1] / uVec_mag;
	uHat[2] = uVec[2] / uVec_mag;

	uHat_dot_delrVec = uHat[0]*delrVec[0] + uHat[1]*delrVec[1] + uHat[2]*delrVec[2];

	// printf("%lf\n", uHat_dot_delrVec);

	vVec[0] = delrVec[0] - uHat[0] * uHat_dot_delrVec;
	vVec[1] = delrVec[1] - uHat[1] * uHat_dot_delrVec;
	vVec[2] = delrVec[2] - uHat[2] * uHat_dot_delrVec;

	vVec_mag = vVec[0]*vVec[0] + vVec[1]*vVec[1] + vVec[2]*vVec[2];
	// printf("%s %f\n", "vVec_mag:", vVec_mag); // Is zero sometimes somehow

	vHat[0] = vVec[0] / vVec_mag;
	vHat[1] = vVec[1] / vVec_mag;
	vHat[2] = vVec[2] / vVec_mag;

	wHat[0] = uHat[1]*vHat[2] - uHat[2]*vHat[1];  // Does a cross product of uVec and vVec, sets wVec parameters to x,y,z values
	wHat[1] = uHat[0]*vHat[2] - uHat[2]*vHat[0];
	wHat[2] = uHat[0]*vHat[1] - uHat[1]*vHat[0];

	wVec[0] = wHat[0] * vVec_mag;
	wVec[1] = wHat[1] * vVec_mag;
	wVec[2] = wHat[2] * vVec_mag;

	vprimeVec[0] = vVec[0] * cos(*phi) + wVec[0] * sin(*phi);
	vprimeVec[1] = vVec[1] * cos(*phi) + wVec[1] * sin(*phi);
	vprimeVec[2] = vVec[2] * cos(*phi) + wVec[2] * sin(*phi);

	// printf("%lf %lf %lf\n", vprimeVec[0], vprimeVec[1], vprimeVec[2]);

	delrVec_dot_uHat = delrVec[0]*uHat[0] + delrVec[1]*uHat[1] + delrVec[2]*uHat[2];

	delrprimeVec[0] = delrVec_dot_uHat * uHat[0] + vprimeVec[0];
	delrprimeVec[1] = delrVec_dot_uHat * uHat[1] + vprimeVec[1];
	delrprimeVec[2] = delrVec_dot_uHat * uHat[2] + vprimeVec[2];

	delx = delrVec[0] - delrprimeVec[0]; // Random number from 0 -> 1 multiplied by maximum movement for each Cartesian coordinate
	dely = delrVec[1] - delrprimeVec[1];
	delz = delrVec[2] - delrprimeVec[2];

	xold = x[k]; // Save current position as old position
	yold = y[k];
	zold = z[k];

	x[k] += delx; // Change position of current particle
	y[k] += dely;
	z[k] += delz;

	/*      if (position_check(x, y, z))
		{ // If positon check returns 1, -> accept
		acptRun += 1;
	// printf("%d\n", k);
	}
	else
	{ // If position check returns 0 -> disregard
	x[k] = xold;
	y[k] = yold;
	z[k] = zold;
	dsrgRun += 1;
	// printf("%d\n", k);
	}*/
}

int main() {
	int n = 1000; // Number of elements to generate
	double *d_output = nullptr; // Device array to store random numbers
	double *h_output = nullptr; // Host array to store random numbers
	double *xhost=nullptr, *yhost=nullptr, *zhost=nullptr;
	double *x=nullptr, *y=nullptr, *z=nullptr;
	long threadsPerBlock, blocksPerGrid, ncyc = 10;

	int imov = 1;

	FILE *fp;

	xhost = (double *)malloc(n * sizeof(double));
	if (xhost == nullptr) {
		printf("Memory allocation error on the host xhost.\n");
		return 1;
	}

	yhost = (double *)malloc(n * sizeof(double));
	if (yhost == nullptr) {
		printf("Memory allocation error on the host yhost.\n");
		return 1;
	}

	zhost = (double *)malloc(n * sizeof(double));
	if (zhost == nullptr) {
		printf("Memory allocation error on the host zhost.\n");
		return 1;
	}

	// Allocate memory on the host
	h_output = (double *)malloc(n * sizeof(double));
	if (h_output == nullptr) {
		printf("Memory allocation error on the host h_output.\n");
		return 1;
	}

	// Allocate memory on the device
	cudaError_t cudaStatus = cudaMalloc((void **)&d_output, n * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		printf("cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
		free(h_output);
		return 1;
	}

	cudaStatus = cudaMalloc((void **)&x, n*sizeof(double));
	if (cudaStatus != cudaSuccess) {
		printf("cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
		free(h_output);
		return 1;
	}

	cudaStatus = cudaMalloc((void **)&y, n*sizeof(double));
	if (cudaStatus != cudaSuccess) {
		printf("cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
		free(h_output);
		return 1;
	}

	cudaStatus = cudaMalloc((void **)&z, n*sizeof(double));
	if (cudaStatus != cudaSuccess) {
		printf("cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
		free(h_output);
		return 1;
	}

	if (imov == 1) 
	{
		fp = fopen("chain.xyz", "w");
	}

	for(long ii=0; ii<ncyc; ii++)
	{
		// Launch the kernel
		threadsPerBlock = 1;
		blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;
		randomKernel<<<blocksPerGrid, threadsPerBlock>>>(d_output, n);

		// Check for kernel launch errors
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			printf("Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			cudaFree(d_output);
			free(h_output);
			return 1;
		}

		threadsPerBlock = 1;
		blocksPerGrid = 1;
		init_pos_dev<<<blocksPerGrid, threadsPerBlock>>>(x, y, z, n);

		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) { 
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
		cudaStatus = cudaMemcpy(h_output, d_output, n * sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			printf("cudaMemcpy failed: %s\n", cudaGetErrorString(cudaStatus));
			cudaFree(d_output);
			free(h_output);
			return 1;
		}

		cudaStatus = cudaMemcpy(xhost, x, n*sizeof(double), cudaMemcpyDeviceToHost);
		cudaStatus = cudaMemcpy(yhost, y, n*sizeof(double), cudaMemcpyDeviceToHost);
		cudaStatus = cudaMemcpy(zhost, z, n*sizeof(double), cudaMemcpyDeviceToHost);

		if (cudaStatus != cudaSuccess) { 
			printf("cudaMemcpy failed: %s\n", cudaGetErrorString(cudaStatus));
			cudaFree(x);
			cudaFree(y);
			cudaFree(z);
			free(xhost);
			free(yhost);
			free(zhost);
			return 1;
		}

		// Print the generated random numbers on the host
		if (imov == 1) 
		{
			for (int i = 0; i < n; ++i) 
			{
				fprintf(fp, "%lf\nPolymer: %ld \n", n, ncyc);
				fprintf(fp, "%lf  %lf  %lf\n", xhost[i], yhost[i], zhost[i]);
			}
		}
	}
	
	if (imov == 1) 
	{
		fclose(fp);
	}
	
	// Free device and host memory
	cudaFree(d_output);
	free(h_output);

	return 0;
}

