#include<stdio.h>
#include <curand_kernel.h> // Include the CUDA Random Number Generation library
#include <curand.h>

// Device function to generate a random number using curand
__device__ double randomNumber(curandState *state) {
    // Generate a random double in the range [0, 1)
    return curand_uniform_double(state);
}

// CUDA kernel that uses the random number generator
__global__ void randomKernel(double *output, int n) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < n) {
        // Initialize the random number generator state for each thread
        curandState state;
        curand_init(clock64(), tid, 0, &state);

        // Generate a random number for each thread and store it in the output array
        output[tid] = randomNumber(&state);
    }
}

int main() {
    int n = 1; // Number of elements to generate
    double *d_output = nullptr; // Device array to store random numbers
    double *h_output = nullptr; // Host array to store random numbers

    // Allocate memory on the host
    h_output = (double *)malloc(n * sizeof(double));
    if (h_output == nullptr) {
        printf("Memory allocation error on the host.\n");
        return 1;
    }

    // Allocate memory on the device
    cudaError_t cudaStatus = cudaMalloc((void **)&d_output, n * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        printf("cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
        free(h_output);
        return 1;
    }

    // Launch the kernel
    int threadsPerBlock = 1;
    int blocksPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;
    randomKernel<<<blocksPerGrid, threadsPerBlock>>>(d_output, n);

    // Check for kernel launch errors
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        printf("Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        cudaFree(d_output);
        free(h_output);
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

    // Print the generated random numbers on the host
    for (int i = 0; i < n; ++i) {
        printf("Random Number %d: %f\n", i, h_output[i]);
    }

    // Free device and host memory
    cudaFree(d_output);
    free(h_output);

    return 0;
}
