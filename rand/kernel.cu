#define COUNT 100

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <curand.h>
#include <curand_kernel.h>
#include <math.h>
#include <assert.h>


unsigned int Rand(unsigned int randx)
{
    randx = randx * 1103515245 + 12345;
    return randx & 2147483647;
}


 int* random(int *c)
{
    int i = 0;
    c[i] = 100;
    for (i = 0;i < COUNT;i++) {
        c[i+1] = Rand(c[i]); 
    }
    return c;
}

int main()
{
    /*int* c;
    cudaMalloc(&c, COUNT * sizeof(int));*/
    int *c;
    int i;
    c = (int*)malloc(COUNT * sizeof(int));
    c = random(c);
    for (i = 0;i < COUNT;i++) {
        std::cout << "rand : " << float(c[i])/ 2147483647 << "\n";
    }
}

// Helper function for using CUDA to add vectors in parallel.

