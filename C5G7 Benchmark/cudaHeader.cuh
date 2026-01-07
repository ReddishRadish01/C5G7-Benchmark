#pragma once

#include <cuda.h>
#include <device_launch_parameters.h>
#include <cuda_runtime.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
//#include <curand.h>
#include "Constants.cuh"
#include "RNG.cuh"
//#include "XSParser.cuh"

#ifdef __CUDACC__
	#define HD __host__ __device__
	#define H  __host__
	#define D  __device__
	#define G  __global__
#else
	#define HD
	#define H
	#define D
	#define G
#endif

