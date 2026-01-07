#pragma once

#include "cudaHeader.cuh"
#include "XSParser.cuh"
#include "Neutron.cuh"
#include "Core.cuh"


class GPU_Manager {
	//GPU_Manager() = default;
public:

	// double pointer needed - or put MatXS*& d_ptr. I kept the double pointer so that it stays safe with CUDA style codes. 
	H static void C5G7DeviceAllocater(MatXS** d_ptr, MatXS& h_instance) {
		cudaMalloc(d_ptr, sizeof(h_instance));
		cudaMemcpy(*d_ptr, &h_instance, sizeof(h_instance), cudaMemcpyHostToDevice);
	}
	// More: if the d_ptr is passed by MatXS* d_ptr, it is passed as value - it doesn't actually change the value(the address) of d_ptr.
	// thus, you need to pass it as the pointer to pointer (MatXS**), or the reference to pointer (MatXS*&)


};
