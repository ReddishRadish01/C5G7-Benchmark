#pragma once

#include "cudaHeader.cuh"
#include "XSParser.cuh"
#include "Neutron.cuh"
#include "Core.cuh"

#include <vector>
#include <algorithm>

struct NeutronBankView {
    Neutron* neutrons;
    Neutron* addedNeutrons;
    int neutronSize;
    int allocatableNeutronNum;
    int addedNeutronSize;
    int addedNeutronIndex;
    unsigned long long seedNo;
};

class GPU_Manager {
	//GPU_Manager() = default;
public:

	// double pointer needed - or put MatXS*& d_ptr. I kept the double pointer so that it stays safe with CUDA style codes. 
	H static void C5G7DeviceAllocater(XSLibrary** d_ptr, XSLibrary& h_instance) {
		cudaMalloc(d_ptr, sizeof(h_instance));
		cudaMemcpy(*d_ptr, &h_instance, sizeof(h_instance), cudaMemcpyHostToDevice);
	}
	// More: if the d_ptr is passed by MatXS* d_ptr, it is passed as value - it doesn't actually change the value(the address) of d_ptr.
	// thus, you need to pass it as the pointer to pointer (MatXS**), or the reference to pointer (MatXS*&)

	

    H static inline void compact_bank_host(NeutronBank& bank) {
        std::vector<Neutron> tmp;
        tmp.reserve(static_cast<size_t>(bank.neutronSize + bank.addedNeutronSize));

        for (int j = 0; j < bank.allocatableNeutronNum; ++j) {
            if (!bank.neutrons[j].isNullified()) {
                bank.neutrons[j].passFlag = false;
                tmp.push_back(bank.neutrons[j]);
            }
            if (!bank.addedNeutrons[j].isNullified()) {
                bank.addedNeutrons[j].passFlag = false;
                tmp.push_back(bank.addedNeutrons[j]);
            }
            bank.neutrons[j].Nullify();
            bank.addedNeutrons[j].Nullify();
        }

        const int cap = bank.allocatableNeutronNum;
        const int total = static_cast<int>(tmp.size());
        const int n0 = std::min(total, cap);
        const int n1 = total - n0;

        for (int j = 0; j < n0; ++j) bank.neutrons[j] = tmp[j];
        for (int j = 0; j < n1; ++j) bank.addedNeutrons[j] = tmp[n0 + j];

        bank.neutronSize = n0;
        bank.addedNeutronSize = n1;
        bank.addedNeutronIndex = n1;
    }

    static inline void compact_bank_device(NeutronBank* d_Bank) {
        NeutronBankView dv;
        cudaMemcpy(&dv, d_Bank, sizeof(NeutronBankView), cudaMemcpyDeviceToHost);

        const int cap = dv.allocatableNeutronNum;

        std::vector<Neutron> hN(static_cast<size_t>(cap));
        std::vector<Neutron> hA(static_cast<size_t>(cap));

        cudaMemcpy(hN.data(), dv.neutrons, sizeof(Neutron) * cap, cudaMemcpyDeviceToHost);
        cudaMemcpy(hA.data(), dv.addedNeutrons, sizeof(Neutron) * cap, cudaMemcpyDeviceToHost);

        std::vector<Neutron> tmp;
        tmp.reserve(static_cast<size_t>(dv.neutronSize + dv.addedNeutronSize));

        for (int j = 0; j < cap; ++j) {
            if (!hN[j].isNullified()) { hN[j].passFlag = false; tmp.push_back(hN[j]); }
            if (!hA[j].isNullified()) { hA[j].passFlag = false; tmp.push_back(hA[j]); }
            hN[j].Nullify();
            hA[j].Nullify();
        }

        const int total = static_cast<int>(tmp.size());
        const int n0 = std::min(total, cap);
        const int n1 = total - n0;

        for (int j = 0; j < n0; ++j) hN[j] = tmp[j];
        for (int j = 0; j < n1; ++j) hA[j] = tmp[n0 + j];

        cudaMemcpy(dv.neutrons, hN.data(), sizeof(Neutron) * cap, cudaMemcpyHostToDevice);
        cudaMemcpy(dv.addedNeutrons, hA.data(), sizeof(Neutron) * cap, cudaMemcpyHostToDevice);

        dv.neutronSize = n0;
        dv.addedNeutronSize = n1;
        dv.addedNeutronIndex = n1;

        cudaMemcpy(d_Bank, &dv, sizeof(NeutronBankView), cudaMemcpyHostToDevice);
    }

};
