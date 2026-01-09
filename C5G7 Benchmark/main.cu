
#include "cudaHeader.cuh"
#include "XSParser.cuh"
#include "Neutron.cuh"
#include "Core.cuh"
#include "GpuManager.cuh"
#include "CoreManager.cuh"
#include "XSManager.cuh"
#include "NeutronBankManager.cuh"
#include "Tally.cuh"
#include "Debug.cuh"



G void GPUTest(int num, XSLibrary* d_MatXS, C5G7Geometry* Core, NeutronBank* bank) {
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < num) {
		//printf("Cross Section: %f\n", d_MatXS->UO2.totalXS[1]);
		Neutron localNeutron = bank->neutrons[idx];
		//MatType meatType = Core->returnAssemblyByPos(bank->neutrons[idx]).returnPincellByPos(bank->neutrons[idx]).meatType;
		MatType meatType = CoreManager::returnPincellByPos(Core, localNeutron).meatType;
		printf("idx: %d, Material type : %s\t at Neutron Pos (%2.3f, %2.3f, %2.3f)\n", idx, to_string(meatType), localNeutron.pos.x, localNeutron.pos.y, localNeutron.pos.z);
	}
}

H void CPUTest(int num, XSLibrary* d_MatXS, C5G7Geometry* Core, NeutronBank* bank) {
	for(int idx = 0; idx < num; idx++) {
		//printf("Cross Section: %f\n", d_MatXS->UO2.totalXS[1]);
		Neutron localNeutron = bank->neutrons[idx];
		//MatType meatType = Core->returnAssemblyByPos(bank->neutrons[idx]).returnPincellByPos(bank->neutrons[idx]).meatType;
		MatType meatType = CoreManager::returnPincellByPos(Core, localNeutron).meatType;
		printf("idx: %d, Material type : %s\t at Neutron Pos (%2.3f, %2.3f, %2.3f)\n", idx, to_string(meatType), localNeutron.pos.x, localNeutron.pos.y, localNeutron.pos.z);
	}
}

int main() {
	int num = 100;
	int threadPerBlock = 32;
	int blockPerDim = (num + threadPerBlock - 1) / threadPerBlock;

	unsigned long long seedNo = 9223592237;

	GnuAMCM h_RNG(seedNo);
	unsigned long long* h_SeedArr = new unsigned long long[num];
	for (int i = 0; i < num; i++) {
		h_SeedArr[i] = (h_RNG.gen() + i) & (0xFFFFFFFFFFFFULL);
	}
	unsigned long long* d_SeedArr = nullptr;
	cudaMalloc(&d_SeedArr, num * sizeof(unsigned long long));
	cudaMemcpy(d_SeedArr, h_SeedArr, num * sizeof(unsigned long long), cudaMemcpyHostToDevice);



	std::vector<MatXS> XS;
	MatXS h_UO2XS("C5txt/UO2.txt", MatType::UO2);
	MatXS h_MOX4_3("C5txt/Mox4_3.txt", MatType::MOX4_3);
	MatXS h_MOX7_0("C5txt/Mox7_0.txt", MatType::MOX7_0);
	MatXS h_MOX8_7("C5txt/Mox8_7.txt", MatType::MOX8_7);
	MatXS h_FC("C5txt/FC.txt", MatType::FC);
	MatXS h_GT("C5txt/GT.txt", MatType::GT);
	MatXS h_Mod("C5txt/Mod.txt", MatType::MOD);

	XS.reserve(7);
	//XS.push_back(h_UO2XS); // uncomfortable
	//XS.push_back(std::move(h_UO2XS)); //h_UO2XS is deprecated
	XS.emplace_back(h_UO2XS);
	XS.emplace_back(h_MOX4_3);
	XS.emplace_back(h_MOX7_0);
	XS.emplace_back(h_MOX8_7);
	XS.emplace_back(h_FC);
	XS.emplace_back(h_GT);
	XS.emplace_back(h_Mod);

	XSLibrary h_MatXS{};
	MatXSFactory::initialize(h_MatXS, XS);
	XSLibrary* d_MatXS = nullptr;
	GPU_Manager::C5G7DeviceAllocater(&d_MatXS, h_MatXS);


	C5G7Geometry h_Core{};
	C5G7GeometryFactory::Initialize(h_Core, "Geometry/C5G7CoreGeometry.txt", "Geometry/UO2Geometry.txt", "Geometry/MOXGeometry.txt");

	C5G7Geometry* d_Core = nullptr;
	cudaMalloc(&d_Core, sizeof(C5G7Geometry));
	
	// initialize vector with nullptr, number is given.
	Assembly* d_bufferAssembly = nullptr;
	cudaMalloc(&d_bufferAssembly, sizeof(Assembly) * h_Core.assemblyNo);

	std::vector<Pincell*> d_bufferPincellVec(h_Core.assemblyNo, nullptr);
	std::vector<Assembly> tmp_Assembly(h_Core.assemblyNo);

	for (int i = 0; i < d_bufferPincellVec.size(); i++) {
		int n = h_Core.assembly[i].totalPincellNo();

		cudaMalloc(&d_bufferPincellVec[i], sizeof(Pincell) * h_Core.assembly[i].totalPincellNo());
		cudaMemcpy(d_bufferPincellVec[i], h_Core.assembly[i].pinCells, sizeof(Pincell) * n, cudaMemcpyHostToDevice);
		
		tmp_Assembly[i] = h_Core.assembly[i];
		tmp_Assembly[i].pinCells = d_bufferPincellVec[i];

		
	}
	
	cudaMemcpy(d_bufferAssembly, tmp_Assembly.data(), sizeof(Assembly) * h_Core.assemblyNo, cudaMemcpyHostToDevice);

	C5G7Geometry tmp_Core = h_Core;
	tmp_Core.assembly = d_bufferAssembly;

	cudaMemcpy(d_Core, &tmp_Core, sizeof(C5G7Geometry), cudaMemcpyHostToDevice);

	Neutron h_Neutron{ {4.0, 4.0, 200.0}, {0.0, 0.0, 0.0}, 1.0, 1.0 };
	Neutron* d_Neutron = nullptr;
	cudaMalloc(&d_Neutron, sizeof(Neutron));
	cudaMemcpy(d_Neutron, &h_Neutron, sizeof(Neutron), cudaMemcpyHostToDevice);

	NeutronBank h_Bank(num, seedNo);

	for (int i = 0; i < h_Bank.neutronSize; i++) {
		vec3 randPos = { h_RNG.uniform(0, h_Core.x), h_RNG.uniform(0, h_Core.y), h_RNG.uniform(0, h_Core.z) };
		h_Bank.neutrons[i] = Neutron(randPos, vec3::randomUnit(h_RNG), static_cast<double>(h_RNG.int_dist(1, 7)), 1.0);
		h_Bank.addedNeutrons[i] = Neutron({0,0,0}, {0,0,0}, 0.0, 0.0);
	}

	NeutronBank* d_Bank = nullptr;
	cudaMalloc(&d_Bank, sizeof(NeutronBank));
	Neutron* d_bufferNeutrons = nullptr; Neutron* d_bufferAddedNeutrons = nullptr;
	cudaMalloc(&d_bufferNeutrons, h_Bank.allocatableNeutronNum * sizeof(Neutron));
	cudaMemcpy(d_bufferNeutrons, h_Bank.neutrons, h_Bank.allocatableNeutronNum * sizeof(Neutron), cudaMemcpyHostToDevice);
	cudaMalloc(&d_bufferAddedNeutrons, h_Bank.allocatableNeutronNum * sizeof(Neutron));
	cudaMemcpy(d_bufferAddedNeutrons, h_Bank.addedNeutrons, h_Bank.allocatableNeutronNum * sizeof(Neutron), cudaMemcpyHostToDevice);
	NeutronBank tmp_Bank = h_Bank;
	tmp_Bank.neutrons = d_bufferNeutrons;
	tmp_Bank.addedNeutrons = d_bufferAddedNeutrons;
	cudaMemcpy(d_Bank, &tmp_Bank, sizeof(NeutronBank), cudaMemcpyHostToDevice);
	

	Neutron h_testNeutron{ {4.0, 4.0, 200.0}, {0.0, 0.0, 0.0}, 1.0, 1.0 };


	/*
	for (int i = 0; i < 10; i++) {
		Debug::fuelLayoutDebug(h_Core.assembly[i]);
	}
	*/
	

	
	//CPUTest(num, &h_MatXS, &h_Core, &h_Bank);
	GPUTest << <blockPerDim, threadPerBlock >> > (num, d_MatXS, d_Core, d_Bank);

}