
#include "cudaHeader.cuh"
#include "XSParser.cuh"
#include "Neutron.cuh"
#include "Core.cuh"
#include "GpuManager.cuh"
#include "CoreManager.cuh"
#include "XSManager.cuh"
#include "Interaction.cuh"
#include "NeutronBankManager.cuh"
#include "Tally.cuh"
#include "Debug.cuh"
#include "Cycle.cuh"

#include <iomanip>
#include <ctime>
#include <sstream>

//#define CPURUN
#define GPURUN

//#define XSRESULTDEBUG

G void GPUTest(int num, XSLibrary* d_MatXS, C5G7Geometry* Core, NeutronBank* bank, unsigned long long* seedArr) {
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < num) {
		GnuAMCM RNG(seedArr[idx]);
		//printf("Cross Section: %.9f\n", d_MatXS->UO2.transXS[0]);
		Neutron localNeutron = bank->neutrons[idx];

		//MatType meatType = Core->returnAssemblyByPos(bank->neutrons[idx]).returnPincellByPos(bank->neutrons[idx]).meatType;
		Pincell currentPincell = CoreManager::returnPincellByPos(Core, localNeutron);
		Assembly currentAssembly = Core->returnAssemblyByNeutron(localNeutron);
		vec3 floorNeutronPos = currentAssembly.returnFlooredNeutronPosInPincell(localNeutron);
		//Pincell thisPincell = currentAssembly.returnPincellByPos(localNeutron);
		double DTC = currentAssembly.DTC(localNeutron, d_MatXS, RNG);
		double DTS = currentAssembly.DTS(localNeutron);
		MatType meatType = currentPincell.meatType;
		MatType modType = currentPincell.modType;

		MatType currentType = currentPincell.meatOrMod(floorNeutronPos);
		double outEnergy = 0.0;
		InteractionType interaction = XSManager::returnInteracitonType(d_MatXS, currentType, RNG, localNeutron.energy, outEnergy);
		//printf("idx: %d, Material type : %s\t at Neutron Pos (%2.3f, %2.3f, %2.3f), DTC:%2.5f \n", idx, to_string(meatType), localNeutron.pos.x, localNeutron.pos.y, localNeutron.pos.z, DTC);
		//printf("idx: %d, pin type: %s, current material: %s, DTC: %f, TransXS: %f\n", idx, to_string(meatType), to_string(currentType), DTC, d_MatXS->returnMatByType(currentType).transXS[static_cast<int>(localNeutron.energy) - 1]);
		/*
		printf("idx: %d, pin type: %s, current material: %s, N Pos (%f, %f, %f), floor Neutron Pos (%2.3f, %2.3f, %2.3f), dirVec: (%f, %f, %f) DTC: %f, DTS: %f\n", idx, to_string(meatType), to_string(currentType),
			localNeutron.pos.x, localNeutron.pos.y, localNeutron.pos.z, floorNeutronPos.x, floorNeutronPos.y, floorNeutronPos.z,
			localNeutron.dirVec.x, localNeutron.dirVec.y, localNeutron.dirVec.z, DTC, DTS
		);
		*/
		printf("idx: %d, pin type: %s, current material: %s, Interaction type: %s\n",
			idx, to_string(meatType), to_string(currentType), to_string(interaction)
		);
	}
}

H void CPUTest(int num, XSLibrary* d_MatXS, C5G7Geometry* Core, NeutronBank* bank, unsigned long long* seedArr) {
	for(int idx = 0; idx < num; idx++) {
		GnuAMCM RNG(seedArr[idx]);
		Neutron localNeutron = bank->neutrons[idx];

		//MatType meatType = Core->returnAssemblyByPos(bank->neutrons[idx]).returnPincellByPos(bank->neutrons[idx]).meatType;
		Pincell currentPincell = CoreManager::returnPincellByPos(Core, localNeutron);
		Assembly currentAssembly = Core->returnAssemblyByNeutron(localNeutron);
		vec3 floorNeutronPos = currentAssembly.returnFlooredNeutronPosInPincell(localNeutron);
		Pincell thisPincell = currentAssembly.returnPincellByPos(localNeutron);
		double DTC = currentAssembly.DTC(localNeutron, d_MatXS, RNG);
		double DTS = currentAssembly.DTS(localNeutron);
		MatType meatType = currentPincell.meatType;
		MatType modType = currentPincell.modType;

		MatType currentType = currentPincell.meatOrMod(floorNeutronPos);
		double outEnergy = 0.0;
		InteractionType interaction = XSManager::returnInteracitonType(d_MatXS, currentType, RNG, localNeutron.energy, outEnergy);
		//printf("idx: %d, Material type : %s\t at Neutron Pos (%2.3f, %2.3f, %2.3f), DTC:%2.5f \n", idx, to_string(meatType), localNeutron.pos.x, localNeutron.pos.y, localNeutron.pos.z, DTC);
		//printf("idx: %d, pin type: %s, current material: %s, DTC: %f, TransXS: %f\n", idx, to_string(meatType), to_string(currentType), DTC, d_MatXS->returnMatByType(currentType).transXS[static_cast<int>(localNeutron.energy) - 1]);

		/*
		printf("idx: %d, pin type: %s, current material: %s, N Pos (%f, %f, %f), floor Neutron Pos (%2.3f, %2.3f, %2.3f), dirVec: (%f, %f, %f) DTC: %f, DTS: %f\n", idx, to_string(meatType), to_string(currentType),
			localNeutron.pos.x, localNeutron.pos.y, localNeutron.pos.z, floorNeutronPos.x, floorNeutronPos.y, floorNeutronPos.z,
			localNeutron.dirVec.x, localNeutron.dirVec.y, localNeutron.dirVec.z, DTC, DTS
		);
		*/

		printf("idx: %d, pin type: %s, current material: %s, Interaction type: %s, curEnergy: %f, outEnergy: %f\n",
			idx, to_string(meatType), to_string(currentType), to_string(interaction), localNeutron.energy, outEnergy
		);
	}
}

int main() {
	int num = 500000;
	int numCycle = 700;
	int inactiveCycle = 250;

	int threadPerBlock = 32;
	int blockPerDim = (num + threadPerBlock - 1) / threadPerBlock;

	double h_multK = 0.0;
	std::cout << "input the initial K muliplication factor:\n";
	//std::cin >> h_multK;
	h_multK = 1.1865;

	double* d_multK = nullptr;
	cudaMalloc(&d_multK, sizeof(double));
	cudaMemcpy(d_multK, &h_multK, sizeof(double), cudaMemcpyHostToDevice);

	unsigned long long seedNo = 92235922381;

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
#ifdef XSRESULTDEBUG
	std::cout << h_UO2XS.transXS[0] << " " << h_Mod.transXS[2] << "\n";
#endif

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

	XSLibrary h_XSLib{};
	MatXSFactory::initialize(h_XSLib, XS);
	XSLibrary* d_XSLib = nullptr;
	GPU_Manager::C5G7DeviceAllocater(&d_XSLib, h_XSLib);


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


	// fetch K to text file

	std::time_t t = std::time(nullptr);
	std::tm tm = *std::localtime(&t);

	std::ostringstream oss;
	oss << "k_history_"
		<< (tm.tm_year + 1900)
		<< std::setw(2) << std::setfill('0') << (tm.tm_mon + 1)
		<< std::setw(2) << std::setfill('0') << tm.tm_mday
		<< "_"
		<< std::setw(2) << std::setfill('0') << tm.tm_hour
		<< std::setw(2) << std::setfill('0') << tm.tm_min
		<< std::setw(2) << std::setfill('0') << tm.tm_sec
		<< ".txt";

	std::ofstream klog(oss.str(), std::ios::out);
	klog << std::fixed << std::setprecision(6);


	for (int i = 0; i < h_Bank.neutronSize; i++) {
		vec3 randPos = { h_RNG.uniform(0, h_Core.x), h_RNG.uniform(0, h_Core.y), h_RNG.uniform(0, h_Core.z) };
		h_Bank.neutrons[i] = Neutron(randPos, vec3::randomUnit(h_RNG), static_cast<double>(h_RNG.int_dist(1, 7)), 1.0);
		//h_Bank.neutrons[i] = Neutron({ 0.1, 0.1, 0.1 }, { -1.0, 0.0, 0.0 }, static_cast<double>(h_RNG.int_dist(1, 7)), 1.0);
		//h_Bank.neutrons[i] = Neutron({ 0.63, 0.63, 0.63 }, { 1.0, 0.0, 0.0 }, static_cast<double>(h_RNG.int_dist(1, 7)), 1.0);
		//h_Bank.neutrons[i] = Neutron(randPos, vec3::randomUnit(h_RNG), 1, 1.0);
		h_Bank.addedNeutrons[i] = Neutron();
		//h_Bank.addedNeutrons[i].status = false;
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

	
	//for (int i = 0; i < 10; i++) {	Debug::fuelLayoutDebug(h_Core.assembly[i]);	}
	
	//CPUTest(num, &h_XSLib, &h_Core, &h_Bank, h_SeedArr);
	//GPUTest << <blockPerDim, threadPerBlock >> > (num, d_MatXS, d_Core, d_Bank, d_SeedArr);

	double tempK = h_multK;
	double previousNumNeutron = h_Bank.getTotalNeutronNum();
	double kAvg = 0.0;

	for (int i = 0; i < numCycle; i++) {
		double absorption = 0.0;
		double fission = 0.0;
		double leak = 0.0;
#ifdef CPURUN
		cycle_Neutron << <blockPerDim, threadPerBlock >> > (d_Bank, d_Core, d_XSLib, d_SeedArr, d_multK, false);
		addedNeutronPassResetter << <blockPerDim, threadPerBlock >> > (d_Bank);
		cycle_addedNeutron << <blockPerDim, threadPerBlock >> > (d_Bank, d_Core, d_XSLib, d_SeedArr, d_multK, true);

		cudaMemcpy(&h_Bank, d_Bank, sizeof(NeutronBank), cudaMemcpyDeviceToHost);
		double currentNumNeutron = h_Bank.getTotalNeutronNum();
		std::cout << "Cycle " << i + 1 << ", currentNum: " << currentNumNeutron;
		h_multK = currentNumNeutron / previousNumNeutron;
		std::cout << "\tk:" << h_multK << "\n";
		cudaMemcpy(d_Bank, &h_Bank, sizeof(NeutronBank), cudaMemcpyHostToDevice);

#ifdef INTERACTIONDEBUG
		std::cout << "addedNeutron: \n";
#endif
		cycle_addedNeutron_CPU(&h_Bank, &h_Core, &h_XSLib, h_SeedArr, &h_multK, true, absorption, fission, leak);
		addedNeutronPassResetter_CPU(&h_Bank);
#ifdef INTERACTIONDEBUG
		std::cout << "\nNeutron:\n";
#endif
		cycle_Neutron_CPU(&h_Bank, &h_Core, &h_XSLib, h_SeedArr, &h_multK, false, absorption, fission, leak);
#endif
		
#ifdef GPURUN
		cycle_addedNeutron << <blockPerDim, threadPerBlock >> > (d_Bank, d_Core, d_XSLib, d_SeedArr, d_multK, true);
		addedNeutronPassResetter<<<blockPerDim, threadPerBlock>>>(d_Bank);
		cycle_Neutron << <blockPerDim, threadPerBlock >> > (d_Bank, d_Core, d_XSLib, d_SeedArr, d_multK, false);

		cudaMemcpy(&h_Bank, d_Bank, sizeof(NeutronBank), cudaMemcpyDeviceToHost);
#endif
		cudaMemcpy(&h_multK, d_multK, sizeof(double), cudaMemcpyDeviceToHost);
		double currentNumNeutron = h_Bank.getTotalNeutronNum();
		double oldK = h_multK;
		std::cout << "Cycle " << i + 1 << ", currentNum: " << currentNumNeutron;
		h_multK = h_multK * currentNumNeutron / previousNumNeutron;
		previousNumNeutron = currentNumNeutron;
		std::cout << "\tk: " << h_multK << "\n";
		//std::cout << "\t n count : " << h_Bank.neutronSize << " addn count : " << h_Bank.addedNeutronSize << " addN addIndex : " << h_Bank.addedNeutronIndex;

		cudaMemcpy(d_multK, &h_multK, sizeof(double), cudaMemcpyHostToDevice);
#ifdef CPURUN
		std::cout << "  capture: " << absorption << " , fission neutron num: " << fission << ", leak: " << leak << "\n";
#endif
		klog << (i + 1) << " " << h_multK << "\n";
		if (i > inactiveCycle) {
			kAvg += h_multK;
		}

		
		if (h_Bank.addedNeutronIndex > h_Bank.allocatableNeutronNum * 0.8) {
			std::cout << "Merging:\n";
#ifdef CPURUN
			std::vector<Neutron> NeutronContainer;
			for (int j = 0; j < h_Bank.allocatableNeutronNum; j++) {
				if (!h_Bank.neutrons[j].isNullified()) {
					NeutronContainer.reserve(1);
					h_Bank.neutrons[j].passFlag = false;
					NeutronContainer.emplace_back(h_Bank.neutrons[j]);
					h_Bank.neutrons[j].Nullify();	// flush all the old neutrons
				}
				if (!h_Bank.addedNeutrons[j].isNullified()) {
					NeutronContainer.reserve(1);
					h_Bank.addedNeutrons[j].passFlag = false;
					NeutronContainer.emplace_back(h_Bank.addedNeutrons[j]);
					h_Bank.addedNeutrons[j].Nullify();
				}
			}

			if (h_Bank.allocatableNeutronNum < NeutronContainer.size()) {
				h_Bank.neutronSize = h_Bank.allocatableNeutronNum;
				h_Bank.addedNeutronSize = NeutronContainer.size() - h_Bank.allocatableNeutronNum;
				h_Bank.addedNeutronIndex = h_Bank.addedNeutronSize;
				std::cout << "Neutron vector size: " << NeutronContainer.size() << "\n";

				for (int j = 0; j < h_Bank.allocatableNeutronNum; j++) {
					h_Bank.neutrons[j] = NeutronContainer[j];
				}
				for (int j = 0; j < h_Bank.addedNeutronSize; j++) {
					h_Bank.addedNeutrons[j] = NeutronContainer[h_Bank.neutronSize + j];
				}

				std::cout << "After sorting: Neutron size: " << h_Bank.getTotalNeutronNum() << "\n";
			}
			else {
				h_Bank.neutronSize = NeutronContainer.size();
				h_Bank.addedNeutronIndex = 0;
				h_Bank.addedNeutronSize = 0;
				std::cout << "Neutron vector size: " << NeutronContainer.size() << "\n";

				for (int j = 0; j < h_Bank.neutronSize; j++) {
					h_Bank.neutrons[j] = NeutronContainer[j];
				}
				std::cout << "After sorting: Neutron size: " << h_Bank.getTotalNeutronNum() << "\n"
			}
		}
#endif

#ifdef GPURUN
			GPU_Manager::compact_bank_device(d_Bank);
		}
#endif
	}

	kAvg /= (numCycle - inactiveCycle);
	
	std::cout << "\n initial Neutron number: " << num << ", for cycle of " << numCycle - inactiveCycle << ", average k = " << kAvg << "\n";

	klog.close();
}