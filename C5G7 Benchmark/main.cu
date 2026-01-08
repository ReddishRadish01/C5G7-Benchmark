
#include "cudaHeader.cuh"
#include "XSParser.cuh"
#include "Neutron.cuh"
#include "Core.cuh"
#include "GpuManager.cuh"
#include "CoreManager.cuh"
#include "XSManager.cuh"
#include "NeutronBankManager.cuh"
#include "Debug.cuh"


G void helloGPU(XSLibrary* d_MatXS) {
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	printf("Cross Section: %f\n", d_MatXS->UO2.totalXS[1]);
	
}

int main() {
	int num = 100;
	int threadPerBlock = 32;
	int blockPerDim = (num + threadPerBlock - 1) / threadPerBlock;

	std::vector<MatXS> XS;
	std::vector<Neutron> Bank;

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

	/*
	Assembly UO2asm{};
	UO2asm.Initialize("Geometry/UO2Geometry.txt");
	Debug::fuelLayoutDebug(UO2asm);
	std::cout << "\n\n";

	Assembly MOXasm{};
	MOXasm.Initialize("Geometry/MOXGeometry.txt");
	Debug::fuelLayoutDebug(MOXasm);
	*/

	C5G7Geometry Core{};
	C5G7GeometryFactory::Initialize(Core, "Geometry/C5G7CoreGeometry.txt", "Geometry/UO2Geometry.txt", "Geometry/MOXGeometry.txt");
	
	for (int i = 0; i < 10; i++) {
		Debug::fuelLayoutDebug(Core.assembly[i]);
	}
	
	//Debug::fuelLayoutDebug(Core.assembly[0]);

	GPU_Manager::C5G7DeviceAllocater(&d_MatXS, h_MatXS);

	//helloGPU << <blockPerDim, threadPerBlock >> > (d_MatXS);


	
}