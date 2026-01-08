#pragma once
#include "cudaHeader.cuh"
#include "XSParser.cuh"
#include "Neutron.cuh"




class Pincell {
public:
	double sideLength = 1.26;	// 1.26 cm. alias with height
	double cellHeight = 1.26;
	double radius = 0.54;	// 0.54 cm radius
	MatType meatType = MatType::Unknown;
	MatType modType = MatType::MOD;
	vec2 centerOffset = { 0.0, 0.0 };

	H Pincell() = default;

	H Pincell(double sideLength, double radius, double height = 0.0, MatType meat = MatType::Unknown, MatType mod = MatType::MOD, vec2 centerOffset = { 0.0, 0.0 })
		: sideLength(sideLength), radius(radius), cellHeight(height), meatType(meat), modType(mod), centerOffset(centerOffset)
	{
		// exceptions: if radius is bigger than sideLength / 2 ? remove it
		if (2 * radius > sideLength) { radius = sideLength / 2.0; }
	}

	HD MatType meatOrMod(vec3 flooredPos) {
		// this is for moderator block
		if (radius == 0.0) { return modType;  }

		// rest is for regular pincell
		vec2 center = { this->sideLength / 2.0 + centerOffset.x, this->sideLength / 2.0 + centerOffset.y };
		if ((flooredPos.x - center.x) * (flooredPos.x - center.x) + (flooredPos.y - center.y) * (flooredPos.y - center.y) > (radius * radius))
			return modType;
		else 
			return meatType;
	}
};

class Assembly {
public:
	vec3 startPos = { 0, 0, 0 };

	vec3 length = { 0, 0, 0 };

	int xNum = 0;
	int yNum = 0;
	int zNum = 0;
	Pincell* pinCells = nullptr;

	H Assembly() = default;

	H void Initialize(std::string assemblyTxt, double cellHeight = 0.0, vec3 startPos = { 0, 0, 0 }, vec3 endPos = { 0, 0, 0 }) {

		std::ifstream assembly(assemblyTxt);

		if (!assembly.is_open()) {
			throw std::runtime_error("Cannot open: " + assemblyTxt);
		}

		if (!(assembly >> this->xNum >> this->yNum >> this->zNum)) {
			throw std::runtime_error("Failed to read xNum yNum zNum from: " + assemblyTxt);
		}

		//assembly >> this->xNum >> this->yNum >> this->zNum;
		int numPincells = this->xNum * this->yNum * this->zNum;

		double pincellLength = 0.0;
		double pincellHeight = 0.0;
		double radius = 0.0;

		assembly >> pincellLength >> pincellHeight >> radius;

		if (cellHeight != 0.0) {
			pincellHeight = cellHeight;
		}
		int mod = 0;
		assembly >> mod;

		MatType modType = (mod >= 0 && mod <= 6) ? static_cast<MatType>(mod) : MatType::MOD;
		
		this->pinCells = new Pincell[numPincells];

		for (int i = 0; i < xNum * yNum; i++) {
			int meat = 0;
			assembly >> meat;
			MatType meatType = (meat >= 0 && meat <= 6) ? static_cast<MatType>(meat) : MatType::MOD;
			for (int j = 0; j < zNum; j++) {
				int index = j * (xNum * yNum) + i;
				this->pinCells[index] = Pincell(pincellLength, radius, pincellHeight, meatType, modType);
				//std::cout << static_cast<int>(this->pinCells[j].meatType) << " ";
			}
			//std::cout << 
		}

		assembly.close();
	}
	
	H void Initialize_MOD(double cellSize, vec3 startPos = { 0, 0, 0 }, vec3 endPos = { 0, 0, 0 }, double cellHeight = 0.0) {
		vec3 length = endPos - startPos;

		double height = 0.0;
		if (cellHeight == 0.0) {
			height = cellSize;
		}
		else { 
			height = cellHeight;
		}

		this->xNum = static_cast<int>(length.x / cellSize);
		this->yNum = static_cast<int>(length.y / cellSize);
		this->zNum = static_cast<int>(length.z / height);

		int numPincells = this->xNum * this->yNum * this->zNum;
		this->pinCells = new Pincell[numPincells];

		for (int i = 0; i < this->xNum * this->yNum; i++) {
			for (int j = 0; j < this->zNum; j++) {
				int index = j * (this->xNum * this->yNum) + i;
				this->pinCells[index] = Pincell(cellSize, 0.0, height, MatType::MOD, MatType::MOD);
				//std::cout << static_cast<int>(this->pinCells[j].meatType) << " ";
			}
		}

	}
	
	HD Pincell& returnPincellByIndex(int x, int y, int z);

	HD int totalPincellNo();

	HD Pincell& returnPincellByPos(Neutron n);
};

class C5G7Geometry {
public:
	double x; // cm
	double y;
	double z;

	//PinCell* pinCells = nullptr;
	Assembly* assembly = nullptr;
	Assembly nullAssembly;

	int assemblyNo = 0;

	// I have to put something here --- ahhh

	C5G7Geometry()
		: x(0.0), y(0.0), z(0.0)
	{
	}

	// lazy initialization - we are not going to fully construct here. Instead, we will use the C5G7GeometryFactory to initialize.
	C5G7Geometry(std::string CoreProfile)
	{
		std::ifstream coreProfile(CoreProfile);
		x = 0.6426;
		y = 0.6426;
		z = 2.1420;
	}


	HD Assembly& returnAssemblyByPos(Neutron& n) {
		for (int i = 0; i < this->assemblyNo; i++) {
			vec3 endPos = this->assembly[i].startPos + this->assembly[i].length;
			if (n.pos.x >= this->assembly[i].startPos.x && n.pos.x < endPos.x) {
				if (n.pos.y >= this->assembly[i].startPos.y && n.pos.y < endPos.y) {
					if (n.pos.z >= this->assembly[i].startPos.z && n.pos.z < endPos.z) {
						return this->assembly[i];
					}
				}
			}
		}

		// this usually means fucked up - 
		// we are never meant to pass the out-of-bounds neutrons, or nullified neutrons in this function.
		// we never want the code to flow into this far, in this function - idK why it ended up here, maybe put some debugger outputs just in case
		n.Nullify();
		printf("Neutron Out-Of-Bounds in Position: (%f, %f, %f)\n", n.pos.x, n.pos.y, n.pos.z);
		// return null assembly 
		return this->nullAssembly;
	}
};


H static void fuelLayoutDebug(Assembly& assembly) {
	std::cout << "\nAssembly layout info: total number of pincells in this assembly: " << assembly.totalPincellNo() << "\n";
	std::cout << "for the z=0 layer:\n";
	for (int k = 0; k < 1; k++) {
		for (int ii = 0; ii < assembly.yNum; ii++) {
			for (int ij = 0; ij < assembly.xNum; ij++) {
				int index = (ii * assembly.xNum + ij);
				std::cout << static_cast<int>(assembly.pinCells[index].meatType) << " ";
			}
			std::cout << "\n";
		}
	}
}


class C5G7GeometryFactory {
public:
	H static void Initialize(C5G7Geometry& Core, std::string totalCoreProfile, std::string UO2Geometry, std::string MOXGeometry, double cellHeight = 0.0) {

		double pincellSize = 0.0;
		
		std::ifstream totalCore(totalCoreProfile);
		std::ifstream lineFetch(totalCoreProfile);
		std::string line;
		int lineNo = 0;
		while (std::getline(lineFetch, line)) {
			bool isBlank = std::all_of(line.begin(), line.end(),
				[](unsigned char c) { return std::isspace(c); });
			if (!isBlank) lineNo++;
		}
		lineFetch.close();

		Assembly asmUO2{};
		asmUO2.Initialize(UO2Geometry);
		Assembly asmMOX{};
		asmMOX.Initialize(MOXGeometry);
		Assembly asmMOD{};
		asmMOD.Initialize_MOD(pincellSize);

		totalCore >> Core.x >> Core.y >> Core.z;
		totalCore >> pincellSize;

		std::vector<Assembly> assemblyVec;
		//std::string line;

		Core.assembly = new Assembly[lineNo - 2];

		for (int i = 0; i < lineNo - 2; i++) {
			assemblyVec.reserve(1);
			
			Assembly assembly{};
			std::string assemblyType;

			totalCore >> assemblyType;

			double startPosX = 0.0; double endPosX = 0.0;
			double startPosY = 0.0; double endPosY = 0.0;
			double startPosZ = 0.0; double endPosZ = 0.0;
			totalCore >> startPosX >> endPosX >> startPosY >> endPosY >> startPosZ >> endPosZ;

			double lengthX = endPosX - startPosX;
			double lengthY = endPosY - startPosY;
			double lengthZ = endPosZ - startPosZ;
			vec3 startPos{startPosX, startPosY, startPosZ};
			vec3 endPos(endPosX, endPosY, endPosZ);
			vec3 length{lengthX, lengthY, lengthZ};

			if (assemblyType == "UO2") {
				assemblyVec.emplace_back(asmUO2);
				assemblyVec[i].startPos = startPos;
				assemblyVec[i].length = length;
			}

			if (assemblyType == "MOX") {
				assemblyVec.emplace_back(asmMOX);
				assemblyVec[i].startPos = startPos;
				assemblyVec[i].length = length;
			}

			if (assemblyType == "MOD") {
				assemblyVec.emplace_back(asmMOD);
				assemblyVec[i].Initialize_MOD(pincellSize, startPos, endPos);
			}
			//fuelLayoutDebug(assemblyVec[i]);

		} 
		std::cout << "line size: " << lineNo - 2 << ", vector size: " << assemblyVec.size() << "\n";

		
		if (lineNo - 2 != assemblyVec.size()) {
			delete[] Core.assembly;
			Core.assembly = new Assembly[assemblyVec.size()];
		}

		std::copy(assemblyVec.begin(), assemblyVec.end(), Core.assembly);

		Core.assemblyNo = assemblyVec.size();

		totalCore.close();
	}
	
};