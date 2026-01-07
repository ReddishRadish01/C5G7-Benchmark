#pragma once
#include "cudaHeader.cuh"
#include "XSParser.cuh"


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

	H void Initialize(std::string assemblyTxt, vec3 startPos = { 0, 0, 0 }, vec3 endPos = { 0, 0, 0 }) {

		std::ifstream assembly(assemblyTxt);

		if (!assembly.is_open()) {
			throw std::runtime_error("Cannot open: " + assemblyTxt);
		}

		if (!(assembly >> xNum >> yNum >> zNum)) {
			throw std::runtime_error("Failed to read xNum yNum zNum from: " + assemblyTxt);
		}

		//assembly >> this->xNum >> this->yNum >> this->zNum;
		int numPincells = xNum * yNum * zNum;

		double pincellLength = 0.0;
		double pincellHeight = 0.0;
		double radius = 0.0;

		assembly >> pincellLength >> pincellHeight >> radius;

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
	
	H void Initialize_MOD(vec3 startPos = {0, 0, 0}, vec3 endPos = {0, 0, 0}) {
			


	}

	HD Pincell& returnByIndex(int x, int y, int z) {
		if (x * y * z > this->xNum * this->yNum * this->zNum) {
			return this->pinCells[0];
		}

		int index = z * (this->xNum * y + x);
		return this->pinCells[index];
	}
	
	HD int pincellNo() {
		return this->xNum * this->yNum * this->zNum;
	}
};

class C5G7Geometry {
public:
	double x; // cm
	double y;
	double z;

	//PinCell* pinCells = nullptr;
	Assembly* assembly = nullptr;

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

};

class C5G7GeometryFactory {

	H static void Initialize(C5G7Geometry& Core, std::string totalCoreProfile, std::string UO2Geometry, std::string MOXGeometry) {

		double pincellSize = 0.0;
		//Core.pinCells = new PinCell[100];


		std::ifstream totalCore(totalCoreProfile);
		std::ifstream lineFetch(totalCoreProfile);
		std::string line;
		int lineNo = 0;
		while (std::getline(lineFetch, line)) {
			bool isBlank = std::all_of(line.begin(), line.end(),
				[](unsigned char c) { return std::isspace(c); });
			if (!isBlank) lineNo++;
		}

		Assembly asmUO2{};
		asmUO2.Initialize(UO2Geometry);
		Assembly asmMOX{};
		asmMOX.Initialize(MOXGeometry);
		Assembly asmMOD{};
		asmMOD.Initialize_MOD();



		totalCore >> Core.x >> Core.y >> Core.z;
		totalCore >> pincellSize;

		std::vector<Assembly> assemblyVec;
		//std::string line;

		
		for (int i = 0; i < lineNo - 2; i++) {
			assemblyVec.reserve(1);
			
			Assembly assembly{};
			vec3 startpos{};
			vec3 endpos{};

			std::string assemblyType;

			totalCore >> assemblyType;
			if (assemblyType == "UO2") {

				//assembly
			}
			




		} 





				// 1. read Total Core principle size.
				// 2. read UO2 Part - put it there. Preferably reuse it.
				// 3. read MOX Part - reuse it.
				// 4. fill the rest part with moderator: 
				// 5. set boundary condition? how?


		totalCore.close();
	}
	
};