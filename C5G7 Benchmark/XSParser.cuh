#pragma once

#include "cudaHeader.cuh"

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

enum class MatType {
	UO2,
	MOX4_3,
	MOX7_0,
	MOX8_7,
	GT,
	FC,
	MOD,
	Unknown
};


static inline std::vector<double> parse_double_from_line(const std::string& line) {
	std::istringstream iss(line);
	std::vector<double> v;
	double x;
	while (iss >> x) v.push_back(x);
	return v;
}


struct G7 {
	MatType matType{};
	double totalXS[7]{};
	double transXS[7]{};
	double absXS[7]{};
	double capXS[7]{};
	double fisXS[7]{};
	double nu[7]{};
	double chi[7]{};

	double elsXS[7][7]{};

	// Default constructor
	H G7() = default;
	H G7(MatType matType) : matType(matType) {}

	H G7(std::string fileName, MatType matType)
		: matType(matType)
	{
		std::ifstream C5Data(fileName);
		if (C5Data.fail()) {
			std::cout << "Error! File name " << fileName << " doesn't Exist, or Intentionally omitted\n";
		}
		else { std::cout << "Loading C5Data from " << fileName << " Complete!\n"; }
		
		std::string dummy;
		std::getline(C5Data, dummy);

		for (int g = 0; g < 7; g++) {
			std::string line;
			do {
				if (!std::getline(C5Data, line)) {
					throw std::runtime_error("Unexpected EOF wile reading XS rows");
				}
			} while (line.find_first_not_of(" \t\r\n") == std::string::npos);

			std::vector<double> vals = parse_double_from_line(line);
			if (vals.size() == 7) {
				totalXS[g] = vals[0];
				transXS[g] = vals[1];
				absXS[g] = vals[2];
				capXS[g] = vals[3];
				fisXS[g] = vals[4];
				nu[g] = vals[5];
				chi[g] = vals[6];
			}
			else {
				totalXS[g] = vals[0];
				transXS[g] = vals[1];
				absXS[g] = vals[2];
				capXS[g] = vals[3];
				fisXS[g] = 0;
				nu[g] = 0;
				chi[g] = 0;
			}

		}

		for (int i = 0; i < 7; i++) {
			for (int j = 0; j < 7; j++) {
				C5Data >> elsXS[i][j];

			}
		}
		
	}

	HD double returnXSbyType(XSType xsType, double currentEnergy, double destination = 0) const {
		switch (xsType) {
		case XSType::tot:		return this->totalXS[static_cast<int>(currentEnergy)];
		case XSType::trans:		return this->transXS[static_cast<int>(currentEnergy)];
		case XSType::abs:		return this->absXS[static_cast<int>(currentEnergy)];
		case XSType::cap:		return this->capXS[static_cast<int>(currentEnergy)];
		case XSType::fis:		return this->fisXS[static_cast<int>(currentEnergy)];
		case XSType::nu:		return this->nu[static_cast<int>(currentEnergy)];
		case XSType::chi:		return this->chi[static_cast<int>(currentEnergy)];
		case XSType::el:		return this->elsXS[static_cast<int>(currentEnergy)][static_cast<int>(destination)];
		}
	}

	H void g7DeviceAllocator(G7*& d_G7);

};


class MatXS {

public:
	G7 UO2{ MatType::UO2 };
	G7 MOX4_3{ MatType::MOX4_3 };
	G7 MOX7_0{ MatType::MOX7_0 };
	G7 MOX8_7{ MatType::MOX8_7 };
	G7 GT{ MatType::GT };
	G7 FC{ MatType::FC };
	G7 MOD{ MatType::MOD };

	MatXS() = default;


	// reference return needed - see the initialize of the MatXSFactory. It requires a reference return.
	HD G7& returnMatByType(MatType matType) {
		switch (matType) {
		case MatType::UO2:		return UO2;
		case MatType::MOX4_3:	return MOX4_3;
		case MatType::MOX7_0:	return MOX7_0;
		case MatType::MOX8_7:	return MOX8_7;
		case MatType::GT:		return GT;
		case MatType::FC:		return FC;
		case MatType::MOD:		return MOD;
		default:				return MOD;
		}
	}

	/*
	// separated device proprietary codes
	D G7 returnMatByType(MatType matType) const {
		switch (matType) {
		case MatType::UO2:		return UO2;
		case MatType::MOX4_3:	return MOX4_3;
		case MatType::MOX7_0:	return MOX7_0;
		case MatType::MOX8_7:	return MOX8_7;
		case MatType::FC:		return FC;
		case MatType::GT:		return GT;
		case MatType::MOD:		return MOD;
		default:				return MOD;
		}
	}
	*/
	
};

class MatXSFactory {

public:
	MatXS matXS;

	H static void initialize(MatXS& matXS, std::vector<G7>& matXS_vec) {
		for (G7& XS : matXS_vec) {
			matXS.returnMatByType(XS.matType) = XS;
		}
	}

	H static void initialize(std::vector<G7>& matXS_vec) {
		
	}


};

