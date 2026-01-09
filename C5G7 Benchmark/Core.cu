#include "cudaHeader.cuh"

#include "XSParser.cuh"
#include "Neutron.cuh"
#include "Core.cuh"

HD double Pincell::DTC(vec3 localPos, Neutron& n, XSLibrary& XSLib, unsigned long long xi) {
	GnuAMCM RNG(xi);
	MatType matType = this->meatOrMod(localPos);

	double transXS = XSLib.returnMatByType(matType).transXS[static_cast<int>(n.energy)];
	return RNG.uniform_open(0.0, 1.0) / transXS;
}


HD double Pincell::DTS(vec3 localPos, Neutron& n) {


	return 0;
}


HD Pincell& Assembly::returnPincellByIndex(int x, int y, int z) {
	if (x >= this->xNum || y >= this->yNum || z >= this->zNum) {
		printf("index [%d][%d][%d] out of bounds\n", x, y, z);
		Pincell OOBPincell = Pincell(0, 0, 0);
		return OOBPincell;
		//return this->pinCells[0];
	}
	
	int index = z * (this->xNum * this->yNum) + (this->xNum * y) + x;
	return this->pinCells[index];
}

HD int Assembly::totalPincellNo() {
	return this->xNum * this->yNum * this->zNum;
}

HD Pincell& Assembly::returnPincellByPos(Neutron n) {
	vec3 localAssemblyPos = n.pos - this->startPos;
	double cellSideLength = this->pinCells[0].sideLength;
	double cellHeight = this->pinCells[0].height;
	int xIdx = static_cast<int>(localAssemblyPos.x / cellSideLength);
	int yIdx = static_cast<int>(localAssemblyPos.y / cellSideLength);
	int zIdx = static_cast<int>(localAssemblyPos.z / cellHeight);

	return this->returnPincellByIndex(xIdx, yIdx, zIdx);

}

HD vec3 Assembly::returnFlooredNeutronPosInPincell(Neutron& n) {
	vec3 localAssemblyPos = n.pos - this->startPos;
	int xIdx = localAssemblyPos.x / (this->length.x / this->xNum);
	int yIdx = localAssemblyPos.y / (this->length.y / this->yNum);
	int zIdx = localAssemblyPos.z / (this->length.z / this->zNum);

	return { localAssemblyPos.x - (this->xNum * xIdx),
			 localAssemblyPos.y - (this->yNum * yIdx),
			 localAssemblyPos.z - (this->zNum * zIdx) };
}


HD double Assembly::DTC(Neutron& n, XSLibrary& XSLib, unsigned long long xi) {
	GnuAMCM RNG(xi);
	Pincell currentPincell = this->returnPincellByPos(n);
	vec3 pincellLocalPos = this->returnFlooredNeutronPosInPincell(n);
	//MatType mat = currentPincell.meatOrMod(pincellLocalPos);
	return currentPincell.DTC(pincellLocalPos, n, XSLib, RNG.gen());
}

HD double Assembly::DTS(Neutron& n) {
	Pincell currentPincell = this->returnPincellByPos(n);
	vec3 pincellLocalPos = this->returnFlooredNeutronPosInPincell(n);

	return currentPincell.DTS(pincellLocalPos, n);
}