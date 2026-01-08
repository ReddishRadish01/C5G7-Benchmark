#include "cudaHeader.cuh"

#include "XSParser.cuh"
#include "Neutron.cuh"
#include "Core.cuh"

HD Pincell& Assembly::returnPincellByIndex(int x, int y, int z) {
	if (x >= this->xNum || y >= this->yNum || z >= this->zNum) {
		printf("index [%d][%d][%d] out of bounds", x, y, z);
		return this->pinCells[0];
	}
	
	int index = z * (this->xNum * this->yNum) + (this->xNum * y) + x;
	return this->pinCells[index];
}

HD int Assembly::totalPincellNo() {
	return this->xNum * this->yNum * this->zNum;
}

HD Pincell& Assembly::returnPincellByPos(Neutron n) {
	vec3 neutronPos = n.pos;
	neutronPos

	
	this.returnByIndex

	return this->pinCells[0];
}