#pragma once

#include "cudaHeader.cuh"
#include "XSParser.cuh"
#include "Neutron.cuh"
#include "Core.cuh"
#include "CoreManager.cuh"


// this is for finding XS from corresponding location
class XSManager {
public:
	HD static double returnXSByMat(MatType matType, MatXS& matXS, InteractionType interactionType, double in, double out = 0.0) {
		//in = static_cast<int>(in);
		//out = static_cast<int>(out);
		switch(interactionType) {
		case InteractionType::ntot: return matXS.returnMatByType(matType).returnXSbyType(XSType::trans, in, out);

		
		}
		
		
	}

};