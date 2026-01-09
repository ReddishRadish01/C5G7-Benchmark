#pragma once

#include "cudaHeader.cuh"
#include "XSParser.cuh"
#include "Neutron.cuh"
#include "Core.cuh"
#include "XSManager.cuh"
#include "CoreManager.cuh"


class Interaction {
public:
	HD static InteractionType returnInteraction(Pincell pincell, Neutron n, XSLibrary XSLib, unsigned long long xi) {


		GnuAMCM RNG(xi);

	}

	HD static double DTS(Pincell pincell, Neutron n) {
		
	}

	HD static double DTC(Pincell pincell, Neutron n, XSLibrary XSLib, unsigned long long xi) {
		GnuAMCM RNG(xi);

	} 


};
