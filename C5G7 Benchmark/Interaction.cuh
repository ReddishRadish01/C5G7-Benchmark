#pragma once

#include "cudaHeader.cuh"
#include "XSParser.cuh"
#include "Neutron.cuh"
#include "Core.cuh"
#include "XSManager.cuh"
#include "CoreManager.cuh"

//#define INTERACTIONDEBUG 
//#define FISSIONPRINT



class Interaction {
public:
	HD static void scatter(Neutron& neutron, double outEnergy, GnuAMCM& RNG) {
		//neutron.updateWithLength(1.0e-10);
		neutron.energy = outEnergy;
		neutron.dirVec = vec3::randomUnit(RNG);
	}

	HD static void absorption(Neutron& n) {
		n.Nullify();
	}

	D static void fission(Neutron& n, NeutronBank* Bank, MatXS& matXS, GnuAMCM& RNG, double* k_mult, bool passFlag) {
		double nu = matXS.nu[static_cast<int>(n.energy) - 1];
		int fissionNum = static_cast<int>(nu / *k_mult + RNG.uniform(0.0, 1.0));

		n.dirVec = vec3::randomUnit(RNG);
		int fissionE = XSManager::returnFissionNeutronEnergy(matXS, RNG);
		int addIndex = atomicAdd(&(Bank->addedNeutronIndex), fissionNum - 1);
		atomicAdd(&(Bank->addedNeutronSize), fissionNum - 1);
		n.energy = fissionE;

		for (int i = 0; i < fissionNum - 1; i++) {
			//Bank->addedNeutrons[addIndex + i].status = true;
			int fissionE = XSManager::returnFissionNeutronEnergy(matXS, RNG);
			Bank->addedNeutrons[addIndex + i].reInitialize(n.pos, vec3::randomUnit(RNG), fissionE, 1.0, passFlag);
		}

	}

	D static void reaction(Neutron& n, NeutronBank* Bank, XSLibrary* XSLib, Pincell currentPincell, vec3 localPos, GnuAMCM& RNG, double* k_mult, bool passFlag, bool add) {
		MatType currentMat = currentPincell.meatOrMod(localPos);
		double outEnergy = 0.0;
		InteractionType interactionT = XSManager::returnInteracitonType(XSLib, currentMat, RNG, n.energy, outEnergy);
		if (interactionT == InteractionType::nel) {
			Interaction::scatter(n, outEnergy, RNG);
		}
		else if (interactionT == InteractionType::ng) {
			Interaction::absorption(n);
			if (add == true) { // absorption in addedneutron
				atomicAdd(&(Bank->addedNeutronSize), -1);
			}
			else {
				atomicAdd(&(Bank->neutronSize), -1);
			}
		}
		else if (interactionT == InteractionType::nf) {
			Interaction::fission(n, Bank, XSLib->returnMatByType(currentMat), RNG, k_mult, passFlag);
		}
	}

	D static void reflection(Neutron& n, double DTS, vec3 updatedSurfacePos, double eps) {
		vec3 reflectNormal = { 0.0, 0.0, 0.0 };
		if (updatedSurfacePos.x <= eps) {
			reflectNormal = { 1.0, 0.0, 0.0 };
		}
		else if (updatedSurfacePos.y <= eps) {
			reflectNormal = { 0.0, 1.0, 0.0 };
		}
		else if (updatedSurfacePos.z <= eps) {
			reflectNormal = { 0.0, 0.0, 1.0 };
		}
		else {
			// this is fucked
			printf("Error - reflectnormal not set\n");
			n.Nullify();
			return;
		}

		vec3 collisionPos = n.pos + n.dirVec * DTS;
		vec3 reflectVec = n.dirVec - reflectNormal * 2 * n.dirVec.dot(reflectNormal);
		collisionPos = collisionPos + reflectVec * eps * 1000;

		n.pos = collisionPos;
		n.dirVec = reflectVec;
	}

	H static void reaction_CPU(int idx, Neutron& n, NeutronBank* Bank, XSLibrary* XSLib, Pincell currentPincell, vec3 localPos, GnuAMCM& RNG, double* k_mult, bool passFlag, bool add, double& absorption, double& fission) {
		MatType currentMat = currentPincell.meatOrMod(localPos);
		double outEnergy = 0.0;
		InteractionType interactionT = XSManager::returnInteracitonType(XSLib, currentMat, RNG, n.energy, outEnergy);
		if (interactionT == InteractionType::nel) {
#ifdef INTERACTIONDEBUG 
			printf("idx %d neutron pos (%f,%f,%f) - n,el reaction, scattered from %1.0f to %1.0f\n", idx, n.pos.x, n.pos.y, n.pos.z, n.energy, outEnergy);
#endif
			Interaction::scatter(n, outEnergy, RNG);
		}
		else if (interactionT == InteractionType::ng) {
#ifdef INTERACTIONDEBUG 
			printf("idx %d neutron pod (%f,%f,%f) - n,g reaction\n", idx, n.pos.x, n.pos.y, n.pos.z);
#endif
			Interaction::absorption(n);
			absorption++;
			if (add == true) { // absorption in addedneutron
				Bank->addedNeutronSize -= 1;
			}
			else {
				Bank->neutronSize -= 1;
			}
		}
		else if (interactionT == InteractionType::nf) {
			// fission - function is written inside here (inlined)
			double nu = XSLib->returnMatByType(currentPincell.meatOrMod(localPos)).nu[static_cast<int>(n.energy) - 1];
			int fissionNum = static_cast<int>(nu / *k_mult + RNG.uniform(0.0, 1.0));
			
			n.dirVec = vec3::randomUnit(RNG);
			int fissionE = XSManager::returnFissionNeutronEnergy(XSLib->returnMatByType(currentPincell.meatOrMod(localPos)), RNG);
			n.energy = fissionE;
			//n.updateWithLength(1.0e-10);
			int addIndex = Bank->addedNeutronIndex;
			Bank->addedNeutronIndex += fissionNum-1;
			Bank->addedNeutronSize += fissionNum-1;
#ifdef FISSIONPRINT 
			printf("Fission on %d, pos: (%f,%f,%f), mat: %s, fission N num: %d, added to index %d.\n", idx, n.pos.x, n.pos.y, n.pos.z, to_string(currentPincell.meatOrMod(localPos)), fissionNum, addIndex);
#endif
			fission += fissionNum;
			for (int i = 0; i < fissionNum - 1; i++) {
				//Bank->addedNeutrons[addIndex + i].status = true;
				int fissionE = XSManager::returnFissionNeutronEnergy(XSLib->returnMatByType(currentPincell.meatOrMod(localPos)), RNG);
				Bank->addedNeutrons[addIndex + i].reInitialize(n.pos, vec3::randomUnit(RNG), fissionE, 1.0, passFlag);
			}
		}
	}

	H static void reflection_CPU(Neutron& n, double DTS, vec3 updatedSurfacePos, double eps) {
		vec3 reflectNormal = { 0.0, 0.0, 0.0 };
		if (updatedSurfacePos.x <= eps) {
			reflectNormal = { 1.0, 0.0, 0.0 };
		}
		else if (updatedSurfacePos.y <= eps) {
			reflectNormal = { 0.0, 1.0, 0.0 };
		}
		else if (updatedSurfacePos.z <= eps) {
			reflectNormal = { 0.0, 0.0, 1.0 };
		}
		else {
			// this is fucked
			printf("Error - reflectnormal not set\n");
			n.Nullify();
			return;
		}

		//vec3 collisionPos = n.pos + n.dirVec * DTS;
		vec3 collisionPos = updatedSurfacePos;
		vec3 reflectVec = n.dirVec - reflectNormal * (2 * n.dirVec.dot(reflectNormal));
		collisionPos = collisionPos + reflectVec * eps * 1000;

		n.pos = collisionPos;
		n.dirVec = reflectVec;

	}
};