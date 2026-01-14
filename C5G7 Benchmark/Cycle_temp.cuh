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

//#define NEUTRONIDXDEBUG
//#define REFLECTIONDEBUG
#define OUTBOUNDDEBUG

constexpr double globaleps = 1.0e-12;

G void cycle_Neutron(NeutronBank* Bank, C5G7Geometry* Core, XSLibrary* XSLib, unsigned long long* seedNo, double* k_mult, bool passFlag) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= Bank->allocatableNeutronNum) {
        return;
    }
    GnuAMCM RNG(seedNo[idx]);
    Neutron localNeutron = Bank->neutrons[idx];
    double eps = 1.0e-10;

    if (!localNeutron.isNullified()) {
        for (int i = 0; i < 100; i++) {
            //vec3 flooredNeutronPos = Core->assembly->returnFlooredNeutronPosInPincell(localNeutron);
            Assembly currentAssembly = Core->returnAssemblyByNeutron(localNeutron);
            Pincell currentPincell = currentAssembly.returnPincellByPos(localNeutron);

            double DTC = currentAssembly.DTC(localNeutron, XSLib, RNG);
            double DTS = currentAssembly.DTS(localNeutron);

            if (DTC < DTS) {    // reaction
                localNeutron.updateWithLength(DTC);
                vec3 flooredNeutronPos = Core->assembly->returnFlooredNeutronPosInPincell(localNeutron);
                Interaction::reaction(localNeutron, Bank, XSLib, currentPincell, flooredNeutronPos, RNG, k_mult, passFlag, false);
                seedNo[idx] = RNG.gen();
                return;
            }
            else {  // do:  boundary check / reflection / position update to DTS and feed it back to main loop
                vec3 updatedPos = localNeutron.pos + localNeutron.dirVec * DTC;
                vec3 updatedSurfacePos = localNeutron.pos + localNeutron.dirVec * DTS;
                // handle vaccum boundary neutrons;
                if (updatedPos.x >= Core->x || updatedPos.y >= Core->y || updatedPos.z >= Core->z) {
                    localNeutron.Nullify();
                    atomicAdd(&(Bank->neutronSize), -1);
                    seedNo[idx] = RNG.gen();
                    return;
                }
                if (updatedSurfacePos.x <= eps || updatedSurfacePos.y <= eps || updatedSurfacePos.z <= eps) {
                    Interaction::reflection(localNeutron, DTS, updatedSurfacePos, eps);
                    if (localNeutron.isNullified()) {
                        // error handling - reflection might return a nullified neutron - thus reduce the size of neutron by 1.
                        atomicAdd(&(Bank->neutronSize), -1);
                        seedNo[idx] = RNG.gen();
                        return;
                    }
                    Bank->neutrons[idx].updateWithLength(DTS * (1.0 + globaleps));
                    continue;
                }
                // neutron is inside the boundary.
                localNeutron.updateWithLength(DTS + eps);
            }
            printf("idx %d neutron on %d loop\n", idx, i);
        }
        //printf("idx %d neutron didn't reacted after 100 loops..?\n", idx);
        return;
    }
    return;

}

G void addedNeutronPassResetter(NeutronBank* Neutrons) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= Neutrons->allocatableNeutronNum) { return; }

    if (!Neutrons->addedNeutrons[idx].isNullified()) {
        Neutrons->addedNeutrons[idx].passFlag = false;
    }
}


G void cycle_addedNeutron(NeutronBank* Bank, C5G7Geometry* Core, XSLibrary* XSLib, unsigned long long* seedNo, double* k_mult, bool passFlag) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= Bank->allocatableNeutronNum) {
        return;
    }
    GnuAMCM RNG(seedNo[idx]);
    Neutron localNeutron = Bank->addedNeutrons[idx];
    double eps = 1.0e-10;

    if (!localNeutron.isNullified()) {
        for (int i = 0; i < 100; i++) {
            //vec3 flooredNeutronPos = Core->assembly->returnFlooredNeutronPosInPincell(localNeutron);
            Assembly currentAssembly = Core->returnAssemblyByNeutron(localNeutron);
            Pincell currentPincell = currentAssembly.returnPincellByPos(localNeutron);

            double DTC = currentAssembly.DTC(localNeutron, XSLib, RNG);
            double DTS = currentAssembly.DTS(localNeutron);

            if (DTC < DTS) {    // reaction
                localNeutron.updateWithLength(DTC);
                vec3 flooredNeutronPos = Core->assembly->returnFlooredNeutronPosInPincell(localNeutron);
                Interaction::reaction(localNeutron, Bank, XSLib, currentPincell, flooredNeutronPos, RNG, k_mult, passFlag, true);
                return;
            }
            else {  // do:  boundary check / reflection / position update to DTS and feed it back to main loop
                vec3 updatedPos = localNeutron.pos + localNeutron.dirVec * DTC;
                vec3 updatedSurfacePos = localNeutron.pos + localNeutron.dirVec * DTS;
                // handle vaccum boundary neutrons;
                if (updatedPos.x >= Core->x || updatedPos.y >= Core->y || updatedPos.z >= Core->z) {
                    localNeutron.Nullify();
                    atomicAdd(&(Bank->neutronSize), -1);
                    return;
                }
                if (updatedSurfacePos.x <= eps || updatedSurfacePos.y <= eps || updatedSurfacePos.z <= eps) {
                    Interaction::reflection(localNeutron, DTS, updatedSurfacePos, eps);
                    if (localNeutron.isNullified()) {
                        // error handling - reflection might return a nullified neutron - thus reduce the size of neutron by 1.
                        atomicAdd(&(Bank->addedNeutronSize), -1);
                    }
                    return;
                }
                // neutron is inside the boundary.
                localNeutron.updateWithLength(DTS * (1 + 1.0e-5));
            }
            printf("idx %d neutron on %d loop\n", idx, i);
        }
        printf("idx %d neutron didn't reacted after 100 loops...?\n", idx);
    }


}