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

#define TALLYINCYCLE


constexpr double globaleps = 1.0e-10;


G void cycle_Neutron(NeutronBank* Bank, C5G7Geometry* Core, TallyC5G7Geometry* CoreTally, XSLibrary* XSLib, unsigned long long* seedNo, double* k_mult, bool passFlag, int iterTimes) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= Bank->allocatableNeutronNum) {
        return;
    }
    GnuAMCM RNG(seedNo[idx]);

    double eps = 1.0e-13;
    int reflectionCounter = 0;
    if (!Bank->neutrons[idx].isNullified()) {
        for (int i = 0; i < iterTimes; i++) {
            //vec3 flooredNeutronPos = Core->assembly->returnFlooredNeutronPosInPincell(localNeutron);
            Assembly currentAssembly = Core->returnAssemblyByNeutron(Bank->neutrons[idx]);
            /*
            if (currentAssembly.length.x == 0.0) {
                Bank->neutrons[idx].Nullify();
                atomicAdd(&(Bank->neutronSize), -1);
                seedNo[idx] = RNG.gen();
                return;
            }
            */
#ifdef TALLYINCYCLE
            TallyAssembly currentTallyAssembly = CoreTally->returnAssemblyByNeutron(Bank->neutrons[idx]);
#endif
            
            
            if (Bank->neutrons[idx].isNullified()) {
                atomicAdd(&(Bank->neutronSize), -1);
#ifdef OUTBOUNDDEBUG
                printf("Neutron idx %d nullfied because Assembly positioning got fucked\n", idx);
#endif
                seedNo[idx] = RNG.gen();
                return;
            }
            
            Pincell currentPincell = currentAssembly.returnPincellByPos(Bank->neutrons[idx]);
            vec3 flooredNeutronPos = currentAssembly.returnFlooredNeutronPosInPincell(Bank->neutrons[idx]);
#ifdef TALLYINCYCLE
            TallyPincell currentTallyPincell = currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]);
#endif
            
            if (currentPincell.sideLength == 0.0 && currentPincell.height == 0.0) {
                // on-surface exceptions - bump it into other cell, and continue
                if ((Bank->neutrons[idx].pos.x >= currentAssembly.startPos.x + currentAssembly.length.x - eps * 100) && (Bank->neutrons[idx].pos.x <= currentAssembly.startPos.x + currentAssembly.length.x + eps * 100) ||
                    (Bank->neutrons[idx].pos.y >= currentAssembly.startPos.y + currentAssembly.length.y - eps * 100) && (Bank->neutrons[idx].pos.y <= currentAssembly.startPos.y + currentAssembly.length.y + eps * 100) ||
                    (Bank->neutrons[idx].pos.z >= currentAssembly.startPos.z + currentAssembly.length.z - eps * 100) && (Bank->neutrons[idx].pos.z <= currentAssembly.startPos.z + currentAssembly.length.z + eps * 100) ) {
                    Bank->neutrons[idx].updateWithLength(eps * 10000);
                    continue;
                }
                printf("error in returning pincell of neutron idx %d, - nullifying this neutron\n", idx);
                Bank->neutrons[idx].Nullify();
                atomicAdd(&(Bank->neutronSize), -1);
                seedNo[idx] = RNG.gen();
                return;
            }
            

            double DTC = currentAssembly.DTC(Bank->neutrons[idx], XSLib, RNG);
            double DTS = currentAssembly.DTS(Bank->neutrons[idx]);

            if (DTC < DTS) {    // reaction

                Bank->neutrons[idx].updateWithLength(DTC);
                //vec3 flooredNeutronPos = Core->assembly->returnFlooredNeutronPosInPincell(Bank->neutrons[idx]);
#ifdef TALLYINCYCLE
                if (currentPincell.meatOrMod(flooredNeutronPos) == MatType::MOD) { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]).modTally), DTC); }
                else { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]).pinTally), DTC); }
#endif
                Interaction::reaction(Bank->neutrons[idx], Bank, XSLib, currentPincell, flooredNeutronPos, RNG, k_mult, passFlag, false);
                seedNo[idx] = RNG.gen();
                return;
            }
            else {  // do:  boundary check / reflection / position update to DTS and feed it back to main loop
                //vec3 updatedPos = Bank->neutrons[idx].pos + Bank->neutrons[idx].dirVec * DTC;
#ifdef TALLYINCYCLE
                if (currentPincell.meatOrMod(flooredNeutronPos) == MatType::MOD) { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]).modTally), DTS); }
                else { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->neutrons[idx]).pinTally), DTS); }
#endif
                vec3 updatedSurfacePos = Bank->neutrons[idx].pos + Bank->neutrons[idx].dirVec * DTS;
                vec3 afterDTCPos = Bank->neutrons[idx].pos + Bank->neutrons[idx].dirVec * DTC;
                // handle vaccum boundary neutrons;
                if (updatedSurfacePos.x >= Core->x  || updatedSurfacePos.y >= Core->y  || updatedSurfacePos.z >= Core->z) {
                    Bank->neutrons[idx].Nullify();
                    atomicAdd(&(Bank->neutronSize), -1);
                    seedNo[idx] = RNG.gen();
                    return;
                }
                if (afterDTCPos.x < 0.0 || afterDTCPos.y < 0.0  || afterDTCPos.z < 0.0) {
                    if (updatedSurfacePos.x <= eps/100 || updatedSurfacePos.y <= eps/100 || updatedSurfacePos.z <= eps/100 ) {

                        Interaction::reflection(Bank->neutrons[idx], DTS, updatedSurfacePos, eps);
#ifdef REFLECTIONDEBUG
                        reflectionCounter++;
                        if (reflectionCounter > 1) {
                            printf("Neutron reflected twice on ");
                            Bank->neutrons[idx].printInfo_Kernel(idx);
                        }
#endif
                        // this reflection already bumps the neutron to the surface location, with some room for error
                        if (Bank->neutrons[idx].isNullified()) {
                            // error handling - reflection might return a nullified neutron - thus reduce the size of neutron by 1.
                            atomicAdd(&(Bank->neutronSize), -1);
                            seedNo[idx] = RNG.gen();
                            return;
                        }
                        continue;
                    }
                }
                // neutron is inside the core
                Bank->neutrons[idx].updateWithLength(DTS + eps * 1000);
                //Bank->neutrons[idx].updateWithLength(DTS);
                
                //Bank->neutrons[idx].dirVec.x >= 0 ? Bank->neutrons[idx].pos.x += globaleps : Bank->neutrons[idx].pos.x -= globaleps;
                //Bank->neutrons[idx].dirVec.y >= 0 ? Bank->neutrons[idx].pos.y += globaleps : Bank->neutrons[idx].pos.y -= globaleps;
                //Bank->neutrons[idx].dirVec.z >= 0 ? Bank->neutrons[idx].pos.z += globaleps : Bank->neutrons[idx].pos.z -= globaleps;
                // ¾ö¸¶¾ø´Â»õ³¢ ¾¾¹ß ÀÌ°Å Ãß°¡ÇÏ´Ï±î ¿°º´³ª´Âµ¥ ¹¹³Ä ¾¾¹ß·Ã
                /*
                vec3 epsVec = {
                    (Bank->neutrons[idx].dirVec.x >= 0 ? +globaleps : -globaleps),
                    (Bank->neutrons[idx].dirVec.y >= 0 ? +globaleps : -globaleps),
                    (Bank->neutrons[idx].dirVec.z >= 0 ? +globaleps : -globaleps)
                };
                Bank->neutrons[idx].pos = Bank->neutrons[idx].pos + epsVec;
                */
                
            }
            //printf("idx %d neutron on %d loop\n", idx, i);

            if (i == iterTimes - 1) {
                Bank->neutrons[idx].printInfo_Kernel(idx);
                Bank->neutrons[idx].Nullify();
                atomicAdd(&(Bank->neutronSize), -1);
            }
        }

        
        //printf("idx %d neutron didn't reacted after %d loops..?\n", idx, iterTimes);
        seedNo[idx] = RNG.gen();
        return;

    }
    seedNo[idx] = RNG.gen();
    return;

}

G void addedNeutronPassResetter(NeutronBank* Neutrons) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= Neutrons->allocatableNeutronNum) { return; }

    if (!Neutrons->addedNeutrons[idx].isNullified()) {
        Neutrons->addedNeutrons[idx].passFlag = false;
    }
}


G void cycle_addedNeutron(NeutronBank* Bank, C5G7Geometry* Core, TallyC5G7Geometry* CoreTally, XSLibrary* XSLib, unsigned long long* seedNo, double* k_mult, bool passFlag, int iterTimes) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx >= Bank->allocatableNeutronNum) {
        return;
    }
    GnuAMCM RNG(seedNo[idx]);
    double eps = 1.0e-13;
    int reflectionCounter = 0;

    if (!Bank->addedNeutrons[idx].isNullified()) {
        if (Bank->addedNeutrons[idx].passFlag) {
            Bank->addedNeutrons[idx].passFlag = false;
            return;
        }
    }

    if (!Bank->addedNeutrons[idx].isNullified()) {
        for (int i = 0; i < iterTimes; i++) {
            //vec3 flooredNeutronPos = Core->assembly->returnFlooredNeutronPosInPincell(localNeutron);
            Assembly currentAssembly = Core->returnAssemblyByNeutron(Bank->addedNeutrons[idx]);
#ifdef TALLYINCYCLE
            TallyAssembly currentTallyAssembly = CoreTally->returnAssemblyByNeutron(Bank->addedNeutrons[idx]);
#endif
            /*
            if (currentAssembly.length.x == 0.0) {
                Bank->addedNeutrons[idx].Nullify();
                atomicAdd(&(Bank->addedNeutronSize), -1);
                seedNo[idx] = RNG.gen();
                return;
            }
            */

            if (Bank->addedNeutrons[idx].isNullified()) {
                atomicAdd(&(Bank->addedNeutronSize), -1);
#ifdef OUTBOUNDDEBUG
                printf("Neutron idx %d nullfied because Assembly positioning got fucked\n", idx);
#endif
                seedNo[idx] = RNG.gen();
                return;
            }
            
            Pincell currentPincell = currentAssembly.returnPincellByPos(Bank->addedNeutrons[idx]);
            vec3 flooredAddedNeutronPos = currentAssembly.returnFlooredNeutronPosInPincell(Bank->addedNeutrons[idx]);
            
            if (currentPincell.sideLength == 0.0 && currentPincell.height == 0.0) {
                if ((Bank->addedNeutrons[idx].pos.x >= currentAssembly.startPos.x + currentAssembly.length.x - eps * 100) && (Bank->addedNeutrons[idx].pos.x <= currentAssembly.startPos.x + currentAssembly.length.x + eps * 100) ||
                    (Bank->addedNeutrons[idx].pos.y >= currentAssembly.startPos.y + currentAssembly.length.y - eps * 100) && (Bank->addedNeutrons[idx].pos.y <= currentAssembly.startPos.y + currentAssembly.length.y + eps * 100) ||
                    (Bank->addedNeutrons[idx].pos.z >= currentAssembly.startPos.z + currentAssembly.length.z - eps * 100) && (Bank->addedNeutrons[idx].pos.z <= currentAssembly.startPos.z + currentAssembly.length.z + eps * 100)) {
                    Bank->addedNeutrons[idx].updateWithLength(eps * 10000);
                    continue;
                }


                printf("error in returning pincell of adddedNeutron idx %d, - nullifying this neutron\n", idx);
                Bank->addedNeutrons[idx].Nullify();
                atomicAdd(&(Bank->addedNeutronSize), -1);
                seedNo[idx] = RNG.gen();
                return;
            }
            
            double DTC = currentAssembly.DTC(Bank->addedNeutrons[idx], XSLib, RNG);
            double DTS = currentAssembly.DTS(Bank->addedNeutrons[idx]);

            if (DTC < DTS) {    // reaction
                Bank->addedNeutrons[idx].updateWithLength(DTC);
                //vec3 flooredNeutronPos = Core->assembly->returnFlooredNeutronPosInPincell(Bank->addedNeutrons[idx]);
#ifdef TALLYINCYCLE
                if (currentPincell.meatOrMod(flooredAddedNeutronPos) == MatType::MOD) { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]).modTally), DTC); }
                else { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]).pinTally), DTC); }
#endif
                Interaction::reaction(Bank->addedNeutrons[idx], Bank, XSLib, currentPincell, flooredAddedNeutronPos, RNG, k_mult, passFlag, true);
                seedNo[idx] = RNG.gen();
                return;
            }
            else {  // do:  boundary check / reflection / position update to DTS and feed it back to main loop
                //vec3 updatedPos = Bank->addedNeutrons[idx].pos + Bank->addedNeutrons[idx].dirVec * DTC;
#ifdef TALLYINCYCLE
                if (currentPincell.meatOrMod(flooredAddedNeutronPos) == MatType::MOD) { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]).modTally), DTS); }
                else { atomicAdd(&(currentTallyAssembly.returnPincellByPos(Bank->addedNeutrons[idx]).pinTally), DTS); }
#endif
                vec3 updatedSurfacePos = Bank->addedNeutrons[idx].pos + Bank->addedNeutrons[idx].dirVec * DTS;
                vec3 afterDTCPos = Bank->addedNeutrons[idx].pos + Bank->addedNeutrons[idx].dirVec * DTC;
                // handle vaccum boundary neutrons;

                if (updatedSurfacePos.x >= Core->x  || updatedSurfacePos.y >= Core->y || updatedSurfacePos.z >= Core->z) {
                    Bank->addedNeutrons[idx].Nullify();
                    atomicAdd(&(Bank->addedNeutronSize), -1);
                    seedNo[idx] = RNG.gen();
                    return;
                }

                if (afterDTCPos.x < 0.0 || afterDTCPos.y < 0.0 || afterDTCPos.z < 0.0) {
                    if (updatedSurfacePos.x <= eps / 10 || updatedSurfacePos.y <= eps / 10 || updatedSurfacePos.z <= eps / 10 ) {
                        Interaction::reflection(Bank->addedNeutrons[idx], DTS, updatedSurfacePos, eps);
#ifdef REFLECTIONDEBUG
                        if (reflectionCounter > 1) {
                            printf("Neutron reflected twice on ");
                            Bank->addedNeutrons[idx].printInfo_Kernel(idx);
                        }
#endif
                        if (Bank->addedNeutrons[idx].isNullified()) {
                            // error handling - reflection might return a nullified neutron - thus reduce the size of neutron by 1.
                            atomicAdd(&(Bank->addedNeutronSize), -1);
                            seedNo[idx] = RNG.gen();
                            return;
                        }
                        continue;
                    }
                }
                // neutron is inside the boundary.
                Bank->addedNeutrons[idx].updateWithLength(DTS + eps * 1000);
                
                //Bank->addedNeutrons[idx].dirVec.x >= 0 ? Bank->addedNeutrons[idx].pos.x += globaleps : Bank->addedNeutrons[idx].pos.x -= globaleps;
                //Bank->addedNeutrons[idx].dirVec.y >= 0 ? Bank->addedNeutrons[idx].pos.y += globaleps : Bank->addedNeutrons[idx].pos.y -= globaleps;
                //Bank->addedNeutrons[idx].dirVec.z >= 0 ? Bank->addedNeutrons[idx].pos.z += globaleps : Bank->addedNeutrons[idx].pos.z -= globaleps;
                /*
                vec3 epsVec = {
                    (Bank->addedNeutrons[idx].dirVec.x >= 0 ? +globaleps : -globaleps),
                    (Bank->addedNeutrons[idx].dirVec.y >= 0 ? +globaleps : -globaleps),
                    (Bank->addedNeutrons[idx].dirVec.z >= 0 ? +globaleps : -globaleps)
                };
                Bank->addedNeutrons[idx].pos = Bank->addedNeutrons[idx].pos + epsVec;
                */
                
            }
            //printf("idx %d neutron on %d loop\n", idx, i);

            if (i == iterTimes - 1) {
                Bank->addedNeutrons[idx].printInfo_Kernel(idx);
                Bank->addedNeutrons[idx].Nullify();
                atomicAdd(&(Bank->addedNeutronSize), -1);
            }
        }
        //printf("idx %d neutron didn't reacted after %d loops..?\n", idx, iterTimes);
        seedNo[idx] = RNG.gen();
        return;
    }
    seedNo[idx] = RNG.gen();
    return;


}

H void cycle_Neutron_CPU(NeutronBank* Bank, C5G7Geometry* Core, TallyC5G7Geometry* CoreTally, XSLibrary* XSLib, unsigned long long* seedNo, double* k_mult, bool passFlag, double& absorption, double& fission, double& leak) {
    for (int idx = 0; idx < Bank->allocatableNeutronNum; idx++) {
        GnuAMCM RNG(seedNo[idx]);
        double eps = globaleps;

        if (Bank->neutrons[idx].isNullified()) {
            seedNo[idx] = RNG.gen();
            continue;
        }

        for (int it = 0; it < 100; it++) {
            Assembly currentAssembly = Core->returnAssemblyByNeutron(Bank->neutrons[idx]);

            if (Bank->neutrons[idx].isNullified()) {
                Bank->neutronSize -= 1;
                seedNo[idx] = RNG.gen();
                break;
            }

            Pincell currentPincell = currentAssembly.returnPincellByPos(Bank->neutrons[idx]);
            vec3 flooredPos = currentAssembly.returnFlooredNeutronPosInPincell(Bank->neutrons[idx]);

            if (currentPincell.sideLength == 0.0 && currentPincell.height == 0.0) {
                Bank->neutrons[idx].Nullify();
                Bank->neutronSize -= 1;
                seedNo[idx] = RNG.gen();
                break;
            }

            double DTC = currentAssembly.DTC(Bank->neutrons[idx], XSLib, RNG);
            double DTS = currentAssembly.DTS(Bank->neutrons[idx]);

            if (DTC < DTS) {
                Bank->neutrons[idx].updateWithLength(DTC);

                if (CoreTally) {
                    TallyAssembly& TA = CoreTally->returnAssemblyByNeutron(Bank->neutrons[idx]);
                    TallyPincell& TP = TA.returnPincellByPos(Bank->neutrons[idx]);
                    if (currentPincell.meatOrMod(flooredPos) == MatType::MOD) TP.modTally += DTC;
                    else TP.pinTally += DTC;
                }

                Interaction::reaction_CPU(idx, Bank->neutrons[idx], Bank, XSLib, currentPincell, flooredPos, RNG, k_mult, passFlag, false, absorption, fission);
                seedNo[idx] = RNG.gen();
                break;
            }
            else {
                if (CoreTally) {
                    TallyAssembly& TA = CoreTally->returnAssemblyByNeutron(Bank->neutrons[idx]);
                    TallyPincell& TP = TA.returnPincellByPos(Bank->neutrons[idx]);
                    if (currentPincell.meatOrMod(flooredPos) == MatType::MOD) TP.modTally += DTS;
                    else TP.pinTally += DTS;
                }

                vec3 updatedSurfacePos = Bank->neutrons[idx].pos + Bank->neutrons[idx].dirVec * DTS;

                if (updatedSurfacePos.x >= Core->x - eps || updatedSurfacePos.y >= Core->y - eps || updatedSurfacePos.z >= Core->z - eps) {
                    Bank->neutrons[idx].Nullify();
                    Bank->neutronSize -= 1;
                    leak += 1.0;
                    seedNo[idx] = RNG.gen();
                    break;
                }

                if (updatedSurfacePos.x <= eps || updatedSurfacePos.y <= eps || updatedSurfacePos.z <= eps) {
                    Interaction::reflection_CPU(Bank->neutrons[idx], DTS, updatedSurfacePos, eps);
                    if (Bank->neutrons[idx].isNullified()) {
                        Bank->neutronSize -= 1;
                        seedNo[idx] = RNG.gen();
                        break;
                    }
                    continue;
                }

                Bank->neutrons[idx].updateWithLength(DTS + eps);
            }
        }

        seedNo[idx] = RNG.gen();
    }
}

H void cycle_addedNeutron_CPU(NeutronBank* Bank, C5G7Geometry* Core, TallyC5G7Geometry* CoreTally, XSLibrary* XSLib, unsigned long long* seedNo, double* k_mult, bool passFlag, double& absorption, double& fission, double& leak) {
    for (int idx = 0; idx < Bank->allocatableNeutronNum; idx++) {
        GnuAMCM RNG(seedNo[idx]);
        double eps = globaleps;

        if (Bank->addedNeutrons[idx].isNullified()) {
            seedNo[idx] = RNG.gen();
            continue;
        }

        if (Bank->addedNeutrons[idx].passFlag) {
            Bank->addedNeutrons[idx].passFlag = false;
            seedNo[idx] = RNG.gen();
            continue;
        }

        for (int it = 0; it < 100; it++) {
            Assembly currentAssembly = Core->returnAssemblyByNeutron(Bank->addedNeutrons[idx]);

            if (Bank->addedNeutrons[idx].isNullified()) {
                Bank->addedNeutronSize -= 1;
                seedNo[idx] = RNG.gen();
                break;
            }

            Pincell currentPincell = currentAssembly.returnPincellByPos(Bank->addedNeutrons[idx]);
            vec3 flooredPos = currentAssembly.returnFlooredNeutronPosInPincell(Bank->addedNeutrons[idx]);

            if (currentPincell.sideLength == 0.0 && currentPincell.height == 0.0) {
                Bank->addedNeutrons[idx].Nullify();
                Bank->addedNeutronSize -= 1;
                seedNo[idx] = RNG.gen();
                break;
            }

            double DTC = currentAssembly.DTC(Bank->addedNeutrons[idx], XSLib, RNG);
            double DTS = currentAssembly.DTS(Bank->addedNeutrons[idx]);

            if (DTC < DTS) {
                Bank->addedNeutrons[idx].updateWithLength(DTC);

                if (CoreTally) {
                    TallyAssembly& TA = CoreTally->returnAssemblyByNeutron(Bank->addedNeutrons[idx]);
                    TallyPincell& TP = TA.returnPincellByPos(Bank->addedNeutrons[idx]);
                    if (currentPincell.meatOrMod(flooredPos) == MatType::MOD) TP.modTally += DTC;
                    else TP.pinTally += DTC;
                }

                Interaction::reaction_CPU(idx, Bank->addedNeutrons[idx], Bank, XSLib, currentPincell, flooredPos, RNG, k_mult, passFlag, true, absorption, fission);
                seedNo[idx] = RNG.gen();
                break;
            }
            else {
                if (CoreTally) {
                    TallyAssembly& TA = CoreTally->returnAssemblyByNeutron(Bank->addedNeutrons[idx]);
                    TallyPincell& TP = TA.returnPincellByPos(Bank->addedNeutrons[idx]);
                    if (currentPincell.meatOrMod(flooredPos) == MatType::MOD) TP.modTally += DTS;
                    else TP.pinTally += DTS;
                }

                vec3 updatedSurfacePos = Bank->addedNeutrons[idx].pos + Bank->addedNeutrons[idx].dirVec * DTS;

                if (updatedSurfacePos.x >= Core->x - eps || updatedSurfacePos.y >= Core->y - eps || updatedSurfacePos.z >= Core->z - eps) {
                    Bank->addedNeutrons[idx].Nullify();
                    Bank->addedNeutronSize -= 1;
                    leak += 1.0;
                    seedNo[idx] = RNG.gen();
                    break;
                }

                if (updatedSurfacePos.x <= eps || updatedSurfacePos.y <= eps || updatedSurfacePos.z <= eps) {
                    Interaction::reflection_CPU(Bank->addedNeutrons[idx], DTS, updatedSurfacePos, eps);
                    if (Bank->addedNeutrons[idx].isNullified()) {
                        Bank->addedNeutronSize -= 1;
                        seedNo[idx] = RNG.gen();
                        break;
                    }
                    continue;
                }

                Bank->addedNeutrons[idx].updateWithLength(DTS + eps);
            }
        }

        seedNo[idx] = RNG.gen();
    }
}

H void addedNeutronPassResetter_CPU(NeutronBank* Bank) {
    for (int idx = 0; idx < Bank->allocatableNeutronNum; idx++) {
        if (!Bank->addedNeutrons[idx].isNullified()) {
            Bank->addedNeutrons[idx].passFlag = false;
        }
    }
}