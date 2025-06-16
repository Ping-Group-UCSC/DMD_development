#pragma once
#include "common_headers.h"
#include "phonon.h"

class ph_amplitudes {
public:
    int modeStart, modeEnd;
    std::vector<complex> Bq;




    ph_amplitudes (parameters *param, phonon *ph): modeStart(param->modeStart), modeEnd(param->modeEnd)
    {
        // OK
        init (ph);
    }


    void init (phonon *ph) {

        if (ionode) std::cout << "PH. AMPLITUDES INITIALIZATION" << std::endl;
        if (ionode) std::cout << "N. MODES: " << ph->nm << std::endl;
        if (ionode) std::cout << "MODES ST. -> END: " << modeStart << " - " << modeEnd << std::endl;

        alloc_Bq_arrays ();
    }

    void alloc_Bq_arrays ();









};