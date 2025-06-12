#pragma once
#include "common_headers.h"

class ph_amplitudes {
public:
    lattice *latt;





    ph_amplitudes (lattice *latt, parameters *param): latt(latt)
    {
        if (ionode) std::cout << "PH. AMPLITUDES INITIALIZATION" << std::endl;
    }














};