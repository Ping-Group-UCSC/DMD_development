#pragma once
#include "common_headers.h"
#include "parameters.h"

class phenom_recomb{

public:
	double tau;
	int nk_glob, nv_dm, nb_dm;
    double **f_dm;
	complex **dm, **ddmdt_phenom;

	phenom_recomb(parameters *param, int nk_glob, int nb_dm, int nv_dm)
		: tau(param->tau_phenom_recomb), nk_glob(nk_glob), nb_dm(nb_dm), nv_dm(nv_dm)
	{
		if (ionode) printf("\n");
		if (ionode) printf("==================================================\n");
		if (ionode) printf("==================================================\n");
		if (ionode) { printf("Initialize phenom_recomb with tau_phenom_recomb = %d\n", tau); fflush(stdout); }
		if (ionode) printf("==================================================\n");
		if (ionode) printf("==================================================\n");
		f_dm = alloc_real_array(nk_glob, nb_dm);
		dm = alloc_array(nk_glob, nb_dm*nb_dm);
		ddmdt_phenom = alloc_array(nk_glob, nb_dm*nb_dm);
	}

    void evolve(complex **dm, double **f_dm, complex **ddmdt_phenom);
};

void phenom_recomb::evolve(complex **dm, double **f_dm, complex **ddmdt_phenom)
{
    zeros(ddmdt_phenom, nk_glob, nb_dm*nb_dm);
    for (int ik_glob = 0; ik_glob < nk_glob; ik_glob++)
    {
        for (int i = 0; i < nb_dm; i++)
        {
            for (int j = 0; j < nb_dm; j++)
            {
                bool within_val_or_cond = (i < nv_dm and j < nv_dm) or (i >= nv_dm and j >= nv_dm);
                if (within_val_or_cond and i != j)
                    continue;
                
                if (i != j)
                    ddmdt_phenom[ik_glob][i*nb_dm + j] += -dm[ik_glob][i*nb_dm + j] / tau;
                else
                    ddmdt_phenom[ik_glob][i*nb_dm + j] += -(dm[ik_glob][i*nb_dm + j] - f_dm[ik_glob][i]) / tau;
            }
        }
    }
}