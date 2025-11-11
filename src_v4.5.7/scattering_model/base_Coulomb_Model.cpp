#include "base_Coulomb_Model.h"

// compute wqmax array

void find_wqmax (int nk, int nb, vector<vector3<double>> &kvec, vector<vector3<double>> &qvec, qIndexMap *qmap, double **e, vector<double> &wqmax) {
    wqmax.resize(qvec.size(), 0);
    for (int ik1=0; ik1<nk; ik1++) {
        for (int ik2=0; ik2<nk; ik2++) {
            size_t iq = qmap->q2iq(kvec[ik1] - kvec[ik2]);
            for (int b1=0; b1<nb; b1++){
                for (int b2=0; b2<nb; b2++){
                    double w = fabs(e[ik1][b1] - e[ik2][b2]);
                    if (w > wqmax[iq]) wqmax[iq] = w;
                }
            }
        }
    }
}

// carrier density correction

double compute_carrier_density_correction(int nk, int nk_full, int bStart, int bEnd, int nv, lattice *latt, electron *elec) {
    double nfreetot_corr = 0;
    if (elec->nk_morek) {
        for (int ik = 0; ik < nk; ik++) {
            for (int i=bStart; i<bEnd; i++) {
                if (i >= nv) nfreetot_corr -= elec->f_dm[ik][i];
                else nfreetot_corr -= (elec->f_dm[ik][i] - 1.); // hole concentration
            }
        }
        for (int ik = 0; ik < elec->nk_morek; ik++) {
            for (int i = bStart; i < bEnd; i++) {
                if (i>=nv) nfreetot_corr += elec->f_dm_morek[ik][i];
                else nfreetot_corr += (elec->f_dm_morek[ik][i] - 1.);
            }
        }
    }
    nfreetot_corr /= (nk_full * latt->cell_size);
    return nfreetot_corr;
}

// overlap matrix calculation < Psi_k | Psi_k' >

void coulomb_model_base::calc_ovlp(int ik, int jk) {
    hermite(elec->U[ik], Uih, elec->nb_wannier, nb);
    zgemm_interface(ovlp, Uih, elec->U[jk], nb, nb, elec->nb_wannier);
}

// print qscr2_static_RPA

void coulomb_model_base::print_qscr2_static_RPA(const std::string &filename) const {
    if (!ionode) return;  // avoid writing from all MPI ranks
    FILE *fpvq = fopen(filename.c_str(), "w");
    fprintf(fpvq, "#|q|^2 |q_scr|^2\n");
    for (size_t iq = 0; iq < qvec.size(); iq++) {
        double q_len_sq = latt->GGT.metric_length_squared(wrap(qvec[iq]));
        fprintf(fpvq, "%14.7le %14.7le\n", q_len_sq, abs(qscr2_static_RPA[iq]));
    }
    fclose(fpvq);
}

// base init functions

void coulomb_model_base::init(double **ft){
    clp.nfreetot = 0;
    for (int ik = 0; ik < nk; ik++)
	for (int b = 0; b < nb; b++){
		if (b < nv) clp.nfreetot += (1 - ft[ik][b]);
		else clp.nfreetot += ft[ik][b];
	}
	clp.nfreetot /= (nk_full * latt->cell_size);
	clp.nfreetot += nfreetot_corr;
	if (ionode) printf("nfree = %lg cm-3 for screening\n", clp.nfreetot / std::pow(bohr2cm, 3));
	init_model(clp.nfreetot);
	if (clp.scrFormula == "RPA") init_RPA();
}

void coulomb_model_base::init(complex **dm){
	double **ft = alloc_real_array(nk, nb);
	for (int ik = 0; ik < nk; ik++)
	for (int b = 0; b < nb; b++)
		ft[ik][b] = real(dm[ik][b*nb + b]);
	init(ft);
	dealloc_real_array(ft);
}