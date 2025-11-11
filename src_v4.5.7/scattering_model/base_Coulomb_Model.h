#pragma once
#include "lattice.h"
#include "electron.h"

// functions declaration

void find_wqmax (int nk, int nb, vector<vector3<double>> &kvec, vector<vector3<double>> &qvec, qIndexMap *qmap, double **e, vector<double> &wqmax);
double compute_carrier_density_correction(int nk, int nk_full, int bStart, int bEnd, int nv, lattice *latt, electron *elec);
void calc_ovlp(int ik, int jk, int nb, complex *Uih, electron *elec, complex *ovlp);

// Base template class

class coulomb_model_base {
public:
    mymp *mp;
    lattice *latt;
    electron *elec;
    int nv;                  // number of valence states
    int nk;                  // number k points
    int nb;                  // number bands
    int bStart, bEnd;
    double prefac_vq_bare;
    double prefac_vq;
    double qmin, qmax;
    double nk_full;
    double nfreetot_corr;
    double T;
    vector<double> wqmax;
    vector<complex> omega;
    int iq_qmin;
    vector<vector3<double>> qvec;
    qIndexMap *qmap = nullptr;
    kIndexMap *kmap = nullptr;
    double **e;              // local electron energies
    double **f;              // local electron occupations
    complex *Uih, *ovlp;
    vector<complex> qscr2_static_RPA;
    vector<vector<complex>> vq_RPA;

    void print_qscr2_static_RPA(const std::string &filename = "qscr2_static_RPA.out") const;
    void calc_ovlp(int ik, int jk);
    void init(double **ft);
    void init(complex **dm);

    coulomb_model_base (lattice *latt, parameters *param, electron *elec, int bStart, int bEnd, double dE)
    : latt(latt), nv(elec->nv_dm - bStart), qmin(0), qmax(0), nk(elec->nk), nb(bEnd - bStart), T(param->temperature),
    nk_full(elec->nk_full), omega(clp.nomega), elec(elec), bStart(bStart), bEnd(bEnd)
    {
        if (ionode) printf("\nInitialize screening formula %s\n", clp.scrFormula.c_str());
		if (ionode) printf("bStart = %d bEnd = %d nv = %d\n", bStart, bEnd, nv);
        if (latt->dim < 3 && clp.scrFormula != "RPA") error_message("ONLY RPA screening implemented for 2D systems");
        prefac_vq = 4 * M_PI / clp.eps / latt->cell_size;
        prefac_vq_bare = 4 * M_PI / latt->cell_size;
        e = trunc_alloccopy_array(elec->e_dm, nk, bStart, bEnd);
		f = trunc_alloccopy_array(elec->f_dm, nk, bStart, bEnd);
        if (clp.scrFormula == "RPA" || clp.scrFormula == "lindhard"){
            Uih = new complex[nb*elec->nb_wannier]{c0};
            ovlp = new complex[nb*nb]{c0};
        }

        //nfreetot

        if (ionode) printf("nfreetot = %f\n", clp.nfreetot);
        nfreetot_corr = compute_carrier_density_correction(nk, nk_full, bStart, bEnd, nv, latt, elec);
        if (ionode) std::cout << "nfreetot corr.: " << nfreetot_corr << std::endl;

        //initialize qmap and qvec,

        if (qmap == nullptr) qmap = new qIndexMap(elec->kmesh);
        qmap->build(elec->kvec, qvec);
        if (ionode) { 
            //string fnameq = dir_debug + "qIndexMap.out";
            string fnameq = "qIndexMap.out";
            qmap->print_map(qvec, fnameq);
        }
		if (ionode) std::cout << "INITIALIZED qmap" << std::endl;
        if (ionode) std::cout << "Q VECTOR SIZE: " << qvec.size() << std::endl;

        //initialize kmap
        
        if (kmap == nullptr) kmap = new kIndexMap(elec->kmesh, elec->kvec);
        if (ionode) {
            string fnamek = "kIndexMap.out";
            kmap->print_map(elec->kvec, fnamek);
        }
        if (ionode) std::cout << "INITIALIZED kmap" << std::endl;
        
        //q vector min/max

        if (qmax == 0 & qmin == 0) latt->find_qmin_qmax(qvec, qmin, qmax, iq_qmin);
        if (ionode) printf("qmin = %lg qmax = %lg\n", qmin, qmax);

        if (clp.dynamic != "static") {
            find_wqmax (nk, nb, elec->kvec, qvec, qmap, e, wqmax);
            if (ionode) {
                string fname = "wmax_q_kbt.out";
                FILE *fp = fopen(fname.c_str(), "w");
                fprintf(fp, "#|q|(au) wmax(kbT)\n");
                for (size_t iq=0; iq<qvec.size(); iq++) {
                    double q_length = sqrt(latt->GGT.metric_length_squared(wrap(qvec[iq])));
                    fprintf(fp, "%10.3le %10.3le\n", q_length, wqmax[iq] / T);
                }
                fclose(fp);
            }
            // for dynamic screening
            if (clp.omegamax == 0) clp.omegamax = dE;
            if (clp.dynamic == "ppa" && clp.ppamodel == "gn") { omega.resize(2); omega[0] = c0; }
            if (!(clp.dynamic == "ppa" && clp.ppamodel == "gn")) clp.smearing = (clp.smearing <= 0) ? 0.5 * T : clp.smearing;
            if (ionode) printf("smearing = %10.3le a.u. (%10.3le meV / %10.3le K)\n", clp.smearing, clp.smearing / eV * 1000, clp.smearing / Kelvin);

        }

        // mpi

        mp = elec->mp;

    }

    virtual void init_model (double n, FILE *fp = stdout) = 0;
    virtual void print () = 0;
    virtual void init_RPA (double **ft = nullptr) = 0;
    virtual complex vq(vector3<double> q, double w = 0) = 0;

};