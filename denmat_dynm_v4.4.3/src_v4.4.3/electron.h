#pragma once
#include "common_headers.h"
#include "lattice.h"
#include "parameters.h"
#include "PumpProbe.h"
#include "ScattModel_Param.h"

class electron{
public:
	mymp *mp;
	lattice *latt;
	double temperature, mu, carrier_density, ne, nh;
	vector3<int> kmesh; double nk_full;
	int ns, nk, nb, nb_wannier, bskipped_wannier; // spin, k point, band
	int nv, nc; // valence and conduction bands
	int bStart_dm, bEnd_dm, nb_dm, nv_dm, nc_dm; // band range related to density matrix
	int bStart_eph, bEnd_eph, nb_eph, nv_eph, nc_eph;
	std::vector<vector3<>> kvec;
	bool print_along_kpath;
	std::vector<vector3<double>> kpath_start, kpath_end;
	int nkpath;
	std::vector<int> ik_kpath;
	double **e, **f, **e_dm, **f_dm; // ocuupation number
	double emin, emax, evmax, ecmin, emid;
	vector3<> B;
	complex ***s, ***v, **U;
	complex **H_BS;
	double degthr;
	complex **ddm_Bpert, **dm_Bpert, **ddm_Bpert_neq, **dm_Bpert_neq; // store ddm_eq and ddm_neq for later analysis of ddm
	bool *has_precess;
	double trace_sq_ddm_tot;
	double *DP_precess_fac;
	double **imsig_eph_kn, *imsig_eph_k, imsig_eph_avg;

	electron(parameters *param)
		:temperature(param->temperature), mu(param->mu), carrier_density(param->carrier_density),
		kmesh(vector3<int>(param->nk1, param->nk2, param->nk3)), nk_full((double)param->nk1*(double)param->nk2*(double)param->nk3), B(param->B),
		print_along_kpath(param->print_along_kpath), kpath_start(param->kpath_start), kpath_end(param->kpath_end), nkpath(param->kpath_start.size()){}
	electron(mymp *mp, lattice *latt, parameters *param)
		:mp(mp), latt(latt), temperature(param->temperature), mu(param->mu), carrier_density(param->carrier_density),
		kmesh(vector3<int>(param->nk1, param->nk2, param->nk3)), nk_full((double)param->nk1*(double)param->nk2*(double)param->nk3), B(param->B),
		print_along_kpath(param->print_along_kpath), kpath_start(param->kpath_start), kpath_end(param->kpath_end), nkpath(param->kpath_start.size()),
		H_BS(nullptr), ddm_Bpert(nullptr), ddm_Bpert_neq(nullptr), dm_Bpert_neq(nullptr)
	{
		if (code == "jdftx"){
			read_ldbd_size();
			alloc_mat(pmp.pumpA0 > 0, eip.ni != 0 || eep.eeMode != "none");
			read_ldbd_kvec();
			if (print_along_kpath && kvec.size() > 0) get_kpath();
			read_ldbd_ek();
			e_dm = trunc_alloccopy_array(e, nk, bStart_dm, bEnd_dm);
			f_dm = trunc_alloccopy_array(f, nk, bStart_dm, bEnd_dm);
			read_ldbd_imsig_eph();
			read_ldbd_smat();
			if (pmp.pumpA0 > 0) read_ldbd_vmat();
			if (eip.ni != 0 || eep.eeMode != "none") read_ldbd_Umat();
		}
		emax = maxval(e_dm, nk, nv_dm, nb_dm);
		ecmin = minval(e_dm, nk, nv_dm, nb_dm);
		if (ionode) printf("\nElectronic states energy range:\n");
		if (nb_dm > nv_dm && ionode) printf("\nemax = %lg ecmin = %lg\n", emax, ecmin);
		if (nb_dm > nv_dm && eip.partial_ionized && eip.Eimp > ecmin) error_message("Impurity level should not be higher than CBM", "electron constructor");
		evmax = maxval(e_dm, nk, 0, nv_dm);
		emin = minval(e_dm, nk, 0, nv_dm);
		emid = (ecmin + evmax) / 2;
		if (nv_dm > 0 && ionode) printf("evmax = %lg emin = %lg\n", evmax, emin);
		if (nv_dm > 0 && eip.partial_ionized && eip.Eimp < evmax) error_message("Impurity level should not be lower than VBM", "electron constructor");
		eip.ni_cell = eip.ni * nk_full * latt->volume;
		if (nv_dm > 0 && nb_dm > nv_dm && eip.partial_ionized){
			if (eip.ni > 0 && eip.Eimp < emid) error_message("impurity level < middle energy in gap in n-type", "electron constructor");
			if (eip.ni < 0 && eip.Eimp > emid) error_message("impurity level > middle energy in gap in p-type", "electron constructor");
		}
	}
	void alloc_mat(bool alloc_v, bool alloc_U);

	void get_kpath();
	vector3<> get_kvec(int&, int&, int&);
	static double fermi(double t, double mu, double e){
		double ebyt = (e - mu) / t;
		if (ebyt < -46) return 1;
		else if (ebyt > 46) return 0;
		else return 1. / (exp(ebyt) + 1);
	}
	void compute_f(double t, double mu){
		zeros(f, nk, nb);
		for (int ik = mp->varstart; ik < mp->varend; ik++)
		for (int i = 0; i < nb; i++)
			f[ik][i] = fermi(t, mu, e[ik][i]);
		mp->allreduce(f, nk, nb, MPI_SUM);
		trunc_copy_array(f_dm, f, nk, bStart_dm, bEnd_dm);
	}
	double find_mu(double carrier_cell, double t, double mu0);
	double compute_nfree(bool isHole, double t, double mu);
	double set_mu_and_n(double carrier_density){
		if (carrier_density != 0){
			mu = find_mu(carrier_density * nk_full * latt->volume, temperature, mu);
			compute_f(temperature, mu);
			if (ionode) print_array_atk(f, nb, "f:");
			if (ionode) print_array_atk(f_dm, nb_dm, "f_dm:");
		}
		ne = compute_nfree(false, temperature, mu) / nk_full / latt->volume;
		nh = compute_nfree(true, temperature, mu) / nk_full / latt->volume;
		double nei = eip.compute_carrier_cell_of_impurity_level(false, temperature, mu) / nk_full / latt->volume;
		double nhi = eip.compute_carrier_cell_of_impurity_level(true, temperature, mu) / nk_full / latt->volume;
		if (ionode) printf("ne = %lg nei = %lg ne+nei = %lg\n", ne, nei, ne + nei);
		if (ionode) printf("nh = %lg nhi = %lg nh+nhi = %lg\n", nh, nhi, nh + nhi);
		eip.carrier_cell_gs = carrier_density * nk_full * latt->volume;
		eip.ne_cell_gs = ne * nk_full * latt->volume; eip.nh_cell_gs = nh * nk_full * latt->volume; clp.nfreetot = fabs(ne) + fabs(nh);
		eip.calc_ni_ionized(temperature, mu);
		return mu;
	}

	void read_ldbd_size();
	void read_ldbd_kvec();
	void read_ldbd_ek();
	void read_ldbd_imsig_eph();
	void read_ldbd_smat();
	void read_ldbd_vmat();
	void read_ldbd_Umat();
	void print_array_atk(bool *a, string s = "");
	void print_array_atk(double *a, string s = "", double unit = 1);
	void print_array_atk(double **a, int n, string s = "", double unit = 1);
	void print_mat_atk(complex **m, int n, string s = "");

	void compute_dm_Bpert_1st(vector3<> Bpert, double t0);
	void compute_DP_related(vector3<> Bpert);

	void set_H_BS(int ik0_glob, int ik1_glob);
};