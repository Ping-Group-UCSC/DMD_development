#include "DensityMatrix_Dynamics.h"
bool DEBUG = false;
string dir_debug = "debug_info/";
bool ionode = false;
algorithm alg;
string code = "";
string material_model = "";
ODEparameters ode;

void dm_dynamics_jdftx(parameters* param);
void dm_dynamics_mos2(parameters* param);

int main(int argc, char **argv)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now(), t2;

	// mpi
	mpkpair.mpi_init(); // kpair parallel
	mpkpair2.mpi_init(); // kpair parallel
	mpk.mpi_init(); // k parallel
	mpk_morek.mpi_init(); // k parallel
	ionode = mpkpair.ionode;

	// read parameters
	parameters* param = new parameters();
	param->read_param();
	if (ionode){
		printf("HELLO THERE\n");
		printf("TESTA DI PIGNA !!\n");
		printf("SI VERGOGNI !!!\n");
	}
	init_model(param);
	if (material_model == "none") dm_dynamics_jdftx(param);

	MPI_Barrier(MPI_COMM_WORLD);
	t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(t2 - t1).count();
	if (ionode) cout << "total time: " << duration / 1.0e6 << " seconds" << endl;
	return 0;
}

void dm_dynamics_jdftx(parameters* param){
	//lattice
	lattice* latt = new lattice(param);
	latt->printLattice();
	// electron
	electron* elec = new electron(&mpk, &mpk_morek, latt, param);
	mpk.distribute_var("dm_dynamics_driver", elec->nk);
	if (elec->nk_morek > 0) mpk_morek.distribute_var("dm_dynamics_driver", elec->nk_morek);
	param->mu = elec->set_mu_and_n(param->carrier_density); // if carrier_density is set (non-zero), recompute mu
	elec->compute_b2(param->de_measure, param->degauss_measure, param->degthr);
	elec->compute_Bin2();
	elec->compute_dm_Bpert_1st(param->Bpert, param->t0);
	if (elec->B.length() > 1e-10 || alg.read_Bso) elec->set_H_BS(mpk.varstart, mpk.varend);
	if (abs(elec->scale_Ez) > 1e-10) elec->set_H_Ez(mpk.varstart, mpk.varend);
	std::cout << "OK2" << std::endl;
	//if ((alg.picture == "schrodinger" || param->t0 == 0) && param->Bpert.length() > 1e-12) elec->compute_DP_related(param->Bpert); // Not needed for this moment. May add back later
	// phonon
	phonon* ph = new phonon(latt, param, elec); // may need elec->kmesh and elec->kvec to construct qvec
	std::cout << "OK3" << std::endl;
	// electron-light or laser
	electronlight* elight;
	if (pmp.active()){
		elight = new electronlight(latt, param, elec, &mpk);
		if (pmp.laserAlg == "perturb") elight->pump_pert();
	}
	std::cout << "OK4" << std::endl;
	// electron-phonon
	electronphonon* eph = new electronphonon(&mpkpair, latt, param, elec, ph, alg.eph_sepr_eh, !alg.eph_need_elec);
	std::cout << "OK4A" << std::endl;
	if (alg.scatt_enable) mpkpair.distribute_var("dm_dynamics_jdftx", eph->nkpair_glob);
	std::cout << "OK4B" << std::endl;
	if (alg.scatt_enable) eph->set_eph();
	std::cout << "OK4C" << std::endl;
	if (alg.scatt_enable && alg.linearize && param->need_imsig) eph->compute_imsig();
	std::cout << "OK4D" << std::endl;
	if (alg.scatt_enable) eph->analyse_g2(param->de_measure, param->degauss_measure, param->degthr);
	if (alg.scatt_enable) eph->analyse_g2_ei(param->de_measure, param->degauss_measure, param->degthr);
	std::cout << "OK5" << std::endl;
	
	MPI_Barrier(MPI_COMM_WORLD);
	dm_dynamics<lattice, electron, electronlight, electronphonon>* dmdyn =
		new dm_dynamics<lattice, electron, electronlight, electronphonon>(latt, param, elec, elight, eph);
	
	//==================================================
	// evolve density matrix
	//==================================================
	if (param->compute_tau_only) return;
	if (alg.ode_method == "rkf45")
		dmdyn->evolve_gsl();
	else if (alg.ode_method == "euler")
		dmdyn->evolve_euler();
}