#pragma once
#include "base_Coulomb_Model.h"

// template value 2

class coulomb_model_2D : public coulomb_model_base {
public:

    vector<complex> cutoff_prefac;

    coulomb_model_2D (lattice *latt, parameters *param, electron *elec, int bStart, int bEnd, double dE)
    : coulomb_model_base (latt, param, elec, bStart, bEnd, dE)
    {
        
        init_model(clp.nfreetot);

        // RPA initialization

        init_RPA ();

    }

    // initialize Coulomb model

    void init_model(double n, FILE *fp = stdout) override {
        if (ionode) std::cout << "2D MODEL INITIALIZER" << std::endl;
        if (ionode) std::cout << "screening formula: " << clp.scrFormula << std::endl;
        if (clp.scrFormula != "RPA") error_message ("ONLY RPA screening implemented in the 2D case");
    }

    void init_RPA (double **ft = nullptr) {
        if (ionode) std::cout << "RPA calculation initialization" << std::endl;
        if (ionode) std::cout << "dynamic calculation: " << clp.dynamic << std::endl;
        if (clp.dynamic != "static") error_message ("ONLY static RPA screening implemented in the 2D case");
        if (ft != nullptr) trunc_copy_array(f, ft, nk, 0, nb);

        // 2D cutoff

        calc_2d_cutoff();

        // qscr2 static RPA

        calc_qscr2_static_RPA ();

        calc_vq_RPA ();
        
    }

    void calc_2d_cutoff () {
        if (ionode) std::cout << "-> start 2D COULOMB CUTOFF CALCULATION" << std::endl;
        cutoff_prefac.resize(qvec.size(), c0);
        // define cutoff distance (a.u.)
        double lz = 0.5 * latt->R(2,2);
        if (ionode) std::cout << "lz -> " << lz << std::endl;
        for (int iq=0; iq < qvec.size(); iq++) {
            vector3<double> q_z {0.0, 0.0, qvec[iq][2]};
            vector3<double> q_plane {qvec[iq][0], qvec[iq][1], 0.0};
            double Qz = sqrt(latt->GGT.metric_length_squared(wrap(q_z)));
            if (Qz > 0) Qz = Qz * q_z[2] / fabs(q_z[2]);
            double Qp = sqrt(latt->GGT.metric_length_squared(wrap(q_plane)));
            double Qplz = Qp * lz;
            double Qzlz = Qz * lz;
            cutoff_prefac[iq] = 1.0 - exp(-Qplz) * cos(Qzlz);
        }
    }

    void calc_qscr2_static_RPA () {
        std::cout << mp->myrank << " -> kvec size: " << elec->kvec.size() << std::endl;
        std::cout << "st " << mp->nprocs << " - " << mp->myrank << " - " << mp->varstart << " -> " << mp->varend << std::endl;
        qscr2_static_RPA.resize(qvec.size(), c0);
        for (int iq = 0; iq < qvec.size(); iq++) {
            if (ionode) std::cout << "iq --------> " << iq << std::endl;
            for (int ik = mp->varstart; ik < mp->varend; ik++){
                size_t jk = 0;
                vector3<> kvj = elec->kvec[ik] - qvec[iq];
                if (kmap->findk(kvj, jk)){
                    calc_ovlp(ik, jk);
                    for (int b1 = 0; b1 < nb; b1++) {
                        for (int b2 = 0; b2 < nb; b2++) {
                            complex de = e[ik][b1] - e[jk][b2], dfde = c0;
                            if (abs(de) < 1e-8) {
                                double favg = 0.5 * (f[ik][b1] + f[jk][b2]);
                                dfde = complex((1 - favg) * favg / T, 0);
                            }
                            else dfde = complex(f[jk][b2] - f[ik][b1], 0) / de;
                            qscr2_static_RPA[iq] += dfde * ovlp[b1*nb + b2].norm();
                        }
                    }
                }
            }
            mp->allreduce(qscr2_static_RPA[iq], MPI_SUM);
            qscr2_static_RPA[iq] = complex(prefac_vq / nk_full, 0) * qscr2_static_RPA[iq];
        }
        if (ionode) print_qscr2_static_RPA();
        if (ionode) printf("\ncalc_qscr2_static_RPA done\n");
        return;
    }

    void calc_vq_RPA () {
        if (ionode) std::cout << " -> start 2D COULOMB VQ CALCULATION" << std::endl;
        if (clp.dynamic == "static") {
            vq_RPA.resize(qvec.size());
            for (size_t iq=0; iq < qvec.size(); iq++) {
                vq_RPA[iq].resize(1, c0);
                double q_length_square = latt->GGT.metric_length_squared(wrap(qvec[iq]));
                vq_RPA[iq][0] = complex(prefac_vq, 0) / (q_length_square + qscr2_static_RPA[iq]) * cutoff_prefac[iq];
            }
        }
        else {
            std::cout << "2D COULOMB INTERACTION -> ONLY STATIC MODE IMPLEMENTED" << std::endl;
            exit(1);
        }
    }

    complex vq(vector3<double> q, double w = 0) override {
		double q_length_square = latt->GGT.metric_length_squared(wrap(q));
        if (q_length_square < 1e-20) return c0;  // skip Gamma
		if (clp.scrFormula == "unscreened") {
			return complex(prefac_vq / q_length_square, 0);
		}
		else if (clp.scrFormula == "RPA") {
			size_t iq = qmap->q2iq(q);
			if (clp.dynamic == "static") return vq_RPA[iq][0];
			else std::cout << "2D ONLY STATIC SCREENING IMPLEMENTED" << std::endl;
		}
		else {
			std::cout << "ONLY RPA 2D SCREENING" << std::endl;
		}
    }

    void print() override {
        std::cout << "2D COULOMB CLASS" << std::endl;
    }

};


// template value 3

class coulomb_model_3D : public coulomb_model_base {
public:

    bool ldebug;
    double kF, vF, kF2, EF, qTF2, qscr2_TF;
    double fac0_Bechstedt, fac2_Bechstedt, fac4_Bechstedt;
    double qscr2_debye;
    homogeneous_electron_gas *heg;
	// frequency variables
	double domega;
	vector<vector<complex>> omegaq;  // frequencies for each q - only set max
	vector<complex> Aq_ppa, Eq2_ppa;
	double wp2;

    coulomb_model_3D (lattice *latt, parameters *param, electron *elec, int bStart, int bEnd, double dE)
    : coulomb_model_base (latt, param, elec, bStart, bEnd, dE)
    {
        init_model(clp.nfreetot);

        // compute electronic screening

        if (clp.scrFormula == "RPA") init_RPA ();
    }

    // initialize Coulomb model

    void init_model(double n, FILE *fp = stdout) override {
        if (ionode) std::cout << "3D MODEL INITIALIZER" << std::endl;
        if (ionode) std::cout << "screening formula: " << clp.scrFormula << std::endl;
        if ((clp.scrFormula == "debye" || clp.scrFormula == "Bechstedt" || clp.scrFormula == "heg") && n <= 0)
			error_message("nfreetot must be postive");
		kF = std::pow(3 * M_PI*M_PI * n, 1. / 3.);
		vF = kF / clp.meff;
		kF2 = kF*kF;
		EF = kF2 / 2. / clp.meff;
		if (ionode) printf("kF = %lg EF = %lg\n", kF, EF);
		qTF2 = 6 * M_PI * n / clp.eps / EF;
		if (clp.scrFormula == "Bechstedt"){
			fac0_Bechstedt = 1. / (clp.eps - 1);
			fac2_Bechstedt = 1. / qTF2;
			fac4_Bechstedt = 3. / 4. / kF2 / qTF2;
		}
		qscr2_debye = 4 * M_PI * n / clp.eps / T; // Eq. 11 in PRB 94, 085204 (2016)
		qscr2_TF = 6 * M_PI * n / clp.eps / EF;
		if (ionode) printf("qscr2_debye = %lg qscr2_TF = %lg\n", qscr2_debye, qscr2_TF);

		double wp = sqrt(4 * M_PI * n / clp.meff / clp.eps);
		if (ionode) printf("wp = %10.3le a.u. (%10.3le meV / %10.3le K)\n", wp, wp / eV * 1000, wp / Kelvin);

		if (clp.scrFormula == "heg" || (ldebug && clp.dynamic == "real-axis")) {
			if (heg != nullptr){
				bool update_heg = fabs(n - heg->n) / std::pow(bohr2cm, latt->dim) > 1; //if n is changed, heg needs to be updated as other parameters are all fixed
				if (update_heg){
					delete heg;
					heg = new homogeneous_electron_gas(n, T, clp.meff, clp.eps, kF, vF, EF, qvec, iq_qmin, qmin, qmax, qmap, latt, wqmax, wp);
				}
			}
			else
				heg = new homogeneous_electron_gas(n, T, clp.meff, clp.eps, kF, vF, EF, qvec, iq_qmin, qmin, qmax, qmap, latt, wqmax, wp);
		}
    }

    void init_RPA (double **ft = nullptr) {
        if (ionode) std::cout << "RPA calculation initialization" << std::endl;
        if (clp.eppa == 0) clp.eppa = sqrt(4 * M_PI * clp.nfreetot / clp.meff / clp.eps);
		if (clp.dynamic == "ppa") omega[1] = ci * clp.eppa;
		if (ft != nullptr) trunc_copy_array(f, ft, nk, 0, nb);
		
		calc_qscr2_static_RPA();

        //determine frequency grids (for each q)		
        if (clp.dynamic == "real-axis"){
			domega = clp.omegamax / (clp.nomega - 1);
			omega[0] = c0;
			for (int io = 1; io < clp.nomega; io++)
				omega[io] = omega[io - 1] + domega;

			double prefac_ratio = (1 + qmin*qmin / qscr2_static_RPA[iq_qmin].abs()) / wqmax[iq_qmin];

			omegaq.resize(qvec.size());
			for (int iq = 0; iq < qvec.size(); iq++){
				//construct wq
				double q2 = latt->GGT.metric_length_squared(wrap(qvec[iq]));
				double qscr2 = qscr2_static_RPA[iq].abs();
				double ratio = q2 < 1e-20 ? 1 : prefac_ratio * wqmax[iq] / (1 + q2 / qscr2);
				int nw = q2 < 1e-20 ? 2 : nw = (int)round(ratio * clp.nomega) + 1; //will not deal with q=0 in this version
				if (nw < 6 && q2 > 1e-20 && fabs(qscr2 / q2) > 0.1) nw = 6;
				if (nw < 2) nw = 2;
				double dw = wqmax[iq] / (nw - 1);
				omegaq[iq].resize(nw);
				omegaq[iq][0] = c0; omegaq[iq][nw - 1] = complex(wqmax[iq], clp.smearing);
				for (int iw = 1; iw < nw - 1; iw++)
					omegaq[iq][iw] = complex(iw*dw, clp.smearing);
			}

			string fname = "wq.out";
			FILE *fp = fopen(fname.c_str(), "w");
			fprintf(fp, "#|q|^2  nw  dw (kBT)\n");
			for (size_t iq = 0; iq < qvec.size(); iq++){
				double q_length_square = latt->GGT.metric_length_squared(wrap(qvec[iq])); // qvec is already wrapped to [-0.5,0.5)
				fprintf(fp, "%10.3le %d %10.3le %10.3le\n", q_length_square, omegaq[iq].size(), omegaq[iq][1].real() / T, omegaq[iq][1].imag() / T);
			}
			fclose(fp);
		}

        calc_vq_RPA();
        
    }

	void calc_qscr2_static_RPA(){
		qscr2_static_RPA.resize(qvec.size(), c0);
		// vq = vq0 / (1 - vq0 * sum_k [(f_k - f_k-q) / (e_k - e_k-q - w - i0)] / nk_full)
		// vq0 = e^2 / V / (eps_r * eps_0) / q^2
		// Therefore, vq = e^2 / V / (eps_r * eps_0) / (q^2 + betas^2)
		// betas^2 = - e^2 / V / (eps_r * eps_0) * sum_k [(f_k - f_k-q) / (e_k - e_k-q - w - i0)] / nk_full
		/*
		vector<complex> qscr2_ref(qvec.size());
		for (int ik = 0; ik < nk; ik++)
		for (int jk = 0; jk < nk; jk++){
			vector3<double> q = elec->kvec[ik] - elec->kvec[jk];
			int iq = qmap->q2iq(q);
			calc_ovlp(ik, jk);
			for (int b1 = 0; b1 < nb; b1++)
			for (int b2 = 0; b2 < nb; b2++){
				if (clp.scrFormula == "lindhard" && b1 != b2) continue;
				complex de = e[ik][b1] - e[jk][b2], dfde = c0;
				if (clp.fderavitive_technique){
					if (abs(de) < 1e-8){
						double favg = 0.5 * (f[ik][b1] + f[jk][b2]);
						dfde = complex((1 - favg) * favg / T, 0); // only true for Fermi-Dirac
					}
					else dfde = complex(f[jk][b2] - f[ik][b1], 0) / de;
				}
				else
					dfde = complex(f[jk][b2] - f[ik][b1], 0) / (de - complex(0, clp.smearing));
				if (clp.scrFormula == "RPA") qscr2_ref[iq] += dfde * ovlp[b1*nb + b2].norm();
				else if (clp.scrFormula == "lindhard") qscr2_ref[iq] += dfde;
			}
		}
		axbyc(qscr2_ref.data(), nullptr, qvec.size(), 0, complex(prefac_vq / nk_full, 0), c0); // y = ax + by + c
		*/
		for (int iq = 0; iq < qvec.size(); iq++){
			if (ionode) std::cout << "iq --------> " << iq << std::endl;
			for (int ik = mp->varstart; ik < mp->varend; ik++){
				size_t jk = 0; vector3<> kj = elec->kvec[ik] - qvec[iq]; // not necessage to wrap k point around Gamma, kmap subroutines will wrap inside
				if (kmap->findk(kj, jk)){
					calc_ovlp(ik, jk);
					for (int b1 = 0; b1 < nb; b1++)
					for (int b2 = 0; b2 < nb; b2++){
						if (clp.scrFormula == "lindhard" && b1 != b2) continue;
						complex de = e[ik][b1] - e[jk][b2], dfde = c0;
						if (abs(de) < 1e-8){
							double favg = 0.5 * (f[ik][b1] + f[jk][b2]);
							dfde = complex((1 - favg) * favg / T, 0); // only true for Fermi-Dirac
						}
						else dfde = complex(f[jk][b2] - f[ik][b1], 0) / de;
						if (clp.scrFormula == "RPA") qscr2_static_RPA[iq] += dfde * ovlp[b1*nb + b2].norm();
						else if (clp.scrFormula == "lindhard") qscr2_static_RPA[iq] += dfde;
					}
				}
			}
			mp->allreduce(qscr2_static_RPA[iq], MPI_SUM);
			qscr2_static_RPA[iq] = complex(prefac_vq / nk_full, 0) * qscr2_static_RPA[iq];
		}

		//2nd implementation of static screening for comparison
		vector<complex> qscr2_2ndway(qvec.size());
		for (int iq = 0; iq < qvec.size(); iq++){
			for (int ik = mp->varstart; ik < mp->varend; ik++){
				size_t jk = 0; vector3<> kj = elec->kvec[ik] - qvec[iq]; // not necessage to wrap k point around Gamma, kmap subroutines will wrap inside
				if (kmap->findk(kj, jk)){
					calc_ovlp(ik, jk);
					for (int b1 = 0; b1 < nb; b1++)
					for (int b2 = 0; b2 < nb; b2++){
						if (clp.scrFormula == "lindhard" && b1 != b2) continue;
						complex de = e[ik][b1] - e[jk][b2] - complex(0, clp.smearing);
						complex dfde = complex(f[jk][b2] - f[ik][b1], 0) / de;
						if (clp.scrFormula == "RPA") qscr2_2ndway[iq] += dfde * ovlp[b1*nb + b2].norm();
						else if (clp.scrFormula == "lindhard") qscr2_2ndway[iq] += dfde;
					}
				}
			}
			mp->allreduce(qscr2_2ndway[iq], MPI_SUM);
			qscr2_2ndway[iq] = complex(prefac_vq / nk_full, 0) * qscr2_2ndway[iq];
		}
		if (ionode){
			string fnamevq = "qscr2_static_RPA.out";
			FILE *fpvq = fopen(fnamevq.c_str(), "w");
			init_model(clp.nfreetot, fpvq);
			fprintf(fpvq, "#|q|^2 |q_scr|^2 |2nd q_scr|^2\n");
			for (size_t iq = 0; iq < qvec.size(); iq++){
				double q_length_square = latt->GGT.metric_length_squared(wrap(qvec[iq])); // qvec is already wrapped to [-0.5,0.5)
				fprintf(fpvq, "%14.7le %14.7le %14.7le\n", q_length_square, abs(qscr2_static_RPA[iq]), abs(qscr2_2ndway[iq]));
			}
			fclose(fpvq);
		}

		if (ionode) printf("\ncalc_qscr2_static_RPA done\n");
	}

	// vq RPA routine

	void calc_vq_RPA(){
		vector<vector<complex>> qscr2_RPA(qvec.size());
		for (int iq = 0; iq < qvec.size(); iq++){
			int nomega = 1;
			if (clp.dynamic == "ppa" && clp.ppamodel == "gn") nomega = 2;
			if (clp.dynamic == "real-axis") nomega = omegaq[iq].size();
			qscr2_RPA[iq].resize(nomega);
			qscr2_RPA[iq][0] = qscr2_static_RPA[iq];
		}
		qscr2_static_RPA = vector<complex>();
		// vq = vq0 / (1 - vq0 * sum_k [(f_k - f_k-q) / (e_k - e_k-q - w - i0)] / nk_full)
		// vq0 = e^2 / V / (eps_r * eps_0) / q^2
		// Therefore, vq = e^2 / V / (eps_r * eps_0) / (q^2 + betas^2)
		// betas^2 = - e^2 / V / (eps_r * eps_0) * sum_k [(f_k - f_k-q) / (e_k - e_k-q - w - i0)] / nk_full
		for (int iq = 0; iq < qvec.size(); iq++){
			int nomega = 1;
			if (clp.dynamic == "ppa" && clp.ppamodel == "gn") nomega = 2;
			if (clp.dynamic == "real-axis") nomega = omegaq[iq].size();
			for (int iw = 1; iw < nomega; iw++){
				for (int ik = mp->varstart; ik < mp->varend; ik++){
					size_t jk = 0; vector3<> kj = elec->kvec[ik] - qvec[iq]; // not necessage to wrap k point around Gamma, kmap subroutines will wrap inside
					if (kmap->findk(kj, jk)){
						calc_ovlp(ik, jk);
						for (int b1 = 0; b1 < nb; b1++)
						for (int b2 = 0; b2 < nb; b2++){
							if (clp.scrFormula == "lindhard" && b1 != b2) continue;
							complex de;
							if (clp.dynamic == "ppa" && clp.ppamodel == "gn" && iw == 1)
								de = e[ik][b1] - e[jk][b2] - omega[iw];
							else
								de = e[ik][b1] - e[jk][b2] - omegaq[iq][iw];
							complex dfde = complex(f[jk][b2] - f[ik][b1], 0) / de;
							if (clp.scrFormula == "RPA") qscr2_RPA[iq][iw] += dfde * ovlp[b1*nb + b2].norm();
							else if (clp.scrFormula == "lindhard") qscr2_RPA[iq][iw] += dfde;
						}
					}
				}
				mp->allreduce(qscr2_RPA[iq][iw], MPI_SUM);
				qscr2_RPA[iq][iw] = complex(prefac_vq / nk_full, 0) * qscr2_RPA[iq][iw];
			}
		}
		
		if (clp.dynamic == "static" || clp.dynamic == "ppa"){
			vq_RPA.resize(qvec.size());
			for (size_t iq = 0; iq < qvec.size(); iq++){
				vq_RPA[iq].resize(1);
				double q_length_square = latt->GGT.metric_length_squared(wrap(qvec[iq])); // qvec is already wrapped to [-0.5,0.5)
				vq_RPA[iq][0] = complex(prefac_vq, 0) / (q_length_square + qscr2_RPA[iq][0]);
			}
		}
		if (clp.dynamic == "ppa"){
			double wp2 = clp.eppa * clp.eppa;

			/*
			//If we use qscr^2(w) = A / (w^2 - wq^2), we can include q=0 in Godby�CNeeds PPA
			Aq_ppa.resize(qvec.size()); Eq2_ppa.resize(qvec.size());

			for (size_t iq = 0; iq < qvec.size(); iq++){
			if (iq == 0 && clp.ppamodel == "hl") continue;

			double q_length_square = latt->GGT.metric_length_squared(wrap(qvec[iq])); // qvec is already wrapped to [-0.5,0.5)

			if (clp.ppamodel == "hl"){
			//Hybertsen - Louie
			double dtmp = wp2 * q_length_square;
			Aq_ppa[iq] = -dtmp * qscr2_RPA[iq][0];
			Eq2_ppa[iq] = dtmp / qscr2_RPA[iq][0];
			}
			else if (clp.ppamodel == "gn"){
			//Godby�CNeeds
			complex ctmp = wp2 * qscr2_RPA[iq][1];
			Aq_ppa[iq] = -ctmp * qscr2_RPA[iq][0];
			Eq2_ppa[iq] = ctmp / (qscr2_RPA[iq][0] - qscr2_RPA[iq][1]);
			}
			}
			*/

			//eps^-1(w) = 1 + A / (w^2 - wq^2)
			Aq_ppa.resize(qvec.size()); Eq2_ppa.resize(qvec.size());

			for (size_t iq = 0; iq < qvec.size(); iq++){
				double q_length_square = latt->GGT.metric_length_squared(wrap(qvec[iq])); // qvec is already wrapped to [-0.5,0.5)
				if (q_length_square < 1e-20) continue; //currently Gamma point is skipped

				if (clp.ppamodel == "hl"){
					//Hybertsen - Louie
					Aq_ppa[iq] = prefac_vq * wp2;
					Eq2_ppa[iq] = wp2 * (1 + complex(q_length_square, 0) / qscr2_RPA[iq][0]);
				}
				else if (clp.ppamodel == "gn"){
					//Godby�CNeeds
					complex eps0inv = c1 / (1 + qscr2_RPA[iq][0] / q_length_square);
					complex epspinv = c1 / (1 + qscr2_RPA[iq][1] / q_length_square);
					Eq2_ppa[iq] = wp2 * (1 - epspinv) / (epspinv - eps0inv);
					Aq_ppa[iq] = prefac_vq * (1 - eps0inv) * Eq2_ppa[iq];
				}
			}
		}
		if (clp.dynamic == "real-axis"){
			vq_RPA.resize(qvec.size());
			for (size_t iq = 0; iq < qvec.size(); iq++){
				vq_RPA[iq].resize(omegaq[iq].size());
				double q_length_square = latt->GGT.metric_length_squared(wrap(qvec[iq])); // qvec is already wrapped to [-0.5,0.5)

				for (size_t iw = 0; iw < omegaq[iq].size(); iw++){
					vq_RPA[iq][iw] = complex(prefac_vq, 0) / (q_length_square + qscr2_RPA[iq][iw]);

					if (ldebug && ionode && q_length_square > 1e-20){
						complex eps = 1 + qscr2_RPA[iq][iw] / complex(q_length_square, 0);
						complex eps_heg = heg->eps_intp(qvec[iq], omegaq[iq][iw].real());
						if (abs(eps) < 0.2 || (abs(eps) / abs(eps_heg) < 0.5)){
							printf("iq = %lu  w = %10.3le T  |q|^2 = %10.3le  qscr2 = %10.3le %10.3le  eps_heg = %10.3le %10.3le\n", 
								iq, omegaq[iq][iw].real() / T, q_length_square, qscr2_RPA[iq][iw].real(), qscr2_RPA[iq][iw].imag(), eps_heg.real(), eps_heg.imag());
						}
					}
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if (ionode){
			//static case is always written
			string fnamevq = "ldbd_vq.out";
			if (exists(fnamevq)) fnamevq = "ldbd_vq_updated.out";
			FILE *fpvq = fopen(fnamevq.c_str(), "w");
			init_model(clp.nfreetot, fpvq);
			fprintf(fpvq, "#|q|^2 |q_scr|^2 vq\n"); fflush(fpvq);
			for (size_t iq = 0; iq < qvec.size(); iq++){
				double q_length_square = latt->GGT.metric_length_squared(wrap(qvec[iq])); // qvec is already wrapped to [-0.5,0.5)
				fprintf(fpvq, "%14.7le %14.7le %14.7le\n", q_length_square, abs(qscr2_RPA[iq][0]), abs(vq_RPA[iq][0])); fflush(fpvq);
			}
			fclose(fpvq);

			if (clp.dynamic == "ppa"){
				string fnamevq = "ldbd_vq_ppa.out";
				if (exists(fnamevq)) fnamevq = "ldbd_vq_ppa_updated.out";
				FILE *fpvq = fopen(fnamevq.c_str(), "w");

				if (clp.ppamodel == "hl")
					fprintf(fpvq, "#|q|^2 |q_scr|^2 wp^2 |Eq|^2\n");
				else if (clp.ppamodel == "gn")
					fprintf(fpvq, "#|q|^2 |q_scr|^2 vq(w=0) |q_scr(i*wp)|^2 Aq Eq^2\n");
				for (size_t iq = 0; iq < qvec.size(); iq++){
					double q_length_square = latt->GGT.metric_length_squared(wrap(qvec[iq])); // qvec is already wrapped to [-0.5,0.5)
					if (clp.ppamodel == "hl"){
						fprintf(fpvq, "%14.7le %14.7le  ", q_length_square, abs(qscr2_RPA[iq][0]));
						fprintf(fpvq, "%14.7le %14.7le\n", wp2, abs(Eq2_ppa[iq])); fflush(fpvq);
					}
					else if (clp.ppamodel == "gn"){
						fprintf(fpvq, "%14.7le %14.7le %14.7le   %14.7le %14.7le   ", q_length_square, abs(qscr2_RPA[iq][0]), abs(vq(qvec[iq], 0)), qscr2_RPA[iq][1].real(), qscr2_RPA[iq][1].imag());
						fprintf(fpvq, "%14.7le %14.7le   %14.7le %14.7le\n", Aq_ppa[iq].real() / prefac_vq, Aq_ppa[iq].imag() / prefac_vq, Eq2_ppa[iq].real(), Eq2_ppa[iq].imag()); fflush(fpvq);
					}
				}
				fclose(fpvq);

				double wmax = 14 * T;
				int nw = 141;
				double dw = wmax / (nw - 1);
				vector<double> w(nw);
				w[0] = 0;
				for (int iw = 1; iw < nw; iw++)
					w[iw] = w[iw - 1] + dw;
				int iq = 0;
				double q_length_square = latt->GGT.metric_length_squared(wrap(qvec[iq]));
				if (q_length_square < 1e-20) iq = 1;
				q_length_square = latt->GGT.metric_length_squared(wrap(qvec[iq]));
				string fnamevqw = "ldbd_vqw_1stq.out";
				FILE *fpvqw = fopen(fnamevqw.c_str(), "wb");
				for (size_t iw = 0; iw < w.size(); iw++)
					fprintf(fpvqw, "%14.7le %14.7le %14.7le\n", w[iw] / T, abs(vq(qvec[iq], w[iw])), prefac_vq / q_length_square / abs(vq(qvec[iq], w[iw])));
				fclose(fpvqw);
			}
			else if (clp.dynamic == "real-axis"){
				vector<int> iq_test_arr{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, (int)round(qvec.size() / 4) - 1, (int)round(qvec.size() / 2) - 1, (int)qvec.size() - 1 };
				for (int iqt = 0; iqt < iq_test_arr.size(); iqt++){
					int iq = iq_test_arr[iqt];
					double q_length_square = latt->GGT.metric_length_squared(wrap(qvec[iq])); // qvec is already wrapped to [-0.5,0.5)

					//directly calculated qscr2 and vqw
					string fnamevqw = "vqw_q" + int2str(iq) + ".out";
					FILE *fpvqw = fopen(fnamevqw.c_str(), "w");
					fprintf(fpvqw, "#w (kBT) qscr2 vqw (smearing = %7.3lf kBT |q|^2 = %14.7le)\n", clp.smearing / T, q_length_square);
					for (size_t iw = 0; iw < omegaq[iq].size(); iw++)
						fprintf(fpvqw, "%10.3le   %10.3le %10.3le   %10.3le %10.3le\n", omegaq[iq][iw].real() / T, qscr2_RPA[iq][iw].real(), qscr2_RPA[iq][iw].imag(), vq_RPA[iq][iw].real(), vq_RPA[iq][iw].imag());
					fclose(fpvqw);

					//interpolated quantities
					double wmax = wqmax[iq] + omegaq[iq][1].real();
					int nw = (int)round(sqrt(3)*omegaq[iq].size());
					double dw = wmax / (nw - 1);
					vector<double> w(nw);
					w[0] = 0; w[nw-1] = wmax;
					for (int iw = 1; iw < nw-1; iw++)
						w[iw] = iw*dw;

					fnamevqw = "vqw_intp_q" + int2str(iq) + ".out";
					fpvqw = fopen(fnamevqw.c_str(), "wb");
					fprintf(fpvqw, "#w (kBT) |vqw| ReEps ImEps |Eps| (smearing = %7.3lf kBT |q|^2 = %14.7le)\n", clp.smearing / T, q_length_square);
					for (size_t iw = 0; iw < w.size(); iw++){
						complex vqw = vq(qvec[iq], w[iw]);
						complex eps = complex(prefac_vq / q_length_square, 0) / vqw;
						fprintf(fpvqw, "%10.3le %10.3le %10.3le %10.3le %10.3le\n", w[iw] / T, abs(vq(qvec[iq], w[iw])), eps.real(), eps.imag(), eps.abs());
					}
					fclose(fpvqw);
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

    // compute vq

    complex vq(vector3<double> q, double w = 0) override {
		double q_length_square = latt->GGT.metric_length_squared(wrap(q));
		complex result;
		if (clp.scrFormula == "unscreened"){
			if (q_length_square < 1e-20) return c0; // skip Gamma point in current version
			return complex(prefac_vq / q_length_square, 0);
		}
		else if (clp.scrFormula == "heg"){
			if (q_length_square < 1e-20) return c0; // skip Gamma point in current version
			//if (ionode) printf("q = %lg  w = %lg in T\n", sqrt(q_length_square), w/T);
			complex eps = heg->eps_intp(q, w);
			//if (ionode) printf("|eps| = %lg\n", abs(eps));
			if (abs(eps) == 0)
				error_message("eps is zero!","vq");
			else{
				if (abs(eps) < 1e-20) printf("|eps| = %10.3le is too tiny!\n", abs(eps));
				return complex(prefac_vq,0) / eps / q_length_square;
			}
		}
		else if (clp.scrFormula == "debye"){
			return complex(prefac_vq / (qscr2_debye + q_length_square), 0);
		}
		else if (clp.scrFormula == "Bechstedt"){
			if (q_length_square < 1e-20) return c0; // skip Gamma point in current version
			double q2 = q_length_square;
			double epsq = 1 + 1. / (fac0_Bechstedt + fac2_Bechstedt * q2 + fac4_Bechstedt * q2 * q2);
			return complex(prefac_vq_bare / epsq / q2, 0);
		}
		else{
			size_t iq = qmap->q2iq(q);
			if (q_length_square < 1e-20) return c0; // skip Gamma point in current version
			if (clp.dynamic == "static" || w == 0)
				return vq_RPA[iq][0].real();
			else if (clp.dynamic == "ppa"){
				double q_length_square = latt->GGT.metric_length_squared(wrap(q));
				//if (clp.ppamodel == "hl"){
				complex cw = complex(fabs(w), clp.smearing);
				result = (prefac_vq + Aq_ppa[iq] / (cw*cw - Eq2_ppa[iq])) / q_length_square;
				//}
				//else if (clp.ppamodel == "gn"){ //for GN ppa, smearing is not needed
				//	result = (prefac_vq + Aq_ppa[iq] / (w*w - Eq2_ppa[iq])) / q_length_square; // Aq_ppa has been multiplied by prefac_vq
				//}
			}
			else if (clp.dynamic == "real-axis"){
				double dw = omegaq[iq][1].real();
				int iw = floor(fabs(w) / dw);
				int nw = omegaq[iq].size();
				if (iw >= nw - 1) result = vq_RPA[iq][nw - 1];
				else result = (vq_RPA[iq][iw] * (omegaq[iq][iw + 1].real() - w) + vq_RPA[iq][iw + 1] * (w - omegaq[iq][iw].real())) / dw; // linear interpolation
			}
			if (w < 0) return result.conj(); // eps(q,-w)=eps(q,w)^*
			else return result;
		}
		return 0;
	}

    void print() override {
        std::cout << "3D COULOMB CLASS" << std::endl;
    }

};