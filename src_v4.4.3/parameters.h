#pragma once
#include "common_headers.h"

class parameters{
public:
	bool restart, compute_tau_only, print_along_kpath;
	std::vector<vector3<double>> kpath_start, kpath_end;
	int freq_measure, freq_measure_ene, freq_compute_tau, freq_update_eimp_model, freq_update_ee_model;
	double de_measure, degauss_measure;
	double t0, tend, tstep, tstep_pump;
	int nk1, nk2, nk3;
	double ewind;
	double temperature;
	double degauss, ndegauss;
	double mu, carrier_density;

	double tau_phenom;
	int bStart_tau, bEnd_tau; // if phenom_tau true

	string pumpMode;

	double Bx, By, Bz; vector3<> B;
	double Bxpert, Bypert, Bzpert; vector3<> Bpert; // to generate an initial spin inbalance

	vector3<double> lattvec1, lattvec2, lattvec3;
	matrix3<> R;
	int dim;
	double thickness;

	parameters(){}

	void read_param();
	void read_jdftx();

	std::string trim(std::string s);
	std::map<std::string, std::string> map_input(fstream& fin);
private:
	double get(std::map<std::string, std::string> map, string key, double defaultVal = NAN, double unit = 1) const;
	vector3<> getVector(std::map<std::string, std::string> map, string key, vector3<> defaultVal = vector3<>(NAN), double unit = 1) const;
	string getString(std::map<std::string, std::string> map, string key, string defaultVal = "") const;
};