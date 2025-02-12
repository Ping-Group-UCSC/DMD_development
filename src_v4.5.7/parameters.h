#pragma once
#include "common_headers.h"

class parameters{
  public:
  	bool restart; //!< restart from previous calculation, needs restart directory
    bool compute_tau_only;
    bool print_tot_band;
    bool print_along_kpath;
  	std::vector<vector3<double>> kpath_start, kpath_end;
  	int freq_measure, freq_measure_ene, freq_compute_tau, freq_update_eimp_model, freq_update_ee_model;
  	double de_measure, degauss_measure;
  	double t0, tend, tstep, tstep_laser;
  	int nk1, nk2, nk3;
  	double ewind;
  	double temperature;
  	double degauss, ndegauss;
  	double mu, carrier_density; 
    bool carrier_density_means_excess_density;
  	double scale_scatt, scale_eph, scale_ei, scale_ee; //!< scaling factors for scattering
  	int modeStart, modeEnd; //!< currectly only affect e-ph analysis part
  	int band_skipped; /*!< dft bands skipped in wannierization
  	                   relevant when computing internal magnetic fields */
  	double degthr;
  	bool rotate_spin_axes;
  	vector3<double> sdir_z, sdir_y, sdir_x; //!< redefined spin directions
  	matrix3<double> sdir_rot;
  	bool need_imsig;
  	double scissor; //!< Apply scissor operation to bandgap units eV
  
  	double tau_phenom; //!< phenominological relaxation
  	int bStart_tau, bEnd_tau; //!< if phenom_tau true
  
  	double Bx, By, Bz, scale_Ez;
    vector3<> B; //!< External B field in Tesla defined by Bx, By, Bz

  	double Bxpert; //!< B perterbation for initial spin inbalance in Tesla x
    double Bypert; //!< B perterbation for initial spin inbalance in Tesla y
    double Bzpert; //!< B perterbation for initial spin inbalance in Tesla z
    vector3<> Bpert; //!< to generate an initial spin inbalance units Tesla defined by B*pert
  	bool gfac_normal_dist; bool gfac_k_resolved;
  	double gfac_mean, gfac_sigma, gfac_cap;
  	bool needL;
  
  	vector3<double> lattvec1, lattvec2, lattvec3;
  	matrix3<> R;
  	int dim;
  	double thickness;
  	std::vector<vector3<>> vpos;
  	std::vector<std::vector<bool>> vtrans;
  	string type_q_ana;
  
  	parameters(){}
  
  	void read_param();
  	void read_jdftx();
  	void get_valley_transitions(string file_forbid_vtrans);
  
  	std::string trim(std::string s);
  	std::map<std::string, std::string> map_input(fstream& fin);
  private:
  	double get(std::map<std::string, std::string> map, string key, double defaultVal = NAN, double unit = 1) const;
  	vector3<> getVector(std::map<std::string, std::string> map, string key, vector3<> defaultVal = vector3<>(NAN), double unit = 1) const;
  	string getString(std::map<std::string, std::string> map, string key, string defaultVal = "") const;
};
