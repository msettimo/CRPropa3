#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <errno.h>
using namespace std;
#include <sstream>

namespace crpropa {

double cPI = 4*atan(1);

 double Bfield;

 int POINTS_VERY_FEW =100;
 string BACKGROUNDRAD;


static const int MC_SAMPLING = 1000;    

double (*MC_Sampling_Hist)[3] = new double [MC_SAMPLING][3];
#ifndef __RANDOM__
#define __RANDOM__
    #include "RandomLib/rvgs.c"
    #include "RandomLib/rngs.c"
#endif /*__RANDOM__*/

const char *VERSION = "0.1";

//#####################
//# Runge-Kutta
//#####################
// Needed for the optimized version of the propagation code
static const int RK_ORDER = 6;
double (*vRKa) = new double [RK_ORDER+1];
double (*vRKc) = new double [RK_ORDER+1];
double (*vRKcs) = new double [RK_ORDER+1];
double (*vRKb)[RK_ORDER+1] = new double [RK_ORDER][RK_ORDER+1];        

static const double M2MPC = 3.240779e-17/1.0e6;
static const double KM2MPC =  3.240779e-20;
static const double S2YR = 3.168808781403e-08;
static const double H0 = 70.4;   
static const double H0y = H0*(double)KM2MPC/S2YR;   // H0 per years
double OC = 0.227;  // Dark matter density
double OB = 0.0456; // Baryon density
double OM = OB+OC;  // Matter density
double OL = 0.728;  // Dark energy density

static const double K_CBR = 1.318684673251832e+13; // =1/pi**2/hbarc**3 [eV^-3 cm^-3]
static const double eV2J = 1.602176487e-19; // from eV to J
static const double ElectronMass = 0.510998918e6; // [eV/c^2]
static const double K_boltz = 8.617342294984e-5;  // [eV/K ] Boltzman constant
static const double C_speed = 299792458;          // [m/s] speed of light
static const double SigmaThompson = 6.6524e-25;

static const double T_CMB = 2.725; // [K] // evolution 2.725*(1-z)  1012.3164 
static const double T_COB = 5270; // [K]  // Visible [380 - 760] nm. Scelgo 550
static const double T_CIB = 1.45e+02; // [k] Middle IR 5 to (25-40) µm according to Nasa. scelgo 20e-6 m
static const double T_CRB = 3e-03; // [k] ~ cm - 10m.  scelgo ~1 m

static const double CMB_en = K_boltz*T_CMB; //2.348e-4;             // [eV]
static const double CRB_en = K_boltz*T_CRB; 
static const double COB_en = K_boltz*T_COB; 
static const double CIB_en = K_boltz*T_CIB; 
static const double CIOB_en = CIB_en+COB_en;    // [eV]

static const double h_Planck = 4.135667e-15; // [eV s]// 
static const double hcut_Planck = h_Planck/2/cPI; // [eV s] hcut = h/2Pi [Js]
static const double LambdaCompton = hcut_Planck/(ElectronMass/C_speed); 

static const double eps_ph_inf_urb = 4.1e-12;   // [eV]
static const double eps_ph_inf_cmb = 0.825e-6;   // [eV]
static const double eps_ph_inf_cib = 2e-3;   // [eV]
static const double eps_ph_inf_cob = 5e-2;   // [eV]
static const double eps_ph_inf_ciob = 2e-3;   // [eV]

static const double eps_ph_sup_urb = eps_ph_inf_cmb;//4e-5;   // [eV]
static const double eps_ph_sup_cmb = eps_ph_inf_cob;   // [eV]
static const double eps_ph_sup_cob = 9.9;   // [eV]
static const double eps_ph_sup_cib = 0.8;   // [eV]
static const double eps_ph_sup_ciob = 9.9;   // [eV]


static const double eps_ph_sup_global = eps_ph_sup_cob;   // [eV] *global
static const double eps_ph_inf_global = eps_ph_inf_urb;   // [eV] *global

int NsecG=0;
vector<double> EGround;
vector<double> zGround;
vector<int> typeGround;

double E0taborg = 0;

double z0ph = 0;
double E0ph = 0;

int particle_type = 0;
double EnergyCM =0;
double GammaEnergy = 0;
double BackGamma = 0;
double PPxsection = 0;
double DPPxsection = 0;
double TPPxsection = 0;
double ICSxsection = 0;
double PPlength = 0;
double DPPlength = 0;
double TPPlength = 0;
double ICSlength = 0;
double n_eps2 = 0;
double eps2 = 0;
double feps_inf =  eps_ph_inf_global;
double feps_sup =  eps_ph_sup_global;


bool debug;

//##########################################################################
//# Functions
//##########################################################################

double z2Mpc(double z){
    if( z<0.4 ){
        return ((double)C_speed/1000./H0)*z;
    }else{
        // AV Uryson, Physics Particles and Nuclei, 2006, Vol. 37, No. 3, pp. 347   67
        // Assuming flat-matter-dominated cosmology. Error is negligible
        return ((double)C_speed/1000./H0)*(2./3.)*(1-pow(1.+z,-1.5));
    }
};
    
double Mpc2z(double D){
    if( D<1700. ){
        return (double)D/((double)C_speed/1000./H0);
    }else{
        // AV Uryson, Physics Particles and Nuclei, 2006, Vol. 37, No. 3, pp. 347   67
        // Assuming flat-matter-dominated cosmology. Error is negligible
        return pow(1-(double)D/((2./3.)*(double)C_speed/1000./H0),-2./3.)-1;
    }
};


}
