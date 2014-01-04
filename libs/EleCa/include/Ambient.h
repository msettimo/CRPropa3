
#include <string>

using namespace std;

class Ambient { 

 public: 
  Ambient();
  Ambient(string);
  ~Ambient(); 

  string fback; 
  double feps;
  double fn_eps;
  double feps_inf;
  double feps_sup;
  double fInteractingEps;

  void SetBackground(string back) { 
    fback = back; 
    SetMin(fback);
  SetMax(fback); }

  string GetBackground() {return fback;}
  double GetDensity() {return fn_eps;}
  double GetEnergy() {return feps; }

  void SetMin(double min);
  void SetMax(double max);

  void SetMin(string BackRad);
  void SetMax(string BackRad);
  void SetLimits(string BackRad);

  double GetMin() {return feps_inf;}
  double GetMax() {return feps_sup;}

  void ExtractDensityBackground();
  void ExtractBackgroundEnergy();
  void ExtractTheta();

 private:
  
};


Ambient::Ambient() {
  fback="ALL";
  feps_inf = eps_ph_inf_global;
  feps_sup = eps_ph_sup_global;
  feps=feps_inf;
  fn_eps=0;
}



Ambient::Ambient(string BackRad) {
  fback=BackRad;
  SetMin(BackRad);
  SetMax(BackRad);
  feps=feps_inf;
  fn_eps=0;
}

Ambient::~Ambient() {
}


void Ambient::SetLimits(string BackRad) {
  SetMin(BackRad);
  SetMax(BackRad);
}

void Ambient::SetMax(string BackRad) {
  fback=BackRad;
  if (BackRad == "ALL" || BackRad == "") 
    feps_sup = eps_ph_sup_global;
  else if (BackRad == "CMB") 
    feps_sup = eps_ph_sup_cmb;
  else if (BackRad == "CIOB") 
    feps_sup = eps_ph_sup_ciob;
  else if (BackRad == "COB") 
    feps_sup = eps_ph_sup_cob;
  else if (BackRad == "CIB") 
    feps_sup = eps_ph_sup_cib;
  else if (BackRad == "URB") 
    feps_sup = eps_ph_inf_cmb;
  else if (BackRad == "CMIRB") 
    feps_sup = eps_ph_sup_ciob;

 }



void Ambient::SetMin(string BackRad) {
  if (BackRad == "ALL" || BackRad == "") 
    feps_inf = eps_ph_inf_global;
  else if (BackRad == "CMB") 
    feps_inf = eps_ph_inf_cmb;
  else if (BackRad == "CIOB") 
    feps_inf = eps_ph_inf_ciob;
  else if (BackRad == "COB") 
    feps_inf = eps_ph_inf_cob;
  else if (BackRad == "CIB") 
    feps_inf = eps_ph_inf_cib;
  else if (BackRad == "URB") 
    feps_inf = eps_ph_inf_urb;
 }

