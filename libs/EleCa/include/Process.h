#ifndef _PROCESS_H_
#define _PROCESS_H_

#include <string.h>

using namespace std;
namespace crpropa {
  class Process {

  public: 
    
    std::string fname;
    double flambda;
    double fsmin; 
    double fsmax;
    double fCMEnergy;
    double fInteractionAngle;
    Particle fPi;
    Particle fPt;
  
    std::string fback;
    double fbackdensity;
    
    Process();
    Process(Process&);
    Process(Particle&, Particle&);
    Process(Particle&, Particle&, string);
    
    ~Process();
    
    void SetName(std::string nm) {fname = nm;}
    std::string GetName() {return fname;}
    
    void SetInteractionAngle(double a) {fInteractionAngle = a;} 
    double GetInteractionAngle() {return fInteractionAngle;} 
    
    void SetLambda(double le) {flambda = le;}
    double GetLambda() {return flambda;} 
    
    void SetLimits(double smin, double smax) {fsmin = smin; fsmax = smax;}
    void SetLimits(Particle& p1, std::string nameproc);
    void SetLimits() {SetLimits(fPi,fname);}
    
    void SetMax(double smax) {fsmax = smax;}
    void SetMin(double smin) {fsmin = smin;}
    double GetMin() {return fsmin;}
    double GetMax() {return fsmax;}
    
    void SetCMEnergy(double s) {fCMEnergy = s;}
    
    void SetCMEnergy(Particle p1, Particle pb) {
      fCMEnergy = 2*p1.GetEnergy()*pb.GetEnergy()*(1-p1.GetBeta()*cos(fInteractionAngle)) + p1.GetMass()*p1.GetMass() + pb.GetMass()*pb.GetMass();
    }
    
    void SetCMEnergy() {
      fCMEnergy = 2*fPi.GetEnergy()*fPt.GetEnergy()*(1-fPi.GetBeta()*cos(fInteractionAngle)) + fPi.GetMass()*fPi.GetMass() + fPt.GetMass()*fPt.GetMass();
      
    }
    
    double GetCMEnergy() {return fCMEnergy;}
    
    void SetIncidentParticle(Particle& p1) {fPi = p1; SetLimits();}
    void SetTargetParticle(Particle& p1) {fPt = p1;}
    
    Particle GetIncidentParticle() {return fPi;}
    Particle GetTargetParticle() {return fPt;}
    
    std::string GetBackground() {return fback;}
    void SetBackground(std::string BackRad);
  
  private:

    
  };
  
//-------------- explicit definitions

Process::Process() {
  fname="";
  SetLimits(0.0,1.0e23);
  flambda =0;
  fCMEnergy =0;
  fInteractionAngle=cPI;
  fback="ALL";
  fbackdensity=0;
}


Process::Process(Particle& p1, Particle& p2) {
  fPi = p1;
  fPt = p2;
  if (p1.GetType()==22)
    fname = "PP";
  else if (abs(p1.GetType())==11) { 
    cerr << "NB: by default process set to ICS" << endl;
    fname="ICS";
  }
  else fname="NONE";
  SetCMEnergy(p1, p2);
  flambda=0;
  fInteractionAngle= cPI;
  fback="ALL";
  SetLimits(p1, fname);
  fbackdensity=0;
}

 Process::Process(Particle& p1, Particle& p2, std::string name) {
  fname = name;
  SetCMEnergy(p1, p2);
  flambda=0;
  fInteractionAngle= cPI;
  fPi = p1;
  fPt = p2;
  fback="ALL";
  SetLimits(p1, fname);
  fbackdensity=0;
}

Process::Process(Process& proc2) {
  fname = proc2.GetName();
  SetLimits(proc2.GetMin(), proc2.GetMax());
  fCMEnergy = proc2.GetCMEnergy();
  fInteractionAngle = proc2.GetInteractionAngle();
  fPi = proc2.GetIncidentParticle();
  fPt = proc2.GetTargetParticle();
  fback = proc2.GetBackground();
  fbackdensity=0;
}

Process::~Process() {
}
 

//-----------

 void Process::SetBackground(std::string BackRad) {

    fback=BackRad;

    double eps_min = eps_ph_inf_global;
    double eps_max = eps_ph_sup_global;

   if( BackRad=="CMB" ){
        eps_min = eps_ph_inf_cmb;
        eps_max = eps_ph_sup_cmb;
	cout << "eps range setted to " <<  eps_min << " , " << eps_max <<endl; 
    }
   else if( BackRad=="COB" ){
        eps_min = eps_ph_inf_cob;
        eps_max = eps_ph_sup_cob;
	cout << "eps range setted to " <<  eps_min << " , " << eps_max <<endl; 
    }
    else if( BackRad=="CIB" ){
        eps_min = eps_ph_inf_cib;
        eps_max = eps_ph_sup_cib;
	cout << "eps range setted to " <<  eps_min << " , " << eps_max <<endl; 
    }
    else if( BackRad=="CIOB" ){
        eps_min = eps_ph_inf_ciob;
        eps_max = eps_ph_sup_ciob;
	cout << "eps range setted to " <<  eps_min << " , " << eps_max <<endl; 
    }
   else if( BackRad=="URB" ){
        eps_min = eps_ph_inf_urb;
        eps_max = eps_ph_sup_urb;
	cout << "eps range setted to " <<  eps_min << " , " << eps_max <<endl; 
    }
  
   feps_inf=eps_min;
   feps_sup=eps_max;

}

 void Process::SetLimits(Particle& p1, std::string nameproc) {
   if (p1.GetType()!=22 && p1.GetType()!=11 && p1.GetType()!=-11) cout << "error in type " << p1.GetType() << " != 11 and !=22 " << endl;

   if (nameproc == "PP") {
     if (abs(p1.GetType())==11) 
       cout << "\nERROR!! wrong particle or process!! " << nameproc << p1.GetType() << "\n" << endl;
     fsmin = 4*ElectronMass*ElectronMass; //min per theta = 0;
     fsmax = 4*p1.GetEnergy()*feps_sup*2.0;//(1.0-cos(fInteractionAngle)); //max per theta=3.14
   }
   if (nameproc == "DPP") {
     if (abs(p1.GetType())==11) 
       cout << "\nERROR!! wrong particle or process!! " << nameproc << p1.GetType() << "\n" << endl;
     fsmin = 16*ElectronMass*ElectronMass;
     fsmax = 4*p1.GetEnergy()*feps_sup*2.0;
   }
    if (nameproc == "ICS") {
      if (p1.GetType()==22 || p1.GetType()==9) 
       cout << "\nERROR!! wrong particle or process!! " << nameproc << p1.GetType() << "\n" << endl;
      fsmin = p1.GetMass()*p1.GetMass() + 2*p1.GetEnergy()*feps_inf*(1 - p1.GetBeta()); 
      fsmax = p1.GetMass()*p1.GetMass() + 2*p1.GetEnergy()*feps_sup*(1 + p1.GetBeta()) ;
    }
    if (nameproc == "TPP") {
      if (p1.GetType()==22 || p1.GetType()==9)
       cout << "\nERROR!! wrong particle or process!! " << nameproc << p1.GetType() << "\n" << endl;
      fsmin = 10*ElectronMass*ElectronMass;
      fsmax = 2*p1.GetEnergy()*feps_sup*(1 + p1.GetBeta()) + p1.GetMass()*p1.GetMass() ;
    }
 }

}


#endif
