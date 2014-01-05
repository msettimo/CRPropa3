#ifndef MPC_PHOTON_ELECA_H_
#define MPC_PHOTON_ELECA_H_

#include <vector>
#include "crpropa/Module.h"
#include "crpropa/Candidate.h"
#include "crpropa/magneticField/MagneticField.h"
#include "../libs/EleCa/include/Constants.h"
#include "../libs/EleCa/include/Particles.h"
#include "../libs/EleCa/include/Process.h"
#include "../libs/EleCa/include/XLoss_CBR.h"
#include "../libs/EleCa/include/EnergyLoss.h"


using namespace std;
using namespace crpropa;

void ReadTables(string bkg);
void InitBkgArray(string BackRad);
  
static vector< vector<double> > IntLengths;
static vector< vector<double> > BkgArray0;

class PhotonEleCa {
    
  private:

    ref_ptr<MagneticField> field; //for future use? 
 
     double vPPle[1101];
     double vDPPle[1101];
     double vTPPle[1101];
     double vICSle[1101];
     double vEtab[1101];

     double BkgArray[100][2];
     double fEthr;
     Candidate c;

  public:    
     vector<int> fdN;
    PhotonEleCa(); 
    PhotonEleCa(Candidate c);
    PhotonEleCa(int type, double E0, double z0) ;
    ofstream outfile;

    vector<int> GetFlux() { return fdN; };    
    void SetEthr(double eth) { fEthr = eth; };
    double GetEthr() { return fEthr; };    
    ~PhotonEleCa() {};
    
    void process(Candidate candidate, double ethr) ;
    void WriteOutput(Particle &p1, vector<Particle> &part, bool spectropt=0);

    std::string getDescription();

    double GetMeanThetaBFDeflection(double Bin, double Ein, int ptype, double Lin);
    double GetLambdaTab(Process proc, string procName);
    double ExtractMinDist(Process &proc, int type, double R, double R2, vector<double> Etarget);
    vector<double> GetEtarget(Process &proc, Particle &particle);     
    void Propagate(Particle &curr_particle, vector<Particle> &ParticleAtMatrix, vector<Particle> &ParticleAtGround);
    double ExtractPhotonEnergyMC(double z, Process &proc);
    double ShootPhotonEnergyMC(double z);
    void SetInitVar(vector< vector<double> > bk, vector< vector<double> > *le);

};


 PhotonEleCa::PhotonEleCa(crpropa::Candidate c){ 
      double zt= c.getRedshift();
      double et = c.current.getEnergy();
      int id = c.current.getId();
     
      Particle p0(id,et,zt);
      for (int j=0; j<1101; ++j) {
	vEtab[j]=IntLengths.at(j).at(0);
	vPPle[j]=IntLengths.at(j).at(1);
	vDPPle[j]=IntLengths.at(j).at(2);
	vTPPle[j]=IntLengths.at(j).at(3);
	vICSle[j]=IntLengths.at(j).at(4);
      }
      
      for (int k=0; k<100; k++) { 
	
	BkgArray[k][0]=BkgArray0[k][0];
	BkgArray[k][1]=BkgArray0[k][1];
      }
      
      fEthr=1e16;
    };



 PhotonEleCa::PhotonEleCa() {
     for (int j=0; j<1101; ++j) {
	vEtab[j]=IntLengths.at(j).at(0);
	vPPle[j]=IntLengths.at(j).at(1);
	vDPPle[j]=IntLengths.at(j).at(2);
	vTPPle[j]=IntLengths.at(j).at(3);
	vICSle[j]=IntLengths.at(j).at(4);
      }
      
      for (int k=0; k<100; k++) { 

	BkgArray[k][0]=BkgArray0[k][0];
	BkgArray[k][1]=BkgArray0[k][1];
      }
     
      fEthr=1e16;
      outfile.open("test1.txt",ios::out);
      if (!outfile.is_open()) {
	cout <<  "Unable to open output file! exiting... " << endl;
	return;
      }
   };


 PhotonEleCa::PhotonEleCa(int type, double E0, double z0) {
      Particle p0(type,E0,z0);
    
      for (int j=0; j<1101; ++j) {
	vEtab[j]=IntLengths.at(j).at(0);
	vPPle[j]=IntLengths.at(j).at(1);
	vDPPle[j]=IntLengths.at(j).at(2);
	vTPPle[j]=IntLengths.at(j).at(3);
	vICSle[j]=IntLengths.at(j).at(4);
      }
      
      for (int k=0; k<100; k++) { 

	BkgArray[k][0]=BkgArray0[k][0];
	BkgArray[k][1]=BkgArray0[k][1];
      }
    
      fEthr=1e16;
    };


void ReadTables(string nameback) {
    string filename= "InteractionLengths_lee.dat";
    string filename2 = "../libs/EleCa/include/data/"+nameback+filename;
  
    cout << filename2.c_str() << endl;
    
    ifstream fin(filename2.c_str());
    
    if (!fin.is_open()) { 
      cout <<  "Unable to open lambda_table file: " << filename2.c_str() << " ! exiting... "; 
      return;
    } else   
    cout << " file " << filename2.c_str() << " opened! " << endl;
    
    int k=0;
    double Etab, PPle, ICSle, DPPle, TPPle;
    vector<double> tmp;
    while (fin.good()) {
      fin >> Etab >> PPle >> ICSle >> DPPle >> TPPle;
      tmp.clear();
      tmp.push_back(Etab); 
      tmp.push_back(PPle); 
      tmp.push_back(DPPle); 
      tmp.push_back(ICSle); 
      tmp.push_back(TPPle); 
      IntLengths.push_back(tmp);
     }
  
   
    double dEtab = log10(IntLengths.at(0).at(0)) - log10(IntLengths.at(1).at(0)); 
    if (fin.is_open()) fin.close(); 
};


void InitBkgArray(string BackRad){
    // Routine to build the array of cumulative distribution of
    // background photons
  
    int i = 0;
    double epsinf=0; double epssup=0;
    vector<double> tmp;
    tmp.resize(2);
    if( BackRad=="CMB" ){
        BACKGROUNDRAD = "CMB";

        double de = pow( (double)eps_ph_sup_cmb/eps_ph_inf_cmb, 1./100.);        
        for(double e=eps_ph_inf_cmb; e<=eps_ph_sup_cmb; e*=de){
	  tmp.at(0) = e;
            if(i==0)
	      tmp.at(1)=CMBR(e);
            else
	      tmp.push_back(BkgArray0.at(i-1).at(1) + CMBR(e));
            
	    BkgArray0.push_back(tmp);
            i++;
        }
    }

    
    if( BackRad=="CIOB" ){
        BACKGROUNDRAD = "CIOB";
        
        double de = pow( (double)eps_ph_sup_ciob/eps_ph_inf_ciob, 1./100.);
        
        for(double e=eps_ph_inf_ciob; e<=eps_ph_sup_ciob; e*=de){
	    tmp.clear();
	    tmp.push_back(e);
            if(i==0)
	      tmp.push_back(CIOBR(e));
            else 
	      tmp.push_back(BkgArray0.at(i-1).at(1) + CIOBR(e));
            
	    BkgArray0.push_back(tmp);
            i++;
        }
    }


    if( BackRad=="URB" ){
        BACKGROUNDRAD = "URB";
        double de = pow( (double)eps_ph_sup_urb/eps_ph_inf_urb, 1./100.);
        for(double e=eps_ph_inf_urb; e<=eps_ph_sup_urb; e*=de){
	  tmp.clear();
	  tmp.push_back(e);

	  if(i==0)
	      tmp.push_back(URB(e));               
	  else
	    tmp.push_back(BkgArray0.at(i-1).at(1) + URB(e));
	
	    BkgArray0.push_back(tmp);
            i++;
        }
    }
    
    if( BackRad!="CMB" && BackRad!="CIOB" ){
        BACKGROUNDRAD = "ALL";
        
        double de = pow( (double)eps_ph_sup_global/eps_ph_inf_global, (double)1./POINTS_VERY_FEW);
        for(double e=eps_ph_inf_global; e<=eps_ph_sup_global; e*=de){
	  tmp.at(0) = e;
	  if(i==0)
	    tmp.at(1)=CBR(e);
	  else
	    tmp.at(1)=BkgArray0.at(i-1).at(1) + CBR(e);
	  BkgArray0.push_back(tmp);
	  i++;
	}
    }


    double a = 1.0/BkgArray0.at(i-1).at(1);

    for(i=0; i<POINTS_VERY_FEW; i++)
       	BkgArray0.at(i).at(1) *=a;
  };


//----- 

double PhotonEleCa::GetMeanThetaBFDeflection(double Bin, double Ein, int ptype, double Lin) {
  //from D. Hooper, S. Sarkar, M. Taylor, arXiv: 0608085, 2006
  
  if (Bin ==0 || Ein==0) return 0;
  if (ptype == 22) return 0;

  double lcoher=1;
 
  return 0.8*(1.0e20/Ein)*sqrt(Lin/10*lcoher)*(Bin/1.0e-9)/180.0*3.1415927;
  };

  
double PhotonEleCa::ExtractMinDist(Process &proc, int type, double R, double R2, vector<double> Etarget) {
    
    double min_dist1=0; double min_dist2=0;
    Process proc1(proc);
    Process proc2(proc);
    double tmp_lambda1=0; double tmp_lambda2=0;
    Particle pt;
    pt.SetType(0);
    pt.Setz(proc.GetIncidentParticle().Getz());
    
    if (type==22) {
    
      proc1.SetName("PP");
      pt.SetEnergy(Etarget[0]);
      proc1.SetTargetParticle(pt);
      proc1.SetCMEnergy();
    
      tmp_lambda1 = GetLambdaTab(proc1,"PP");
    
      min_dist1 = -tmp_lambda1*log(R);
             
      proc2.SetName("DPP");
      pt.SetEnergy(Etarget[1]);
      proc2.SetTargetParticle(pt);
      proc2.SetCMEnergy();
      tmp_lambda2 = GetLambdaTab(proc2,"DPP");
      min_dist2 = -tmp_lambda2*log(R2);
      
      if (debug) 
	cerr << "comparing 2 mindists: " << min_dist1 << "(" << tmp_lambda1 << ") vs " << min_dist2 << " ( " << tmp_lambda2 << ") " << endl;  

       if (min_dist2 < min_dist1) {
         min_dist1 = min_dist2;
         proc.SetName("DPP");
         pt.SetEnergy(Etarget[1]);
         proc.SetTargetParticle(pt);
         proc.SetCMEnergy();
       }
       else {
         proc.SetName("PP");
         pt.SetEnergy(Etarget[0]);
         proc.SetTargetParticle(pt);
         proc.SetCMEnergy();
       }
     }//end if type 0
    else if (abs(type)==11) {

       proc1.SetName("ICS");
       pt.SetEnergy(Etarget[0]);
       proc1.SetTargetParticle(pt);
       tmp_lambda1 = GetLambdaTab(proc1,"ICS");
       min_dist1 = -tmp_lambda1*log(R);

       proc2.SetName("TPP");
       pt.SetEnergy(Etarget[1]);
       proc2.SetTargetParticle(pt);
       tmp_lambda2 = GetLambdaTab(proc2,"TPP");
       min_dist2 = -tmp_lambda2*log(R2);

     if (min_dist2 < min_dist1) {
         min_dist1 = min_dist2;
         proc.SetName("TPP");
	 pt.SetEnergy(Etarget[1]);
	 proc.SetTargetParticle(pt);
         proc.SetCMEnergy();
     }
     else {
       proc.SetName("ICS");
       pt.SetEnergy(Etarget[0]);
       proc.SetTargetParticle(pt);
       proc.SetCMEnergy();
       }
   }//else e+/e-
   else cerr << "something wrong in particle type ( " << type << ". Propagation of photons and e+/e- is the only allowed.)" << endl;
   
   return min_dist1;
};


  
  //-------
  
double PhotonEleCa::GetLambdaTab(Process proc, string procName) {
    
    double E1 = proc.GetIncidentParticle().GetEnergy();
    double z = proc.GetIncidentParticle().Getz();
    double res=0;

   
    E0taborg = vEtab[0];
   
    double dEtab = log10(vEtab[0]) - log10(vEtab[1]); 
    double evolution = GetEvolution(proc.GetTargetParticle().GetEnergy(), z);
    int i = (int)(( log10(E0taborg) -log10(E1*(1+z)) )/dEtab);
    
    if (i<0) { 
      cout << "WARNING!! GetLambdaTab in " << procName << " : i= " << i << " <0! E1*(1+z) =   " << E1 << "* (1 + " << z <<  ") < " <<  E0taborg << ".. returning lambda[0];"<< endl;
    }
  
    else  if (i>=1001) { 
      cout << "WARNING!! GetLambdaTab in " << procName << " : i>= " << 1001 <<" ! E1*(1+z) =   " << E1 << "* (1 + "<< z << ") .. returning lambda[nentries];" << endl; 
     
    }
    else {
      if (procName == "PP") res=vPPle[i];
      else if (procName == "DPP") res=vDPPle[i];
      else if (procName == "ICS") res=vICSle[i];
      else if (procName == "TPP") res=vTPPle[i];
    }
    
    
    if (evolution!=0)  {
      if (res/evolution < 0) cerr << "ERROR UNPHYSICAL SOLUTION!! CHECK HERE LAMBDA OR EVOLUTION!!" <<endl;
      return res/evolution;
    }
    cerr << "warning!! evolution ==0 " <<endl;
    return 0; 
  };

double PhotonEleCa::ShootPhotonEnergyMC(double z){
    // Routine for the MC sampling of background photon energy
        
    double h = Uniform(0,1);
  
    for(int i=0; i<POINTS_VERY_FEW; i++){
        if( h<BkgArray[i][1] ){
            return BkgArray[i][0]*(1.+z);
            break;
        }
    }
    
    return 0.;
}


vector<double> PhotonEleCa::GetEtarget(Process &proc, Particle &particle) {

  vector<double> Etarget;
  double Etarget_tmp=0;
  double smintmp=0;
  double z_curr=particle.Getz();
  double Energy=particle.GetEnergy();
  int pType = particle.GetType();
  if (pType==22) { 
    proc.SetName("PP");
    proc.SetLimits();
    smintmp=proc.GetMin();
    Etarget_tmp=0;
    while (Etarget_tmp < ElectronMass*ElectronMass/Energy ) {
      Etarget_tmp =  ShootPhotonEnergyMC(z_curr);
    }
    Etarget.push_back(Etarget_tmp);

    proc.SetName("DPP");
    proc.SetLimits();
    smintmp=proc.GetMin();
    Etarget_tmp=0;

    while (Etarget_tmp < smintmp/(4.0*Energy) ) {
      Etarget_tmp =  ShootPhotonEnergyMC(z_curr);
    }
    Etarget.push_back(Etarget_tmp);
  }

  else if (abs(pType)==11){
    proc.SetName("ICS");
    proc.SetLimits();
    smintmp=proc.GetMin();
    Etarget_tmp=0;
    while (Etarget_tmp < smintmp/(4.0*Energy) ) {
      Etarget_tmp =  ShootPhotonEnergyMC(z_curr);
    }
    Etarget.push_back(Etarget_tmp);

    proc.SetName("TPP");
    proc.SetLimits();
    smintmp=proc.GetMin();
    Etarget_tmp=0;

    while (Etarget_tmp<smintmp/(4.0*Energy) ) {
      Etarget_tmp =  ShootPhotonEnergyMC(z_curr);
    }
    Etarget.push_back(Etarget_tmp);
  }//end e/e
  else cerr << "something wrong in particle type ( " << pType << ". Propagation of photons and e+/e- is the only allowed.)" << endl;

  if (Etarget.size()!=2) {cout << "something wrong with the Etarget!! " << endl; exit(0);}
 
  return Etarget;
}

double PhotonEleCa::ExtractPhotonEnergyMC(double z, Process &proc){
  double esoft=0;
  double snew=0;
  double emin = proc.GetMin();
  Particle pi=proc.GetIncidentParticle();
  Particle pb=proc.GetTargetParticle();
  
  double Epi= pi.GetEnergy();
  double m= pi.GetMass();
   while (esoft<emin/(4.0*Epi)) {
   
    double h = Uniform(0,1);

    for(int i=0; i<POINTS_VERY_FEW; i++){
      if( h<BkgArray[i][1] ) {	
         esoft = BkgArray[i][0]*(1.+z);
	break; 
      }
    } 

    snew = 4*Epi*esoft + m*m; 
  }
  pb.SetEnergy(esoft);
  proc.SetTargetParticle(pb);
  proc.SetCMEnergy();
  return esoft;
};


void PhotonEleCa::WriteOutput(Particle & p1, vector<Particle> &ParticleAtGround, bool spectrum) {
  double Ethr=1e16;
  double Bfield=0;
  double E0nucl=p1.GetEnergy(); double z0nucl=p1.Getz(); 

   int NsecG=ParticleAtGround.size();
   vector<double> EGround;
   vector<int> wGround; 
   vector<int> typeGround;
   int cpart = 0;

   for (int ig=0; ig<NsecG;++ig) {
      EGround.push_back(ParticleAtGround.at(ig).GetEnergy());
      typeGround.push_back(ParticleAtGround.at(ig).GetType());
      wGround.push_back(ParticleAtGround.at(ig).GetWeigth());
   
      cpart+=wGround.at(ig);
    }

   if (spectrum) {
      double emin=7.0;
      double dE=(24.0-7.0)/170.0;
      int ipos=0;

      for (int h=0; h<NsecG; ++h) {
	ipos = (int)((log10(EGround.at(h))-emin)/dE);
	if (typeGround.at(h)==22) 
	  fdN[ipos]+= wGround.at(h); 
      }
    }//end opt_spectrum 
    else {
      outfile << Ethr << " " << Bfield/1e-9 << " " << E0nucl << " " << z0nucl << " " << NsecG << "  "; 
      for (int h=0; h<NsecG; ++h) {

	outfile << wGround.at(h) << " " << EGround.at(h) << " " << typeGround.at(h) << "  " ;
      }
    }
    outfile << endl;
};


#endif /* MPC_PHOTON_ELECA_H_ */
