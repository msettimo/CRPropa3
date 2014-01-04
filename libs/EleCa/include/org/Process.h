class Process {

 public: 

  string fname;
  double flambda;
  double fsigma; 
  double fsmin; 
  double fsmax;
  double fCMEnergy;
  double fInteractionAngle;  // in rad
  double fInteractionDistance;
  double fEnergyLoss;
  Particle fPi; //incident particle
  Particle fPt; //target particle
  Ambient fambient;
  string fback;
  double fbackdensity;
  double fBfield;

  /*  to move somewhere else!! the class process is intended to describe only a process
      vector<string> fallnames;
      vector<double> fallLambdas;
      vector<double> fallSigmas;
      void SetName(vector<string> names) {fallnames.Clear(); fallnames = names;}
      vector<string> GetNames() {return fallnames;}
      void SetLambda(vector<double> lle) {fallLambdas.Clear(); fallLambdas = lle;} 
      void SetSigma(vector<double> ss) {fallSigmas.Clear(); fallSigmas = ss;}
  */

  Process();
  Process(Process&);
  Process(Particle&, Particle&);
  Process(Particle&, Particle&, string);

  ~Process();

  void SetName(string nm) {fname = nm;}
  string GetName() {return fname;}

  void SetInteractionAngle(double a) {fInteractionAngle = a;} 
  double GetInteractionAngle() {return fInteractionAngle;} 

  void SetInteractionDistance(double d1) {fInteractionDistance = d1;} 
  double GetInteractionDistance() {return fInteractionDistance;} 

 void SetEnergyLoss(double eta1) {fEnergyLoss = eta1;} 
  double GetEnergyLoss() {return fEnergyLoss;} 


  void SetLambda(double le) {flambda = le;}
  //  void SetLambda(Particle& p1, Particle& pb) {flambda = ComputeLambda(p1, pb);}
  double GetLambda() {return flambda;} 

  void SetSigma(double s1) {fsigma = s1;}
  //  void SetSigma(Particle& p1, Particle& pb) {fsigma = ComputeSigma(p1, pb);}
  double GetSigma() {return fsigma;}

  void SetLimits(double smin, double smax) {fsmin = smin; fsmax = smax;}
  void SetLimits(Particle& p1, string nameproc);
  void SetLimits() {SetLimits(fPi,fname);}
  
  // void SetThreshold(double smin) {fsmin = smin;}
  // void SetThreshold(Particle p1, string nameproc, string BackRad);
  //  double GetThreshold() {return fsmin;}
 
  // void SetMax(Particle p1, string nameproc, string BackRad);
  void SetMax(double smax) {fsmax = smax;}
  void SetMin(double smin) {fsmin = smin;}
  double GetMin() {return fsmin;}
  double GetMax() {return fsmax;}

  void SetB(double b) {fBfield=b;}
  double GetB() {return fBfield;}

  void SetCMEnergy(double s) {fCMEnergy = s;}
  void SetCMEnergy(Particle p1, Particle pb) {
    fCMEnergy = 2*p1.GetEnergy()*pb.GetEnergy()*(1-p1.GetBeta()*cos(fInteractionAngle)) + p1.GetMass()*p1.GetMass() + pb.GetMass()*pb.GetMass(); }
  void SetCMEnergy() {
//   cout << " debug: " << fPi.GetEnergy() << "  " << fPt.GetEnergy() << "  " << fPi.GetBeta() << " " << fInteractionAngle << " " << fPi.GetMass() << "  "  << fPt.GetMass() << endl; 
   fCMEnergy = 2*fPi.GetEnergy()*fPt.GetEnergy()*(1-fPi.GetBeta()*cos(fInteractionAngle)) + fPi.GetMass()*fPi.GetMass() + fPt.GetMass()*fPt.GetMass();
//   cout << fCMEnergy << "  ["  << fsmin <<" , " << fsmax << "] vs " <<  1-4.0*ElectronMass*ElectronMass/fCMEnergy << endl;

   }

  double GetCMEnergy() {return fCMEnergy;}

  //  double ComputeLambda(Particle& p1, Particle& p2);
  // double ComputeSigma(Particle& p1, Particle& p2);

  void SetIncidentParticle(Particle& p1) {fPi = p1; SetLimits();}
  void SetTargetParticle(Particle& p1) {fPt = p1;}

  Particle GetIncidentParticle() {return fPi;}
  Particle GetTargetParticle() {return fPt;}

  string GetBackground() {return fback;}
  void SetBackground(string BackRad);

  void SetBackgroundDensity(double n_eps) {fbackdensity = n_eps;}
  double GetBackgroundDensity() {return fbackdensity;}

  Ambient GetAmbient() {return fambient;}
  void SetAmbient(Ambient ab) {fambient=ab;}

 private:
   

};

//-------------- explicit definitions
Process::Process() {
  fname="";
  SetLimits(0.0,1.0e23);
  fsigma=0;
  flambda =0;
  fCMEnergy =0;
  fInteractionAngle=cPI;
  fInteractionDistance=0;
  fEnergyLoss =0;
  fback="ALL";
  fbackdensity=0;
}


Process::Process(Particle& p1, Particle& p2) {
  fPi = p1;
  fPt = p2;
  if (p1.GetType()==0)
    fname = "PP";
  else if (abs(p1.GetType())==1) { 
    cerr << "NB: by default process set to ICS" << endl;
    fname="ICS";
  }
  else fname="NONE";
  SetCMEnergy(p1, p2);
  fsigma=0;
  flambda=0;
  fInteractionAngle= cPI;
  fInteractionDistance = 0; 
  fEnergyLoss= 0;
  fback="ALL";
  SetLimits(p1, fname);
  fbackdensity=0;
}


Process::Process(Particle& p1, Particle& p2, string name) {
  fname = name;
  SetCMEnergy(p1, p2);
  fsigma=0;
  flambda=0;
  //ComputeSigma(p1, p2, name);
  // ComputeLambda(p1, p2, name);
  fInteractionAngle= cPI; // impostare che venga gia calcolato... ma interaction angle e' estratto in funzone del processo.. mmh.. vedere casomai sei ricorsiva.. 
  fInteractionDistance = 0; //ExtractInteractionDistance(flambda);
  fEnergyLoss= 0;
  fPi = p1;
  fPt = p2;
  fback="ALL";
  SetLimits(p1, fname);
  fbackdensity=0;
}

Process::Process(Process& proc2) {
  fname = proc2.GetName();
  SetLimits(proc2.GetMin(), proc2.GetMax());
  flambda = proc2.GetLambda();
  fsigma = proc2.GetSigma();
  fCMEnergy = proc2.GetCMEnergy();
  fInteractionAngle = proc2.GetInteractionAngle();
  fInteractionDistance = proc2.GetInteractionDistance();
  fEnergyLoss = proc2.GetEnergyLoss();
  fPi = proc2.GetIncidentParticle();
  fPt = proc2.GetTargetParticle();
  fback = proc2.GetBackground();
  fbackdensity=0;
}

Process::~Process() {
}
 

//-----------

void Process::SetBackground(string BackRad) {

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

void Process::SetLimits(Particle& p1, string nameproc) {

   if (nameproc == "PP") {
     if (abs(p1.GetType())==1) 
       cout << "\nERROR!! wrong particle or process!! " << nameproc << p1.GetType() << "\n" << endl;
     fsmin = 4*ElectronMass*ElectronMass; //min per theta = 0;
     fsmax = 4*p1.GetEnergy()*feps_sup*2.0;//(1.0-cos(fInteractionAngle)); //max per theta=3.14
   }
   if (nameproc == "DPP") {
     if (abs(p1.GetType())==1) 
       cout << "\nERROR!! wrong particle or process!! " << nameproc << p1.GetType() << "\n" << endl;
     fsmin = 16*ElectronMass*ElectronMass;
     fsmax = 4*p1.GetEnergy()*feps_sup*2.0;//(1.0-cos(fInteractionAngle));
   }
    if (nameproc == "ICS") {
      if (p1.GetType()==0 || p1.GetType()==9) 
       cout << "\nERROR!! wrong particle or process!! " << nameproc << p1.GetType() << "\n" << endl;
      fsmin = p1.GetMass()*p1.GetMass() + 2*p1.GetEnergy()*feps_inf*(1 - p1.GetBeta()); //min for theta =0;
      fsmax = p1.GetMass()*p1.GetMass() + 2*p1.GetEnergy()*feps_sup*(1 + p1.GetBeta()) ;
    }
    if (nameproc == "TPP") {
      if (p1.GetType()==0 || p1.GetType()==9)
       cout << "\nERROR!! wrong particle or process!! " << nameproc << p1.GetType() << "\n" << endl;
      fsmin = 10*ElectronMass*ElectronMass;
      fsmax = 2*p1.GetEnergy()*feps_sup*(1 + p1.GetBeta()) + p1.GetMass()*p1.GetMass() ;
    }
  }
