class Particle { 

 public: 
  int fpType;  // possible values: 22 = photons, +/- 11: electrons; 
  double fz;
  double ftheta;
  double fbeta;
  double fphi;
  double fmass;
  double fEnergy;
  vector<double> fposition;
  vector<double> fvelocity;
  double fvx;
  double fvy;
  double fvz;

  double fposx;
  double fposy;
  double fposz;

  double fdistance;
  int fId;
  int fParentId;

  Particle();
  Particle(vector<double>); // 8 entries

  /* Particle(id, parentid, type, z, E, theta, phi); */
  Particle(int, int, int, double, double, double, double); 

 /* Particle(type, z, E,theta); */
  Particle(int, double, double,  double); 
  ~Particle();

  void SetParticle(vector<double>);
  void SetType(int type) {fpType = type; SetMassAndBeta(type); }
  void SetId(int id) {fId = id;}
  void SetParentId(int id) {fParentId = id;}

  void SetEnergy(double e1) {fEnergy = e1;}
  void Setz(double z1) {fz = z1;}
  void SetTheta(double t1) {ftheta = t1;}
  void SetMass(double m) {fmass = m;}
  void SetBeta(double b) {fbeta = b;}
  void SetPhi(double phi) {fphi = phi;}
  void SetPhi() {fphi = 360*Uniform(0.0,1.0);}


  void SetVelocity(vector<double> velocity) {fvelocity=velocity;}
  void SetVelocity();

  int GetType() {return fpType;}
  int GetId() {return fId;}
  int GetParentId() {return fParentId;}
  double GetEnergy() {return fEnergy; }
  double Getz() {return fz;}
  double GetTheta() {return ftheta;}
  double GetPhi() {return fphi;}
  double GetMass() {return fmass;}
  double GetBeta() {return fbeta;}

 private:

  void SetMassAndBeta(int);
  
};


Particle::Particle () {
  fpType = 0;
  fId = -9;
  fParentId = -1;
  fEnergy = 0;
  fposx=fposy=fposz=0;
  fvx=fvy=fvz=0;
  ftheta =0;
  fphi=2*cPI*Uniform(0.0,1.0); //2*cPI*Uniform(0.0,1.0)
  fz = 0;
  fmass=0;
  fbeta=1;
}



Particle::Particle (vector<double> tmp) {
  fpType = (int)tmp.at(0);
  fEnergy = tmp.at(1);
  //  fposx = tmp.at(2); 
  // fposy= tmp.at(3); 
  // fposz=tmp.at(4);
  SetMassAndBeta((int)tmp.at(0));
  ftheta=tmp.at(2);
  fphi=tmp.at(3);
  fz = tmp.at(4);
  SetVelocity();
}

Particle::Particle(int p1, double z1, double E1, double t1) {
  fpType = p1;
  fEnergy = E1;
  fz = z1;
  ftheta =t1;
  fposx=0;
  fposy=0;
  fposz=0;

  fphi=2*cPI*Uniform(0.0,1.0);
  if (fpType == 0 || fpType == 9) { fmass=0; fbeta = 1;}
  else { 
    fmass = ElectronMass;
    if (fEnergy == 0) { 
      cout << "ERROR! Set Energy before Particle Type and Mass/Beta values" <<endl;
      fbeta=1;}
    else 
      fbeta = sqrt(1 - (fmass*fmass)/(fEnergy*fEnergy));
    if (fbeta>1)  cout << "ERROR -- beta can not be larger than 1!!!" << endl;
  }
 SetVelocity();
}

Particle::~Particle() {
  /*  delete fpType; 
  delete fEnergy;
  delete fz;
  delete ftheta; 
  delete fmass;
  delete fbeta;
  delete fphi; */
}


void Particle::SetVelocity() { 
  //in parsec per year
  // particle velocity (with random direction on a sphere)
   fvx = fbeta*c_speed_pc_yr*sin(ftheta)*cos(fphi);
   fvy  = fbeta*c_speed_pc_yr*sin(ftheta)*sin(fphi);
   fvz = fbeta*c_speed_pc_yr*cos(ftheta);  
}

void Particle::SetMassAndBeta(int ptype) {
   if (ptype == 0 || ptype == 9) { 
     SetMass(0); 
     SetBeta(1);}
  else { 
    fmass=ElectronMass;
  if (fEnergy == 0) { 
      cout << "ERROR! Set Energy before Particle Type and Mass/Beta values" <<endl;
      fbeta=1;}
    else 
      fbeta = (double)sqrt(1 - (fmass*fmass)/(fEnergy*fEnergy));
  }
} 



void Particle::SetParticle(vector<double> tmp) {
  fId = (int)tmp.at(0);
  fParentId = (int)tmp.at(1);
  fpType = (int)tmp.at(2);
  fz = tmp.at(3);
  fEnergy = tmp.at(4);
  ftheta=tmp.at(5);
  fphi=tmp.at(6);
  //  fposx = tmp.at(2); 
  // fposy= tmp.at(3); 
  // fposz=tmp.at(4);
  SetMassAndBeta(fpType);
  SetVelocity();
}
