vector<double> GetEtarget(Process &proc, Particle &particle) {
  
  if (particle.GetType()!=22 || abs(particle.GetType())!=11) {
    cerr << "something wrong with particle type" <<endl;
    return;
  }

  vector<double> Etarget;
  double Etarget_tmp=0;
  double smintmp=0;
  double z_curr=particle.Getz();
  double Energy=particle.GetEnergy();
  int pType = particle.GetType();
  if (pType==0) { 
    proc.SetName("PP");
    proc.SetLimits();
    smintmp=proc.GetMin();

    while (Etarget_tmp<smintmp/(4.0*Energy) ) {
      Etarget_tmp =  ShootPhotonEnergyMC(z_curr);
    }

    Etarget.push_back(Etarget_tmp);
    proc.SetName("DPP");
    proc.SetLimits();
    smintmp=proc.GetMin();

    while (Etarget_tmp<smintmp/(4.0*Energy) ) {
      Etarget_tmp =  ShootPhotonEnergyMC(z_curr);
    }
    Etarget.push_back(Etarget_tmp);
  }//end if pType

  else if (abs(pType)==1){
    proc.SetName("ICS");
    proc.SetLimits();
    smintmp=proc.GetMin();
    while (Etarget_tmp<smintmp/(4.0*Energy) ) {
      Etarget_tmp =  ShootPhotonEnergyMC(z_curr);
    }
    Etarget.push_back(Etarget_tmp);

    proc.SetName("TPP");
    proc.SetLimits();
    smintmp=proc.GetMin();
    while (Etarget_tmp<smintmp/(4.0*Energy) ) {
      Etarget_tmp =  ShootPhotonEnergyMC(z_curr);
    }
    Etarget.push_back(Etarget_tmp);
  }//end e/e
  else cerr << "something wrong in particle type ( " << pType << ". Propagation of photons and e+/e- is the only allowed.)" << endl;

  if (Etarget.size()!=2) {cout << "something wrong with the Etarget!! " << endl; exit(0);}
  return Etarget;
}

//======================

double ExtractMinDist(Process &proc, Particle &pt, int type, double R, double R2, vector<double> Etarget) {

  double min_dist1=0; double min_dist2=0;
  Process proc1(proc);
  Process proc2(proc);
  double tmp_lambda1=0; double tmp_lambda2=0;

   if (type==0) {
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
   else if (abs(type)==1) {
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
}


double ExtractMinDist(Process &proc, Particle &pt, int type, double R, double R2, vector<double> Etarget, double corrB) {
  if (corrB==0) corrB=1;
  double dist=ExtractMinDist(proc, pt, type, R, R2, Etarget);
  return dist/corrB;
}

