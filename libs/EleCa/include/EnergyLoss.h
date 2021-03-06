namespace crpropa {

  double dTdZ(double z){
    // See High Energy Cosmic Rays, Todor Stanev, Pag. 232 (2009)
    return -1./( (1+z)*H0y*sqrt( pow(1.+z,3.)*OM + OL + (1-OM-OL)*pow(1.+z,2.) ) );
  };
  
  double betaRsh(double z){
    // Energy loss term due to cosmological redshift    
    return H0y*sqrt( pow(1.+z,3.)*OM + OL + (1-OM-OL)*pow(1.+z,2.) );
  };
  
  
  double fLossAdiabatic(double E, double z){
    return  -dTdZ(z) * betaRsh(z)*E;
  };
  
  double AdiabaticELoss(double z0, double z, double E0){
    return E0*(double)(1.+z)/(1.+z0);
  };
  
  
  // ##################################################
  // # Synchrotron rate of Energy loss averaged over angles
  // ##################################################
  // for an ensemble of electrons that are scattered randomly in all directions:
  double MeanRateSynchrotronLoss(double E, double B) {
    double dEdt = 0; 
    if (B>0) dEdt = 3.79e-6*pow(E/1e9*B,2)*1e9/E;
    
    return dEdt;
  };
  
  
  // ####################################
  // # Synchrotron Energy loss
  // ####################################
  
  double ESynchrotronAdiabaticLoss(double z, double E, double B) {
    double dEdt =  MeanRateSynchrotronLoss(E, B);
    
    return E*(-dTdZ(z)*(dEdt+betaRsh(z)) );
  };

    
void InitRK(){
    // Current Runge-Kutta method for solving ODE
    vRKa[0] = 0;
    vRKa[1] = 0;
    vRKa[2] = 1./5.;
    vRKa[3] = 3./10.;
    vRKa[4] = 3./5.;
    vRKa[5] = 1.;
    vRKa[6] = 7./8.;

    for(int i=0; i<RK_ORDER; i++){
        for(int j=0; j<RK_ORDER+1; j++) vRKb[i][j] = 0.;
    }
    
    vRKb[2][1] = 1./5.;
    vRKb[3][1] = 3./40.;
    vRKb[3][2] = 9./40.;
    vRKb[4][1] = 3./10.;
    vRKb[4][2] = -9./10.;
    vRKb[4][3] = 6./5.;
    vRKb[5][1] = -11./54.;
    vRKb[5][2] = 5./2.;
    vRKb[5][3] = -70./27.;
    vRKb[5][4] = 35./27.;
    vRKb[6][1] = 1631./55296.;
    vRKb[6][2] = 175./512.;
    vRKb[6][3] = 575./13824.;
    vRKb[6][4] = 44275./110592.;
    vRKb[6][5] = 253./4096.;

    vRKc[0] = 0.;
    vRKc[1] = 37./378.;
    vRKc[2] = 0.;
    vRKc[3] = 250./621.;
    vRKc[4] = 125./594.;
    vRKc[5] = 0.;
    vRKc[6] = 512./1771.;

    vRKcs[0] = 0.;
    vRKcs[1] = 2825./27648.;
    vRKcs[2] = 0.;
    vRKcs[3] = 18575./48384.;
    vRKcs[4] = 13525./55296.;
    vRKcs[5] = 277./14336.;
    vRKcs[6] = 1./4.;
};


//===================================

double EnergyLoss1D(double Energy, double z0, double zfin, double B){
  
  double zStep = 2.5e-5;

  InitRK();
  double k1, k2, k3, k4, k5, k6;
  
  bool FLAG_PROPAG = 1;

  while(z0>zfin && FLAG_PROPAG){    

    k1 = -zStep * ESynchrotronAdiabaticLoss( z0-zStep, Energy, B );
    k2 = -zStep * ESynchrotronAdiabaticLoss( z0-zStep*vRKa[2], Energy + k1*vRKb[2][1],B );

    k3 = -zStep * ESynchrotronAdiabaticLoss( z0-zStep*vRKa[3], Energy + k1*vRKb[3][1] + k2*vRKb[3][2], B );
    k4 = -zStep * ESynchrotronAdiabaticLoss( z0-zStep*vRKa[4], Energy + k1*vRKb[4][1] + k2*vRKb[4][2] + k3*vRKb[4][3], B );
    k5 = -zStep * ESynchrotronAdiabaticLoss( z0-zStep*vRKa[5], Energy + k1*vRKb[5][1] + k2*vRKb[5][2] + k3*vRKb[5][3] + k4*vRKb[5][4] , B); 
    k6 = -zStep * ESynchrotronAdiabaticLoss( z0-zStep*vRKa[6], Energy + k1*vRKb[6][2] + k2*vRKb[6][2] + k3*vRKb[6][3] + k4*vRKb[6][4] + k5*vRKb[6][5], B );
    
    Energy = Energy + ( k1*vRKc[1] + k2*vRKc[2] + k3*vRKc[3] + k4*vRKc[4] + k5*vRKc[5] + k6*vRKc[6] );

    z0 -= zStep;  

    if(fabs(z0)<1e-8 || z0<0) { 
      z0 = 0.; 
      FLAG_PROPAG =0;}
    if (fabs(z0-zfin)<1e-8 || z0<zfin) z0 =zfin;

    if( Energy<1e9) {
      FLAG_PROPAG = 0;
    }

    if( isnan(Energy) ) return 0;
  }

    return Energy;
};
    

double dSigmadE_ICS(double Ee, double Eer, double s, double theta) {
    /*!
    Differential cross-section for inverse Compton scattering. from lee, eq. 23
    */ 

    double beta = (s-ElectronMass*ElectronMass)/(s+ElectronMass*ElectronMass);
  
    if (Eer/Ee < (1-beta)/(1+beta) || Eer/Ee>1) { 
        cerr << "ERROR, Energy outside limits for ICS [Lee96]! " <<endl; 
        return 0.;   
    }else{
        double q = ((1-beta)/beta)*(1-Ee/Eer);
        double A = Eer/Ee + Ee/Eer;
        double k = (3.0/8.0)*(SigmaThompson*ElectronMass*ElectronMass)/(s*Ee); 
	double dsigmadE = k*((1+beta)/beta)*(A + 2*q + q*q);
        
        return dsigmadE;
    }
};


double dSigmadE_PP(double Ee, double E0, double eps, double theta) {
    /*!
        Differential cross-section for pair production.
    */ 
    double s = ElectronMass*ElectronMass + 2*eps*E0*(1-cos(theta));
    double beta = sqrt(1-4*ElectronMass*ElectronMass/s);

    if (Ee/E0 <= 0.5*(1-beta) || Ee/E0>= 0.5*(1+beta)) { 
        cerr << "ERROR, Energy outside limits for PP [Lee96]! " <<endl; 
        return 0.;   
    }else{
        double q = E0-Ee;
        double k = (3.0/4.0)*(SigmaThompson*ElectronMass*ElectronMass)/(s*E0);
        double A = Ee/q + q/Ee;
        double B = E0*(1-beta*beta)*(1./Ee + 1./q);
        double C = -( (1-beta*beta) * (1-beta*beta) * E0*E0/4.0 ) * pow(1./Ee + 1./q,2);
        
        double dsigmadE = k*(A+B+C);
    
        return dsigmadE;
    }
};

double ExtractPPSecondariesEnergy(Particle pi, Particle pt) {
    /*!
        Input: incident gamma Energy E0, background photon energy eps, 
        incidence angle theta.
        Returns the energy of the produced e+ (e-)
    */
    double E0 = pi.GetEnergy();
    double eps = pt.GetEnergy();
    double beta = pi.GetBeta();
    double theta = cPI;
    double s2 = ElectronMass*ElectronMass + 2*eps*E0*(1-(beta)*cos(theta));

     bool failed=1;
    
    for(int i=0; i<MC_SAMPLING; i++){
        for(int j=0; j<3; j++) MC_Sampling_Hist[i][j] = 0.;
    }

    double f = pow( (double)(1+beta)/(1-beta) ,(double)1./(double)MC_SAMPLING);
    int cnt = 0;
    double NormFactor = 0;
    
   

    for (double Ee=f*0.5*(1-beta)*E0; Ee<0.5*(1+beta)*E0; Ee*=f){
        MC_Sampling_Hist[cnt][0] = Ee;
        MC_Sampling_Hist[cnt][1] = dSigmadE_PP(Ee, E0, eps, theta);

        NormFactor += MC_Sampling_Hist[cnt][1];
        MC_Sampling_Hist[cnt][2] = NormFactor;
            
        if(MC_Sampling_Hist[cnt][1] > 0.){
            cnt++;
        } else {
            break;
        }
    }
    
    NormFactor = (double)1./(double)NormFactor;
    
    for(int i=0; i<cnt; i++) MC_Sampling_Hist[i][2] *= NormFactor;
    
    double rnd; 
    double Ee = 0;
    int k=0;
   while (failed) { 
     k++;

    rnd = Uniform(0,1);
    Ee = 0;
    for(int i=0; i<cnt-1; i++){
        if( MC_Sampling_Hist[i][2]<= rnd <= MC_Sampling_Hist[i+1][2] ){ 
            Ee = MC_Sampling_Hist[i][0];
            failed=0;
            break;
        } 
     	if (failed) 
	  cout << "cnt: " << cnt << " failed " << k << endl;
     }

    } //end while
   
   if( Uniform(0,1)<0.5 ) 
        return Ee;
    else
        return E0-Ee;
};

double ExtractPPSecondariesEnergy(Process proc) {
    /*!
        Input: incident gamma Energy E0, background photon energy eps, 
        incidence angle theta.
        Returns the energy of the produced e+ (e-)
    */

    double E0 = proc.GetIncidentParticle().GetEnergy();
    double s = proc.GetCMEnergy();
    double eps = proc.GetTargetParticle().GetEnergy();
    double theta = cPI;

    double beta = sqrt(1.-4.0*ElectronMass*ElectronMass/s);
    double s2 = ElectronMass*ElectronMass + 2*eps*E0*(1-(beta)*cos(theta));

    bool failed=1;
  
    for(int i=0; i<MC_SAMPLING; i++){
        for(int j=0; j<3; j++) MC_Sampling_Hist[i][j] = 0.;
    }

    double f = pow( (double)(1+beta)/(1-beta) ,(double)1./(double)MC_SAMPLING);
    int cnt = 0;
    double NormFactor = 0;
   
    for (double Ee=f*0.5*(1-beta)*E0; Ee<0.5*(1+beta)*E0; Ee*=f){
        MC_Sampling_Hist[cnt][0] = Ee;
        MC_Sampling_Hist[cnt][1] = dSigmadE_PP(Ee, E0, eps, theta);
        NormFactor += MC_Sampling_Hist[cnt][1];
        MC_Sampling_Hist[cnt][2] = NormFactor;
            
        if(MC_Sampling_Hist[cnt][1] > 0.){
            cnt++;
        } else {
            break;
        }
    }
    
    NormFactor = (double)1./(double)NormFactor;
    
    for(int i=0; i<cnt; i++) MC_Sampling_Hist[i][2] *= NormFactor;
    
    double rnd; 
    double Ee = 0;
    int k=0;

   while (failed) { 
    rnd = Uniform(0.,1.0);
    Ee = 0;
    k++;
    double min=1e6;
    double max=-1;

    for(int i=0; i<cnt-1; i++){

	    if (MC_Sampling_Hist[i][2]<min) min=MC_Sampling_Hist[i][2];
	    if (MC_Sampling_Hist[i+1][2]>max) max=MC_Sampling_Hist[i+1][2];

           if( MC_Sampling_Hist[i][2]<= rnd <= MC_Sampling_Hist[i+1][2] ){ 
            Ee = MC_Sampling_Hist[i][0];
            failed=0;
            break;
         }
    }
    if (failed) {
	  cout << "failed in extractPP " << Ee << " " << beta  << " * s: " <<  s << " E: " <<  E0 << " eps : " << eps << " me: " << ElectronMass*ElectronMass/E0 << "  ) " << " cnt : " << cnt << endl;
    cout << " Limits  " << proc.GetMin() << endl;
    if (cnt==0)   exit(0);
    }

  } //end while
   
   if( Uniform(0,1.0)<0.5 ) 
        return Ee;
    else
        return E0-Ee;
};


double ExtractICSSecondariesEnergy(Particle pi, Particle pt) {
    /*!
        Input: incident electron energy Ee, background photon energy eps, 
        incidence angle theta.
        Returns the energy of the recoiled e+ (e-)
    */
  if (abs(pi.GetType())!=11) {
    cerr << "something wrong in type ExtractICSEnergy " << endl;
    return 0.; 
  }

  double Ee = pi.GetEnergy();
  double eps=pt.GetEnergy();
  double theta=cPI;
  double s  = 2*Ee*eps*(1-pi.GetBeta()*cos(cPI)) + pi.GetMass()*pi.GetMass();
  double beta = (s-ElectronMass*ElectronMass)/(s+ElectronMass*ElectronMass);
  bool failed=1;       


    for(int i=0; i<MC_SAMPLING; i++){
        for(int j=0; j<3; j++) MC_Sampling_Hist[i][j] = 0.;
    }

    double f = pow( (double)(1+beta)/(1-beta) ,(double)1./MC_SAMPLING);
    int cnt = 0;
    double NormFactor = 0;


    for (double Eer=f*((1-beta)/(1+beta))*Ee; Eer<=Ee; Eer*=f){
        MC_Sampling_Hist[cnt][0] = Eer;
        MC_Sampling_Hist[cnt][1] = dSigmadE_ICS(Ee, Eer, s, theta);

        NormFactor += MC_Sampling_Hist[cnt][1];
        MC_Sampling_Hist[cnt][2] = NormFactor;
            
        if(MC_Sampling_Hist[cnt][1] > 0.){
            cnt++;
        }else{
            break;
        }
    }
    
    NormFactor = (double)1./(double)NormFactor;

    for(int i=0; i<cnt; i++) MC_Sampling_Hist[i][2] *= NormFactor;
    
    double rnd = 0; 
    double Eer = 0;

   while (failed) {
   rnd = Uniform(0,1);
   Eer=0;
    for(int i=0; i<cnt-1; i++){
        if( MC_Sampling_Hist[i][2]<= rnd <= MC_Sampling_Hist[i+1][2] ){ 
            Eer = MC_Sampling_Hist[i][0];
            failed=0;
            break;
        }
    }

  }
     return Eer;
};



double ExtractICSSecondariesEnergy(Process proc) {
 /*!
   Input: incident electron energy Ee, background photon energy eps, 
   incidence angle theta.
   Returns the energy of the recoiled e+ (e-)
 */
  double Ee = proc.GetIncidentParticle().GetEnergy();
  double s = proc.GetCMEnergy();
  double theta = proc.GetInteractionAngle();
  double beta = (s-ElectronMass*ElectronMass)/(s+ElectronMass*ElectronMass);
  bool failed=1;       

    // reInitialization to zero.. 
    for(int i=0; i<MC_SAMPLING; i++){
        for(int j=0; j<3; j++) MC_Sampling_Hist[i][j] = 0.;
    }

    double f = pow( (double)(1+beta)/(1-beta) ,(double)1./MC_SAMPLING);
    int cnt = 0;
    double NormFactor = 0;


    for (double Eer=f*((1-beta)/(1+beta))*Ee; Eer<=Ee; Eer*=f){
        MC_Sampling_Hist[cnt][0] = Eer;
        MC_Sampling_Hist[cnt][1] = dSigmadE_ICS(Ee, Eer, s, theta);

        NormFactor += MC_Sampling_Hist[cnt][1];
        MC_Sampling_Hist[cnt][2] = NormFactor;
            
        if(MC_Sampling_Hist[cnt][1] > 0.){
            cnt++;
        }else{
            break;
        }
    }
    
    NormFactor = (double)1./(double)NormFactor;

    for(int i=0; i<cnt; i++) MC_Sampling_Hist[i][2] *= NormFactor;
    
    double rnd = 0; 
    double Eer = 0;

    while (failed) {
      rnd = Uniform(0,1.0);
      Eer=0;
      for(int i=0; i<cnt-1; i++){
        if( MC_Sampling_Hist[i][2]<= rnd <= MC_Sampling_Hist[i+1][2] ){ 
	  Eer = MC_Sampling_Hist[i][0];
	  failed=0;
	  break;
        }
    }

}
     return Eer;
};


double ExtractTPPSecondariesEnergy(Particle pi, Particle pt) {
  /* approximation based on A. Mastichiadis et al., 
    Astroph. Journ. 300:178-189 (1986), eq. 30. 
    This approx is valid only for   alpha >=100 
    where alpha = p0*eps*costheta - E0*eps; 
    for our purposes, me << E0 --> p0~ E0 --> 
    alpha = E0*eps*(costheta - 1) >= 100; 
  */

  double E0 = pi.GetEnergy();
  double eps = pt.GetEnergy();
  double s = 2*E0*eps*(1-pi.GetBeta()*cos(cPI)) + pi.GetMass()*pi.GetMass();
  double Epp= 5.7e-1 * pow(eps,-0.56)*pow(E0,0.44);
  double Epp2= E0*(1 - 1.768*pow(s/ElectronMass/ElectronMass,-3.0/4.0))/2.0;
  //return the energy of each e+/e- in the pair. 
  return Epp;

};



double ExtractTPPSecondariesEnergy(Process proc) {
  /* approximation based on A. Mastichiadis et al., 
     Astroph. Journ. 300:178-189 (1986), eq. 30. 
     This approx is valid only for   alpha >=100 
     where alpha = p0*eps*costheta - E0*eps; 
     for our purposes, me << E0 --> p0~ E0 --> 
     alpha = E0*eps*(costheta - 1) >= 100; 
  */

  double E0 = proc.GetIncidentParticle().GetEnergy();
  double eps = proc.GetTargetParticle().GetEnergy();
  double Epp= 5.7e-1 * pow(eps,-0.56)*pow(E0,0.44);
  double s = proc.GetCMEnergy();
  double Epp2= E0*(1 - 1.768*pow(s/ElectronMass/ElectronMass,-3.0/4.0))/2.0;
  return Epp;
};


double ExtractDPPSecondariesEnergy(double E0)  {
  /*
    we use the same assumption of lee (i.e., all the energy goes equaly shared between only 1 couple of e+e-. 
    In DPPpaper has been shown that this approximation is valid within -1.5%
  */
  if (E0==0) cout << "error in extracting DPP: can not be =0 " <<endl;
  return (double)E0/2.0;
}    

}
