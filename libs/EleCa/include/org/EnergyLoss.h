double MaxLossSigmaPP(Process proc) {
  // to normalize dSigma/dEe I did a check of dsigma/dEe for a scan of Eg and eps extracted in background ALL
  // for each value of Egamma I did 100 measurement of dSigma/dEe and I took the maximum.
  // the log10(max/sigma(s)) is parameterized as a function of log10(s1) with a pol3 function. 
  //see plot test1.root and newdSigmadEePP.root
  double s = proc.GetCMEnergy();
  double  sigma = ComputeSigma(proc);
  double par[4] = {4.8962,-1.4028,0.01437,-0.000192};
  double x = log10(s);
  double log10maxrelsigma = par[0]+par[1]*x + par[2]*x*x + par[3]*x*x*x;
  return pow(10,log10maxrelsigma)*sigma;
}


// ####################################################
// # Metric element
// ####################################################
double dTdZ(double z){
    // Metric element
    
    // See High Energy Cosmic Rays, Todor Stanev, Pag. 232 (2009)
    return -1./( (1+z)*H0y*sqrt( pow(1.+z,3.)*OM + OL + (1-OM-OL)*pow(1.+z,2.) ) );
}


// ####################################################
// # Adiabatic expansion
// ####################################################

double betaRsh(double z){
    // Energy loss term due to cosmological redshift    
    return H0y*sqrt( pow(1.+z,3.)*OM + OL + (1-OM-OL)*pow(1.+z,2.) );
}


// ###################################################
// # Adiabatic losses
// ###################################################
   /* Energy loss Model */
 /* 
     This gives a rate: 1/tau
     To obtain attenuation distance, just invert and
     multiply by (C_speed/H0/1000.)
    */

double fLossAdiabatic(double E, double z){
  return  -dTdZ(z) * betaRsh(z)*E;
}


    //AdiabaticELoss modified the order of z0 and z
double AdiabaticELoss(double z0, double z, double E0){
  // z is the current redshift
  // z0 is the last redshift where particle has been placed
  // E0 is the energy of the particle at z0
  
  // This analytic function can be found by solving
  // the energy loss ODE
  
  return E0*(double)(1.+z)/(1.+z0);
    }

// ##################################################
// # Synchrotron rate of Energy loss (Stanev, book)
// ##################################################
// alpha is the angle of the electron with respect to the direction of B (i.e. alpha is the pitch angle)

double RateSynchrotronLoss(double E, double B, double alpha) {
  double gamma = E/ElectronMass;
  double beta = 1- 1/pow(gamma,2);
  double sigma = SigmaThompson*1e-4;//conv in 1/m^2
  double dEdt =0;
  if (B>0) dEdt= 2*sigma*C_speed*gamma*gamma*B*B/(8.0*cPI)*beta*beta*pow(sin(alpha),2);
  return dEdt;
}



// ##################################################
// # Synchrotron rate of Energy loss averaged over angles
// ##################################################
// for an ensemble of electrons that are scattered randomly in all directions: 

double MeanRateSynchrotronLoss(double E, double B) {
  double dEdt = 0; 
 if (B>0) dEdt = 3.79e-6*pow(E/1e9*B,2)*1e9/E;

  return dEdt;
}


// ####################################
// # Synchrotron Energy loss
// ####################################

double ESynchrotronAdiabaticLoss(double z, double E, double B) {
  double dEdt =  MeanRateSynchrotronLoss(E, B);

  return E*(-dTdZ(z)*(dEdt+betaRsh(z)) );
}



// ####################################################
// #from hermes
// ####################################################
    
void InitRK(){
    // Current Runge-Kutta method for solving ODE

    // Cash-Karp parameters for embedded RK(4,5)
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
}


//===================================

double EnergyLoss1D(double Energy, double z0, double zfin, double B, double zStep){
  
  // double zStep = 2.5e-5;//5*(double)z0/N_STEPS_Z_PROPAGATION;
  //  if (zStep>z0) zStep=z0/10.;
  InitRK();
  double k1, k2, k3, k4, k5, k6;
  
  bool FLAG_PROPAG = 1;

  while(z0>zfin && FLAG_PROPAG){    
    
    // RK(4,5) with Cash-Karp parameters. See NumRec chap 16 pag.718
    // Evolution of E only for protons
    //  cout <<  " k1 : " << Energy <<   "  " << z0 - zStep << "  " << B  << endl;
    k1 = -zStep * ESynchrotronAdiabaticLoss( z0-zStep, Energy, B );
    k2 = -zStep * ESynchrotronAdiabaticLoss( z0-zStep*vRKa[2], Energy + k1*vRKb[2][1],B );
    // cout << "vRKA: " << vRKa[2] << "  " << k1 << "  " <<  vRKb[2][1] << endl;
    k3 = -zStep * ESynchrotronAdiabaticLoss( z0-zStep*vRKa[3], Energy + k1*vRKb[3][1] + k2*vRKb[3][2], B );
    k4 = -zStep * ESynchrotronAdiabaticLoss( z0-zStep*vRKa[4], Energy + k1*vRKb[4][1] + k2*vRKb[4][2] + k3*vRKb[4][3], B );
    k5 = -zStep * ESynchrotronAdiabaticLoss( z0-zStep*vRKa[5], Energy + k1*vRKb[5][1] + k2*vRKb[5][2] + k3*vRKb[5][3] + k4*vRKb[5][4] , B); 
    k6 = -zStep * ESynchrotronAdiabaticLoss( z0-zStep*vRKa[6], Energy + k1*vRKb[6][2] + k2*vRKb[6][2] + k3*vRKb[6][3] + k4*vRKb[6][4] + k5*vRKb[6][5], B );
    
    Energy = Energy + ( k1*vRKc[1] + k2*vRKc[2] + k3*vRKc[3] + k4*vRKc[4] + k5*vRKc[5] + k6*vRKc[6] );
    //   cout <<" in RK: " <<  k1 << "  " << k2 << "   " << k3 << "  " << k4 << "  " << k5 << "  " << k6 << " " <<   Energy << endl; 
 
    
    z0 -= zStep;  

    if(fabs(z0)<1e-8 || z0<0) { 
      z0 = 0.; 
      FLAG_PROPAG =0;}
    if (fabs(z0-zfin)<1e-8 || z0<zfin) z0 =zfin;

    if( Energy<1e9) {
      FLAG_PROPAG = 0;
    }
    //  cout << " in runge kutta " << Energy << endl;                       
    if( isnan(Energy) ) return 0;
  }

    return Energy;
 }    
    


//--------

double dSigmadE_ICS(double Ee, double Eer, double s, double theta) {
    /*!
        Differential cross-section for inverse Compton scattering. from lee, eq. 23
    */ 

    double beta = (s-ElectronMass*ElectronMass)/(s+ElectronMass*ElectronMass);
        
    // for kinematic reasons..
    if (Eer/Ee < (1-beta)/(1+beta) || Eer/Ee>1) { 
        cerr << "ERROR, Energy outside limits for ICS [Lee96]! " <<endl; 
        return 0.;   
    }else{
        double q = ((1-beta)/beta)*(1-Ee/Eer);
        double A = Eer/Ee + Ee/Eer;
    
        // NOTE: see eq 9 in documentation, a c^4 factor is missing..
        double k = (3.0/8.0)*(SigmaThompson*ElectronMass*ElectronMass)/(s*Ee);
            
        double dsigmadE = k*((1+beta)/beta)*(A + 2*q + q*q);
        
        return dsigmadE;
    }
}


double dSigmadE_PP(double Ee, double E0, double eps, double theta) {
    /*!
        Differential cross-section for pair production.
    */ 
    
    //double s=proc.GetCMEnergy();
    //double E0=proc.GetIncidentParticle().GetEnergy();
    
    double s = ElectronMass*ElectronMass + 2*eps*E0*(1-cos(theta));
    double beta = sqrt(1-4*ElectronMass*ElectronMass/s);
    
    // for kinematic reasons..
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
}


double ExtractPPSecondariesEnergy(Process proc) {
    /*!
        Input: incident gamma Energy E0, background photon energy eps, 
        incidence angle theta.
        Returns the energy of the produced e+ (e-)
    */

    double E0 = proc.GetIncidentParticle().GetEnergy();
    double s = proc.GetCMEnergy();
    double eps = proc.GetTargetParticle().GetEnergy();
    double theta = proc.GetInteractionAngle();
    double beta = sqrt(1-4*ElectronMass*ElectronMass/s);
    double s2 = ElectronMass*ElectronMass + 2*eps*E0*(1-(beta)*cos(theta));
    bool failed=1;
    // reInitialization to zero..
    for(int i=0; i<MC_SAMPLING; i++){
        for(int j=0; j<3; j++) MC_Sampling_Hist[i][j] = 0.;
    }


    /*
        Exponential binning of [a,b]:
        
        b = a*f^n => f = (b/a)^(1/n)
        
        For a=0.5*(1-beta)*E0 and b=0.5*(1+beta)*E0 we have f = [(1+beta)/(1-beta)]^(1/n)
    */

    double f = pow( (double)(1+beta)/(1-beta) ,(double)1./(double)MC_SAMPLING);
    int cnt = 0;
    double NormFactor = 0;

   // debug 
/*
    cerr << setprecision(10) << "E0: " << E0 << endl;
    cerr << setprecision(100) << "s: " << s << " vs " << s2 << endl;
    cerr << setprecision(100) << "betacalc: " << 4.0*ElectronMass*ElectronMass/s << " vs " <<  4.0*ElectronMass*ElectronMass/s2 <<  endl;
    cerr << setprecision(100) << "beta: " << beta << endl;
    cerr << setprecision(10) << "theta: " << theta*RAD2DEG << endl;
    cerr << setprecision(10) << "Extraction in ["<< 0.5*(1-beta)*E0 << " , " << 0.5*(1+beta)*E0 << "]" << endl;
  */  
    
    for (double Ee=f*0.5*(1-beta)*E0; Ee<0.5*(1+beta)*E0; Ee*=f){
        MC_Sampling_Hist[cnt][0] = Ee;
        MC_Sampling_Hist[cnt][1] = dSigmadE_PP(Ee, E0, eps, theta);

        NormFactor += MC_Sampling_Hist[cnt][1];
        MC_Sampling_Hist[cnt][2] = NormFactor;
            
        if(MC_Sampling_Hist[cnt][1] > 0.){
            cnt++;
        }else{
            break;
        }
    }
    
    NormFactor = (double)1./(double)NormFactor;
    //Normalize the cdf to be used for MC extraction
    for(int i=0; i<cnt; i++) MC_Sampling_Hist[i][2] *= NormFactor;
    
    /* //debug
    for(int i=0; i<cnt; i++){
       cout << setprecision(20) << MC_Sampling_Hist[i][0] << "  " << MC_Sampling_Hist[i][1] << "  " << MC_Sampling_Hist[i][1]*NormFactor << "  " << MC_Sampling_Hist[i][2]*NormFactor << endl;
    }
    cerr << "Norm: " << NormFactor << endl;
    exit(0);
    */
    
    // Sample a random between 0 and 1 and get the energy corresponding to the cdf

    double rnd; 
    double Ee = 0;

   while (failed) { 
    rnd = Uniform(0,1);
    Ee = 0;

  //  cout << rnd << " in " << MC_Sampling_Hist[0][2] << "  and " << MC_Sampling_Hist[cnt][2] << " because of " << cnt << " loops " << endl ; 
    for(int i=0; i<cnt-1; i++){
        if( MC_Sampling_Hist[i][2]<= rnd <= MC_Sampling_Hist[i+1][2] ){ 
            Ee = MC_Sampling_Hist[i][0];
            failed=0;
            break;
        }
    }
} //end while
    //cerr << "\t\treturn " << Ee << endl;
    //exit(0);
    
    /*
        Generally, one need only the Ee to be returned from this function. In fact,
        it follows exactly what it is expected (i.e. the normalized integral of
        d sigma/dE).
        However, in our code we require that, with equal probability, sometimes 
        it returns Ee other times the E0-Ee (corresponding to the energy of the 
        other particle in the pair)
    */
    if (Ee==0 || E0-Ee==0) {

    cerr << endl << setprecision(10) << "E0: " << E0 <<  " results:   " << Ee << endl;
    cerr << setprecision(100) << "s: " << s << " vs " << s2 << endl;
    cerr << setprecision(100) << "betacalc: " << 4.0*ElectronMass*ElectronMass/s << " vs " <<  4.0*ElectronMass*ElectronMass/s2 <<  endl;
    cerr << setprecision(100) << "beta: " << beta << endl;
    cerr << setprecision(10) << "theta: " << theta*RAD2DEG << endl;
    cerr << setprecision(10) << "Extraction in ["<< 0.5*(1-beta)*E0 << " , " << 0.5*(1+beta)*E0 << "]" << endl;

    cout << Ee << endl;

    for(int i=0; i<cnt-1; i++){
    cout << i << " : " <<  MC_Sampling_Hist[i][2] << " <= " << rnd <<  " <= " << MC_Sampling_Hist[i+1][2] <<endl;  
      if( MC_Sampling_Hist[i][2]<= rnd <= MC_Sampling_Hist[i+1][2] ){
         Ee = MC_Sampling_Hist[i][0];
           cout <<  i  << "   Ee = " <<  MC_Sampling_Hist[i][0] << endl;
            break;
          }
       }
       exit(0);
 
    }

    
   if( Uniform(0,1)<0.5 ){
//     cout << "ene e- = " << Ee << " and e+: " << E0-Ee << " from E0: " << E0 <<  endl;
        return Ee;
    }else{
  //    cout << "ene e- = " << E0-Ee << " and e+: " << Ee << " from E0: " << E0 << endl;
        return E0-Ee;
    }

}


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

    /*
        Exponential binning of [a,b]:
        
        b = a*f^n => f = (b/a)^(1/n)
        
        For a=(1-beta)/(1+beta)*Ee and b=Ee we have f = [(1+beta)/(1-beta)]^(1/n)
    */
    double f = pow( (double)(1+beta)/(1-beta) ,(double)1./MC_SAMPLING);
    int cnt = 0;
    double NormFactor = 0;

    /* // debug
    cerr << setprecision(10) << "Ee: " << Ee << endl;
    cerr << setprecision(100) << "s: " << s << endl;
    cerr << setprecision(100) << "beta: " << beta << endl;
    cerr << setprecision(10) << "theta: " << theta*RAD2DEG << endl;
    cerr << setprecision(10) << "Extraction in ["<< ((1-beta)/(1+beta))*Ee << " , " << Ee << "]" << endl;
    */

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
    //Normalize the cdf to be used for MC extraction
    for(int i=0; i<cnt; i++) MC_Sampling_Hist[i][2] *= NormFactor;
    
    /* //debug
    for(int i=0; i<cnt; i++){
        cout << setprecision(20) << MC_Sampling_Hist[i][0] << "  " << MC_Sampling_Hist[i][1] << "  " << MC_Sampling_Hist[i][1]*NormFactor << "  " << MC_Sampling_Hist[i][2] << endl;
    }
    cerr << "Norm: " << NormFactor << endl;
    exit(0);
    */
    
    // Sample a random between 0 and 1 and get the energy corresponding to the cdf
    double rnd; 
    double Eer = 0;

 while (failed) {
   rnd = Uniform(0,1);
   Eer=0;
    for(int i=0; i<cnt-1; i++){
        if( MC_Sampling_Hist[i][2]<= rnd <= MC_Sampling_Hist[i+1][2] ){ 
            Eer = MC_Sampling_Hist[i][0];
            failed=1;
            break;
        }
    }

}
    //cerr << "\t\treturn " << Ee << endl;
    //exit(0);
    
    return Eer;
}


double ExtractTPPSecondariesEnergy(Process proc) {
  /* approximation based on A. Mastichiadis et al., 
    Astroph. Journ. 300:178-189 (1986), eq. 30. 
This approx is valid only for   alpha >=100 
where alpha = p0*eps*costheta - E0*eps; 
for our purposes, me << E0 --> p0~ E0 --> 
alpha = E0*eps*(costheta - 1) >= 100; 

p0 = sqrt(E^2 - m2c4) --> 

alpha = sqrt(E0^2 - m^2)*eps - E0*eps;  
es. 1e18 --> sqrt(1e36 - 25.1 * 1e4) - 1e18
*/

  double E0 = proc.GetIncidentParticle().GetEnergy();
  double eps = proc.GetTargetParticle().GetEnergy();
  double Epp= 5.7e-1 * pow(eps,-0.56)*pow(E0,0.44);
//option 2: according to Lee 98: 
//the  inelasticity is 1.768*pow(s/me2)^{-3/4} that should be energy remaining to the incident particle. 
//thus the other 2 electrons equally share the remaining energy (E0 - n*E0)
  double s = proc.GetCMEnergy();
  double Epp2= E0*(1 - 1.768*pow(s/ElectronMass/ElectronMass,-3.0/4.0))/2.0;
//  cout << " in TPP function " << s << " (ene: " <<   E0<< " , eps:  " << eps << " ) " <<  Epp << "  vs " << Epp2 << endl;
  return Epp;
}


double ExtractDPPSecondariesEnergy(double E0)  {
if (E0==0) cout << "error in extracting DPP: can not be =0 " <<endl;
  //we use the same assumption of lee (i.e., all the energy goes equaly shared between only 1 couple of e+e-.
  // In DPPpaper has been shown that this approximation is valid within -1.5%
  return (double)E0/2.0;
}    
    

