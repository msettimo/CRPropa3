

// ######################################################
// # Photo-Pair production cross-section (Lee 1996). 
// ######################################################

double PhotonPPCrossSection(double eCM) {
  if (eCM<4*ElectronMass*ElectronMass) return 0; 
  //qsta cross section e' ottenuta integrando su eps
  double betapar = sqrt(1-4*ElectronMass*ElectronMass/eCM);
  return 3./16. *SigmaThompson * (1-betapar*betapar)*( (3-pow(betapar,4))*log((1+betapar)/(1-betapar)) -2*betapar*(2-betapar*betapar) );
}



// ######################################################
// # Triple Pair production cross-section (Lee 1996)
// ######################################################

double TPPCrossSection(double s) {
  if (s<4*ElectronMass*ElectronMass) return 0;
  else return SigmaThompson*3*AlphaFineStruct/(8.0*cPI)*(28/9.0*log(s/(ElectronMass*ElectronMass)) - 218/27.0);
  
}


// ######################################################
// # Compton scattering from Lee --- OK!! nn toccare piu! 
// ######################################################

double PhotonICSCrossSection(double s) {
  double sigma=SigmaThompson;
  if (s > (ElectronMass*ElectronMass + 2*ElectronMass*feps_sup)) {  
    double betapar = (s - ElectronMass*ElectronMass)/(s + ElectronMass*ElectronMass); 
    double k = 3/8.*SigmaThompson*ElectronMass*ElectronMass/s/betapar;
    double A = 2.0/(betapar*(1+betapar)) *(2.0 + 2.0*betapar - betapar*betapar -2.0*betapar*betapar*betapar);
    double B = - 1.0/(betapar*betapar)*(2.0-3.0*betapar*betapar - betapar*betapar*betapar)*log((1+betapar)/(1-betapar));
   sigma= k*(A + B);
   } 
  return sigma;
}

