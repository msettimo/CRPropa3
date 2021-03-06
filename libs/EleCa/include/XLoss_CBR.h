namespace crpropa {
/*=====================================================================

	License
	=======
		
    This file is part of xHERMES () and EleCa packages for 
    multi messenger data analysis.

	(C) Copyright 2009, 2010, 2011  Manlio De Domenico and Mariangela Settimo
	
	Author: Manlio De Domenico
	Lab. for Complex Systems, Scuola Superiore di Catania
	Universita' degli Studi di Catania, Italy
	Mail: manlio.dedomenico@ct.infn.it

	(C) Copyright 2011, 2012 Mariangela Settimo	
	Author: Mariangela Settimo
	Universaat Siegen, Germany, now at LPNHE Paris
	Mail: mariangela.settimo@gmail.com
	
=====================================================================*/

//##########################################################################
//# Cosmic Background Radiations
//##########################################################################

double CIB_Evolution_Baseline(double z){
    // Function for the CIB baseline evolution.
    // Stecker, Malkan, Scully (2006) arXiv:astro-ph/0510449v4

    double tmp = 0;
    double m = 3.1;
    double z_flat = 1.3;
    
    if( z<=z_flat ) tmp = pow(1.+z,m);
    if( z_flat<z && z<6 ) tmp = pow(1.+z_flat,m);
    if( z>6 ) tmp = 0;
    
    return tmp;
}

double CIB_Evolution_Fast(double z){
    // Function for the CIB fast evolution.
    // Stecker, Malkan, Scully (2006) arXiv:astro-ph/0510449v4

    double tmp = 0;
    double m = 4.;
    double z_flat = 1.;
    
    if( z<=z_flat ) tmp = pow(1.+z,m);
    if( z_flat<z && z<6 ) tmp = pow(1.+z_flat,m);
    if( z>6 ) tmp = 0;
    
    return tmp;
}

double CMB_Evolution(double z){ return pow(1.+z,3.); }

double CIB_Evolution(double z){ return CIB_Evolution_Fast(z); }

double CIOB_Evolution(double z){ return CIB_Evolution_Fast(z); }

double COB_Evolution(double z){ return pow(1.+z,3.); }


double URB_Evolution(double z){ 
  //from Protheroe - Bierman astro-ph:9605119 
  if (z< 0.8) return pow(1.+z,4.); 
  return pow(1+0.8,4);   // z>= z0 
}

double CMBR(double eps){
    double tmp = 0;
    
    if( eps>eps_ph_inf_cmb && eps<eps_ph_sup_cmb ){
        tmp = (K_CBR*eps*eps)/( exp((double)eps/(K_boltz*T_CMB)) -1 );
    }else{
        tmp = 0;    
    }
    
    if( isnan(tmp) ) tmp = 0;
    
    return tmp;
};

double CIBR(double eps){
    double tmp = 0;
    
    if( eps>eps_ph_inf_cib && eps<=eps_ph_sup_cib ){
        tmp = 5e-1*((2.2e-6*K_CBR*eps*eps)/( exp((double)eps/(K_boltz*T_CMB)/9.17) -1 ) + (2.e-11*K_CBR*eps*eps)/( exp((double)eps/(K_boltz*T_CMB)/128.44) -1 ));
    }else{
        tmp = 0;    
    }

    if( isnan(tmp) ) tmp = 0;
    
    return tmp;
};

double CIOBR(double eps){
    // parametrization for infrared/optical by 
    // Hopkins, A. M. & Beacom, J. F. 2006, ApJ, 651, 142
    // See Model D Finke et al, arXiv:0905.1115v2
    
    double tmp = 0;
    
     if( eps>eps_ph_inf_ciob && eps<eps_ph_sup_ciob ){
        double x = log(eps);    
        tmp= -5.32524895349885 -0.0741140642891119*x -0.252586527659431*pow(x,2.) +0.234971297531891*pow(x,3.) -0.217014471117521*pow(x,4.) -0.364936722063572*pow(x,5.) +0.0880702191711222*pow(x,6.)+0.221947767409286*pow(x,7.) +0.0445499623085708*pow(x,8.) -0.0517435600939147*pow(x,9.)-0.0295646851279071*pow(x,10.)-0.00011943632049331*pow(x,11.) +0.00461621589174355*pow(x,12.) +0.00150906100702171*pow(x,13.) +1.91459088023263e-05*pow(x,14.) -0.000110272619218937*pow(x,15.) -3.45221358079085e-05*pow(x,16.) -5.42000122025042e-06*pow(x,17.) -4.90862622314226e-07*pow(x,18.) -2.45145316799091e-08*pow(x,19.) -5.25792204884819e-10*pow(x,20.);
        tmp = 0.4*(double)exp(tmp)/eps/eps;
    }else{
        tmp = 0;    
    }

    if( isnan(tmp) ) tmp = 0;
    
    return tmp;
};

double COBR(double eps){
    double tmp = 0;
    
    if( eps>eps_ph_inf_cob && eps<eps_ph_sup_cob ){
        tmp = 1.2e-15*(K_CBR*eps*eps)/( exp((double)eps/(K_boltz*T_COB)) -1 );
    }else{
        tmp = 0;    
    }

    if( isnan(tmp) ) tmp = 0;
    
    return tmp;
};



// Universal Radio Background from Protheroe, Bierman 1996. 

double URB(double eps){
  if (eps<eps_ph_inf_urb || eps>eps_ph_sup_urb) return 0;
 
  double v=eps/h_Planck;
  double x=log10(v/1e9);
 
  double p0= -2.23791e+01;
  double p1= -2.59696e-01;
  double p2= 3.51067e-01; 
  double p3= -6.80104e-02;
  double p4= 5.82003e-01; 
  double p5= 2.00075e+00; 
  double p6= 1.35259e+00;
  double p7= -7.12112e-01;  //xbreak

  double intensity = 0;
  if (x>p7) intensity = p0 + p1*x+ p2*x*x + p3*x*x*x/(exp(p4*x) - 1) + p6+p5*x;
  else intensity =  p0 + p1*x+ p2*x*x + p3*x*x*x/(exp(p4*x) - 1);
  intensity = pow(10,intensity);
  double n_eps=0;
  n_eps = 4*cPI/(h_Planck*C_speed)*(intensity/eps);
  return n_eps/eV2J/1.0e6;

};

 
double CMIBR(double eps){
    /*! 
        Cosmic background radiation photon number density (eV^-1 cm^-3)
        as a function of ambient photon energy (eV), from CMB to Optical (COB)
        
        Ref:
        
        Funk et al, Astropart.Phys. 9 (1998) 97-103
        J.L. Puget, F.W. Stecker and J. Bredekamp, Astroph. J. 205 (1976) 638–654. 
    */    
  return CMBR(eps) + CIBR(eps);
};

double CMIOBR(double eps){
    /*! 
        Cosmic background radiation photon number density (eV^-1 cm^-3)
        as a function of ambient photon energy (eV), from CMB to Optical (COB)
        
        Ref:
        
        Funk et al, Astropart.Phys. 9 (1998) 97-103
        J.L. Puget, F.W. Stecker and J. Bredekamp, Astroph. J. 205 (1976) 638–654. 
    */    
  return CMBR(eps) + CIOBR(eps);
};


double CBR(double eps, double z){
    /*! 
        Cosmic background radiation photon number density (eV^-1 cm^-3)
        as a function of ambient photon energy (eV), from CMB to Optical (COB)
        
        Ref:
        
        Funk et al, Astropart.Phys. 9 (1998) 97-103
        J.L. Puget, F.W. Stecker and J. Bredekamp, Astroph. J. 205 (1976) 638–654. 
    */    
  return CMBR(eps)*CMB_Evolution(z) + CIOBR(eps)*CIOB_Evolution(z) + URB(eps)*URB_Evolution(z);
};

double CBR(double eps){
  return CMBR(eps) + CIOBR(eps) + URB(eps);
};


double GetEvolution(double eps, double z) {

  if (eps>= eps_ph_inf_urb  && eps <= eps_ph_sup_urb ) return URB_Evolution(z);
  if (eps>= eps_ph_inf_ciob && eps <= eps_ph_sup_ciob ) return CIOB_Evolution(z);
  if (eps>= eps_ph_inf_cmb  && eps <= eps_ph_sup_cmb) return CMB_Evolution(z);
  return 1;
};

double GetEvolution(double z, string background) {
  if (background=="CMB") return CMB_Evolution(z);
  if (background=="CIB") return CIB_Evolution(z);
  if (background=="CMIOB") return CIOB_Evolution(z);
  if (background=="CIOB") return CIOB_Evolution(z);
  if (background=="COB") return COB_Evolution(z);
  if (background=="URB") return URB_Evolution(z);

  if (background=="ALL") { // come trattare questo caso???? 
    return CMB_Evolution(z); 
 }
  return 1;
};

//------

double CBR(double eps, double z, string background) {
  double evolution=1;
  if (z==0) evolution=1;
  else evolution=GetEvolution(eps,z);

  if (background=="CMB") return CMBR(eps)*evolution;
  else if (background=="URB") return URB(eps)*evolution;
  else if (background=="CIOB") return CIOBR(eps)*evolution;
  else if (background=="CMIB") return CMIBR(eps)*evolution;
  else if (background=="CMIOB") return CMIOBR(eps)*evolution;
  else return (CMBR(eps)*CMB_Evolution(z) + CIOBR(eps)*CIOB_Evolution(z) + URB(eps)*URB_Evolution(z));
};

}
