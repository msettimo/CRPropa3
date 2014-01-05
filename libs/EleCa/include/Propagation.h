using namespace crpropa;
void PhotonEleCa::Propagate(Particle &curr_particle, vector<Particle> &ParticleAtMatrix, vector<Particle> &ParticleAtGround) {
 
  double Ethr = fEthr;
  double theta_deflBF=0.0;
  double BNorm =curr_particle.GetB();
  
  double zin = curr_particle.Getz();
  double Ein = curr_particle.GetEnergy();
  int type = curr_particle.Getz();

  int wi_last = curr_particle.GetWeigth();

  double z_curr=zin; 
  double Ecurr=Ein; 

  bool interacted=0;
  double min_dist=1e12;
  double walkdone = 0;

  double E1 = 0;
  double E2 = 0; 
  double E3 = 0; 

  double stepsize =0;
  double Elast=0;

  double R = Uniform(0.0,1.0);
  double R2 = Uniform(0.0,1.0);
  bool fast = 1;

  Process proc;
  proc.SetIncidentParticle(curr_particle);
 
  vector<double> EtargetAll=GetEtarget(proc,curr_particle);
  min_dist = ExtractMinDist(proc, curr_particle.GetType(), R, R2, EtargetAll);

  interacted=0;
  double dz=0;
  double zpos=zin;

  double corrB_factor=0;
  double realpath=0;
  
  double min_dist_last = min_dist;
  //   debug = 1;
  if (debug) cout << "starting propagation... min_dist_last: " << min_dist_last << endl; 
  
  while (!interacted) {

     proc.SetInteractionAngle(cPI);
     theta_deflBF=0;
     realpath=0.1*min_dist;
 
      theta_deflBF=GetMeanThetaBFDeflection(BNorm,curr_particle.GetEnergy(),curr_particle.GetType(),min_dist);
      corrB_factor = cos(theta_deflBF);

      stepsize=realpath*corrB_factor;
      dz=Mpc2z(stepsize);
     
      if (debug && BNorm>0)
       cout << "z_curr " << z_curr << ", type: " << curr_particle.GetType() << ", Ecurr " << curr_particle.GetEnergy() << "  = " << min_dist << " Mpc -->  theta: " << theta_deflBF*180.0/3.1415927 << "  " << corrB_factor << endl;

      if (debug) cout << "done : " << walkdone << " + " << realpath << " = " << walkdone+realpath << " has to be:  " << min_dist << " Mpc. Now @ z= " << zpos << " vs nexz = " <<  zpos-Mpc2z(stepsize) <<endl;

      if ((walkdone+realpath) > min_dist) { 
       if (debug) 
          cout << " walkdone + realpath > min_dist: correcting realpath from" << realpath << " to "<< min_dist-walkdone  << " and the stepsize is changed from " << dz << " to " ;
         realpath=min_dist-walkdone;
         stepsize=realpath*corrB_factor;
         dz=Mpc2z(stepsize);
	 if (debug) cout << dz << endl;
         interacted=1;
      }

      if (zpos-dz <= 0) {
        dz=zpos;
        stepsize=z2Mpc(dz);
        realpath=stepsize/corrB_factor;
      }

      zpos-=dz;
      walkdone+=realpath;
      Elast=Ecurr;

      double adiab = Ecurr - AdiabaticELoss(zpos+Mpc2z(realpath), zpos+dz, Ecurr); 
     
      if (type==0 || type==22)
        Ecurr = EnergyLoss1D(Ecurr,zpos+Mpc2z(realpath),zpos,0);  
      else 
	Ecurr = EnergyLoss1D(Ecurr,zpos+Mpc2z(realpath),zpos,BNorm);  

      z_curr = zpos; 

      curr_particle.Setz(z_curr);
      curr_particle.SetEnergy(Ecurr);

      if (z_curr>0 && Ecurr<Ethr) {
	return;     
      }
      if (z_curr<=0) {
	ParticleAtGround.push_back(curr_particle);	
	return;
      } 

      proc.SetIncidentParticle(curr_particle);
      proc.SetCMEnergy();
      proc.SetLimits();
      //      vector<double> EtargetAll=GetEtarget(proc,curr_particle);
      min_dist = ExtractMinDist(proc, curr_particle.GetType(), R, R2,EtargetAll);
   }//end while


  if (interacted==1) {
    if (debug)  cerr << "******** producing secondary particles according to " << proc.GetName() << " process *********** " <<endl;

    if (proc.GetName()=="PP") {

      E1=ExtractPPSecondariesEnergy(proc);

      if(E1==0 || E1==Ecurr) cerr << "ERROR in PP process:  E : "<< Ecurr << "  " <<  E1<<" "  << endl;
      
      if (E1>Ethr) {
	Particle pp(11,E1,z_curr);
	pp.SetWeigth(wi_last);
	ParticleAtMatrix.push_back(pp);
      }
      if (Ecurr-E1 > Ethr) {
	Particle pe(-11,Ecurr-E1,z_curr);
	pe.SetWeigth(wi_last);
	ParticleAtMatrix.push_back(pe);
      }

      return;
     } //if PP
    else if (proc.GetName()=="DPP") {
      E1 = (Ecurr - 2*ElectronMass)/2.0; 
       if(E1==0 )cerr<<"ERROR in DPP process E : "<<E1<<endl;

       if (E1 > Ethr) {
	 Particle pp(11,E1,z_curr);
	if (fast==1) 
	  pp.SetWeigth(wi_last*2);
	else {
	  pp.SetWeigth(wi_last);
	  Particle pe(-11,E1,z_curr);
	  pe.SetWeigth(wi_last);
	  ParticleAtMatrix.push_back(pe);
	}
	ParticleAtMatrix.push_back(pp);
       }
   
       return;
    }//end if DPP 
    else if (proc.GetName()=="ICS") { 

      E1 = ExtractICSSecondariesEnergy(proc); 
      E2 = Ecurr-E1;  
      if(E1==0 || E2==0) cerr<<"ERROR in ICS process E : "<<E1<<" "<<E2<<endl;

      if (E1>Ethr) {
	Particle pp(curr_particle.GetType(),E1,z_curr);
	pp.SetWeigth(wi_last);
	ParticleAtMatrix.push_back(pp);
      } 
      if (E2 > Ethr) { 
	Particle pg(22,E2,z_curr);
	pg.SetWeigth(wi_last);
	ParticleAtMatrix.push_back(pg);
	}

      return;
    } //if ics
    else if (proc.GetName()=="TPP") {
      E1 = E2 = ExtractTPPSecondariesEnergy(proc);
      E3 = Ecurr - E1 - E2; 
     if(E1==0 || E2==0 || E3==0)cerr<<"ERROR in TPP process E : "<<E1<<" "<<E2<<endl;
    
      if (E1 > Ethr) {
	Particle pp(11,E1,z_curr);
	if (fast==1) 
	  pp.SetWeigth(wi_last*2);
	else { 
	  pp.SetWeigth(wi_last);
	  Particle pe(-11,E1,z_curr);
	  pe.SetWeigth(wi_last);
	  ParticleAtMatrix.push_back(pe);
	}
	ParticleAtMatrix.push_back(pp);
      }
      if (E3>Ethr) {
	Particle psc(curr_particle.GetType(),E3,z_curr);
	psc.SetWeigth(wi_last);
	ParticleAtMatrix.push_back(psc);
      }
      return;
    }
  } 
 
  return;

};  
