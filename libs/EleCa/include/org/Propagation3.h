#include "testlambda.h"

vector<double> SetTmprow(Particle p0) {
  vector<double> tmprow;
  tmprow.push_back(p0.GetId());
  tmprow.push_back(p0.GetParentId());
  tmprow.push_back(p0.GetType());
  tmprow.push_back(p0.Getz());
  
  tmprow.push_back(p0.GetEnergy());
  tmprow.push_back(p0.GetTheta());
  tmprow.push_back(p0.GetPhi());
  return tmprow;
}

vector<double> FillParticleRow(int id, int parentid, int type, double z, double E, double Theta, double Phi) {

  vector<double> tmprow01;
  tmprow01.push_back(id);
  tmprow01.push_back(parentid);
  tmprow01.push_back(type);
  tmprow01.push_back(z);
  tmprow01.push_back(E);
  tmprow01.push_back(Theta);
  tmprow01.push_back(Phi);
  return tmprow01;
}

//small angle approximation

  double GetMeanThetaBFDeflection(double Bin, double Ein, int ptype, double Lin) {
  //from D. Hooper, S. Sarkar, M. Taylor, arXiv: 0608085, 2006
  
  if (Bin ==0 || Ein==0) return 0;
  if (ptype == 0) return 0;

  double lcoher=1;
  // equation:   theta = 0.8*(1.0e20/Ein)*sqrt(Lin/10 Mpc)*sqrt(lcoher/1Mpc)*(Bin/nG);
  return 0.8*(1.0e20/Ein)*sqrt(Lin/10*lcoher)*(Bin/1.0e-9)/180.0*3.1415927;
}


//----
void Propagate(vector<double> tmprow, Process &proc) {

  double zBFTot = 0.0;
  double thetaBF=0.0;
  double CohLength=1.;  //Mpc
  const double xacc = 1.0E-6;
  const int MAXIT   = 10000;     // Maximum allowed number of iterations.

  //Secondary particles are returned as a vector.
  double BNorm =proc.GetB();
  
  double zin = tmprow.at(3);
  double z_curr=zin; 

  int sid=0;
  int sid2=0;

  bool interacted=0;
  double min_dist=1e12;
  double walkdone = 0;

  double phi1, phi2, phi3;
  phi1=phi2=phi3=0;   

  double T1=0; 
  double T2 = 0;
  double T3 = 0;
  double E1= 0;
  double E2 = 0; 
  double E3 = 0; 


  double stepsize =0; //in Mpc
  double Ecurr = tmprow.at(4);
  int type=(int)tmprow.at(2);
  //  vector< vector<double> > outrow;
  vector<double> tmprow1;
  //  vector<double> tmprow2;
  // vector<double> tmprow3;

   //some init values
  Particle curr_particle;
  curr_particle.SetParticle(tmprow);  


  Particle pt;
  double tmp_lambda=0;
  double Etarget=0;
  double Etarget2=0;
  double min_dist2=0;
  double R=0;
  double R2=0;
  double Elast=0;

  zin=curr_particle.Getz();
  z_curr = zin;
  int iloc=0;
  int iloc0=(int)(abs(zin-z0ph)/dprof);

  double theta_deflBF=0; 
  double smintmp=1.0e3;

  R = Uniform(0.0,1.0);
  R2 = Uniform(0.0,1.0);
  proc.SetIncidentParticle(curr_particle);
  vector<double> EtargetAll;


  if (type==0 || type==1 || type==-1)
     EtargetAll=GetEtarget(proc,curr_particle);
  else cout << "In Propagation 3: something wrong with ptype" << endl;
 
   min_dist = ExtractMinDist(proc, pt, curr_particle.GetType(), R, R2, EtargetAll);
 
   interacted=0;
   double dz=0;
   double zpos=zin;
   //cout << "\n\nreading now the particle!! it's at z0 = " << zpos <<  endl;
   double corrB_factor=0;
   double realpath=0;

   double min_dist_last = min_dist;
   if (debug) cout << "starting propagation... min_dist_last: " << min_dist_last << endl; 

   while (!interacted) {

     proc.SetInteractionAngle(cPI);
     theta_deflBF=0;

     realpath=0.1*min_dist;
        //stepsize=0.1*min_dist; 
       // dz=Mpc2z(stepsize); 
 
     //MAGNETIC FIELD!!
      theta_deflBF=GetMeanThetaBFDeflection(BNorm,curr_particle.GetEnergy(),curr_particle.GetType(),min_dist);
      corrB_factor = cos(theta_deflBF);

      
      stepsize=realpath*corrB_factor;
      //realpath=stepsize/corrB_factor;
      dz=Mpc2z(stepsize);

      if (debug && BNorm>0)
       cout << "z_curr " << z_curr << ", type: " << curr_particle.GetType() << ", Ecurr " << curr_particle.GetEnergy() << "  = " << min_dist << " Mpc -->  theta: " << theta_deflBF*180.0/3.1415927 << "  " << corrB_factor << endl;

      // realpath=stepsize/corrB_factor;
      // cout << stepsize << " ( = " << dz << " ), " << realpath  <<endl; 

      //stepsize = path undeflected --> dz
      //realpath = path with B deflection  -->dz_defl
 
      if (debug) cout << "done : " << walkdone << " + " << realpath << " = " << walkdone+realpath << " has to be:  " << min_dist << " Mpc. Now @ z= " << zpos << " vs nexz = " <<  zpos-Mpc2z(stepsize) <<endl;

      //some checks
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
	// cout << " zpos - dz = " << zpos << " - " << dz << " < 0!  correcting path from " <<  dz << " to ";
        dz=zpos;
        stepsize=z2Mpc(dz);
        realpath=stepsize/corrB_factor;

	//  cout << dz << ". (theta, corrB_factor,path) are correspondingly corrected from " << endl;
	//  cout << theta_deflBF << "," << corrB_factor << "," << realpath << " Mpc to: " << endl;
       // theta_deflBF=GetMeanThetaBFDeflection(BNorm,curr_particle.GetEnergy(),curr_particle.GetType(), stepsize); 
       // theta_deflBF=GetMeanThetaBFDeflection(BNorm,curr_particle.GetEnergy(),curr_particle.GetType(), min_dist); 
       // corrB_factor = cos(theta_deflBF);// 1 + theta_deflBF*theta_deflBF/2.0;
       // realpath=stepsize/corrB_factor;
	// cout << theta_deflBF << "," << corrB_factor << "," << realpath << endl << endl;
      } 

      // travell the calculated path!
      zpos-=dz;
      walkdone+=realpath;

      //loose energy
      Elast=Ecurr;
      RKstepsize=2.5e-5;

      double adiab = Ecurr - AdiabaticELoss(zpos+Mpc2z(realpath), zpos+dz, Ecurr);  //AdiabaticLoss(znew, z0,E0);
      Ecurr = EnergyLoss1D(Ecurr,zpos+Mpc2z(realpath),zpos,BNorm*abs(type),RKstepsize);  
    //here we need to use zfin = path to calculate the real energy loss. 
      

    /*in the function EnergyLoss1D I should loose energy according to adiabatic loss only for the path dz' while the synch should be considered for the realpath
     zpos e' gia' aggiornato. qsta cosa perche' siccome calcolo l'energia su un dz' > dz voglio evitare che dz' < 0 allora parto dalla posizione corrente e assumo 
che la particella torni indietro nel tempo ad una posizione leggermente maggiore di qlla vera. l'altra opzione sarebbe di fare l'integrale in r, ma se nn cambio la 
funzione frs e' piu' facile la futura implementazione con Hermes. 
    */

//     cout << "zpos,walkdone,diff" <<type <<  "  " << zpos <<"  " << Mpc2z(walkdone) << "  " << Mpc2z(walkdone)-dz <<  "  " << Elast << " " << Ecurr << "  " << adiab << endl;
      Eoutput+=(Elast-Ecurr);
      EoutputLoss+=(Elast-Ecurr);
     // if (type!=0) cout << "type : " << type << " B= " << BNorm << " Elast " << Elast << " - Ecurr: " << Ecurr << " = " << Elast-Ecurr <<  endl;
      z_curr = zpos; 

    if (z_curr>0 && Ecurr<Ethr) {       
         NsecBT++;
         Eoutput+=Ecurr;
         EoutputBT+=Ecurr;
         return;
     }

     curr_particle.Setz(z_curr);
     curr_particle.SetEnergy(Ecurr);
     proc.SetIncidentParticle(curr_particle);
     proc.SetCMEnergy();
     proc.SetLimits();

 //update of lambda to be done only very close to the interaction point
  //cout << " Updating MinDist for particle: " << type << endl;
  min_dist = ExtractMinDist(proc, pt, curr_particle.GetType(), R, R2,EtargetAll);
         
  if (interacted == 1) {
       Einput+=proc.GetTargetParticle().GetEnergy(); 
       EinputEps+=proc.GetTargetParticle().GetEnergy(); 
       //cout << " Einput: " << Einput << endl;

      //filling the shw profile - fill up to (ipos-1) because in ipos the new secondaries will be produced and counted. 
      iloc=(int)(abs(z_curr-z0ph)/dprof);
      for (int intg=iloc0; intg<iloc; intg++) { 
       if (curr_particle.GetEnergy()>=Ethr) { 
        if (curr_particle.GetType()==0) Nph[intg]++; else Nem[intg]++;}
      }//fill the shw profile

     }

     if (z_curr<=0) {
//    cout << "particle reaches the ground. setting zfin and Efin" << endl; 
      tmprow.push_back(z_curr);
      tmprow.push_back(curr_particle.GetEnergy());
      ParticleAtGround.push_back(tmprow);

      //cout << "outrow has: " <<  outrow.size() <<  " entries with last tmprow.size = " << tmprow.size() << endl;
      if (curr_particle.GetEnergy()==0) cout << "WARNINGGGG!! for particle " << type << "reaching the ground E = 0" << endl;  
      Eoutput+=curr_particle.GetEnergy();
      EoutputG+=curr_particle.GetEnergy();

      //filling the shw profile
      iloc=(int)(abs(z_curr-z0ph)/dprof);
      for (int intg=iloc0; intg<=iloc; intg++) { 
       if (curr_particle.GetEnergy()>=Ethr) { 
        if (curr_particle.GetType()==0) Nph[intg]++; else Nem[intg]++;}
      }//fill the shw profile
      return;
    } 

   }//end while


  if (interacted==1) {
//    cerr << "******** producing secondary particles according to " << proc.GetName() << " process *********** " <<endl;

    if (proc.GetName()=="PP") {
       T1=T2=cPI; 
       E1= ExtractPPSecondariesEnergy(proc);  //energy of e- 
       E2 = Ecurr-E1; 
       //cout << "CheckEsec: " << proc.GetName() << "@ E0 = " << Ecurr << " and products : " << (double)E1/Ecurr*100 << " and " << (double)E2/Ecurr*100 << "[%]"; 
      //if (E1<1e16 || E2<1e16 ) cout << " ---> at least 1below Ethr " <<endl; else cout << endl; 

      if(E1==0 || E2==0) cerr << "ERROR in PP process:  E : "<<Ecurr << "  " <<  E1<<" "<<E2<<endl;

      // struttura di tmprow: ID | ParentId | type | zin | Ein | theta | phi | zfin | Efin
      //zfin ed Efin vengono assegnate al prox step, quando la particella viene processata
       //la prima entry dovrebbe essere un id della matrice ParticleMatrix. Lo riassegno dopo, a questo livello l'sdid = 
       //1 e 2 per le nuove particelle e -1 se la particella e' recoiled. 
       //nel main riassegno poi i valori giusti = lastid + tmprow.at(0) e 
       //if (tmprow.at(0) < 0) tmprow.at(0)=tmprow.at(1);

       //per riconoscere le particelle recoiled controllo che tmprow.at(0) == tmpr.at(1);

       if (E1<Ethr) {
	 NsecBT++;
	 Eoutput+=E1;
	 EoutputBT+=E1;
       }
      else {
	tmprow1 = FillParticleRow(1, curr_particle.GetId(), -1, z_curr, E1, T1, phi1);
	ParticleMatrix.push_back(tmprow1);
      }
       if (E2<Ethr) { 
	 NsecBT++; 
	 Eoutput+=E2;
	 EoutputBT+=E2;
       }
     else {
       tmprow1 = FillParticleRow(2, curr_particle.GetId(), +1, z_curr, E2, T2, phi2);
       ParticleMatrix.push_back(tmprow1);
       }
     } //if PP
    else if (proc.GetName()=="DPP") {
       T1=T2=cPI; 
      // E1 = E2 = Ecurr/4.0;
       E1 = Ecurr/2.0; //from lee approx: all the E goes in one pair e+e-
       E2 = 0; // that correspond to (Ecurr - E1 - E2)/2.0;
       NsecBT+=2; //because of this second pair
       if(E1==0 )cerr<<"ERROR in DPP process E : "<<E1<<endl;

       if (E1<Ethr) { NsecBT+=2; 
         Eoutput+=2*E1;//4*E1;
         EoutputBT+=2*E1;//4E1
       } else {
       tmprow1 = FillParticleRow(1, curr_particle.GetId(), -1, z_curr, E1, T1, phi1);
       ParticleMatrix.push_back(tmprow1);
       tmprow1 = FillParticleRow(2, curr_particle.GetId(), +1, z_curr, E1, T2, phi2);
       ParticleMatrix.push_back(tmprow1);
      // tmprow1 = FillParticleRow(3, curr_particle.GetId(), -1, z_curr, E2, T1, phi1);
      // ParticleMatrix.push_back(tmprow1);
      // tmprow1 = FillParticleRow(4, curr_particle.GetId(), +1, z_curr, E2, T2, phi2);
      // ParticleMatrix.push_back(tmprow1);
       }
    } 
    else if (proc.GetName()=="ICS") { 

      T1= T2 =cPI; 
      E2= ExtractICSSecondariesEnergy(proc); //the function return the energy of the recoil electron
      E1= Ecurr-E2;   //energy of photon is the total initial energy - the recoiled one
      if(E1==0 || E2==0) cerr<<"ERROR in ICS process E : "<<E1<<" "<<E2<<endl;

      //cout << "CheckEsec: " << proc.GetName() << "@ E0 = " << Ecurr << " and products : " << (double)E1/Ecurr*100 << " and " << (double)E2/Ecurr*100 << " [%] " ;
     // if (E1<1e16 || E2<1e16 ) cout << " ---> at least 1below Ethr " <<endl; else cout << endl; 
      if (E1<Ethr) {
	NsecBT++;
	 Eoutput+=E1;
	 EoutputBT+=E1; 
      }
      else {
	tmprow1 = FillParticleRow(1,curr_particle.GetId(), 0, z_curr, E1, T1, phi1);  
	ParticleMatrix.push_back(tmprow1);
      } 
      if (E2<Ethr) { 
	NsecBT++;
	 Eoutput+=E2;
	 EoutputBT+=E2;
      }
	else {
	  tmprow1 = FillParticleRow(curr_particle.GetId(), curr_particle.GetId(), curr_particle.GetType() , z_curr, E2, T2, phi2);
	  ParticleMatrix.push_back(tmprow1);
	}
    } //if ics
    else if (proc.GetName()=="TPP") {
      E1 = E2 = ExtractTPPSecondariesEnergy(proc);
      E3 = Ecurr - E1 - E2; 
     if(E1==0 || E2==0)cerr<<"ERROR in TPP process E : "<<E1<<" "<<E2<<endl;

      //cout << "CheckEsec: " << proc.GetName() << "@ E0 = " << Ecurr << " and products : " << (double)E1/Ecurr*100 << " and " << (double)E2/Ecurr*100 << " and " <<  E3/Ecurr*100 << "[\%]";
 //     if (E1<1e16 || E2<1e16 || E3<1e16 ) cout << " ---> at least 1below Ethr " <<endl; else cout << endl; 
     
      T1 = T2 = T3= cPI;//tmpvalues  calcolare gli altri

      sid=0;
      sid2=0;
      if (E1<Ethr) {
	NsecBT++;
	Eoutput+=E1; 
	EoutputBT+=E1; 
      } 
      else {
	tmprow1 = FillParticleRow(1,curr_particle.GetId(),-1, z_curr, E1,T1,phi1); 
	ParticleMatrix.push_back(tmprow1);
      }
      if (E2<Ethr) {
	NsecBT++;
	Eoutput+=E2; 
	EoutputBT+=E2; 
      } else {
	tmprow1 = FillParticleRow(2,curr_particle.GetId(), 1, z_curr, E2,T2,phi2);
	ParticleMatrix.push_back(tmprow1);
      }
      if (E3<Ethr) {
	NsecBT++; 
	Eoutput+=E3; 
	EoutputBT+=E3; 
      } 
      else {
	tmprow1 = FillParticleRow(curr_particle.GetId(),curr_particle.GetId(),curr_particle.GetType(),z_curr,E3,T3,phi3);
	ParticleMatrix.push_back(tmprow1);
      }
    }//if tpp

  } //end if interacted
  
  
  return;

} //end propagate();
