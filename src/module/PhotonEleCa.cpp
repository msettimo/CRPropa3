#include "crpropa/module/PhotonEleCa.h"
#include "../libs/EleCa/include/Propagation.h"

using namespace crpropa;

void PhotonEleCa::process(Candidate candidate, double Ethr)  {
 
  Particle p0(candidate);
  Particle p1;

  p0.SetB(0);
  z0ph=candidate.getRedshift();
  E0ph=candidate.current.getEnergy();

  vector<Particle> ParticleAtMatrix;
  vector<Particle> ParticleAtGround;
  ParticleAtMatrix.push_back(p0);

  EGround.clear();
  typeGround.clear();

  int pcont=0;
 
  while (ParticleAtMatrix.size()>0) {
    
    p1 = ParticleAtMatrix.at(0);

   if (!p1.IsGood()) {
     ParticleAtMatrix.erase(ParticleAtMatrix.begin()); 
    } else {
      pcont++;
      Propagate(p1, ParticleAtMatrix, ParticleAtGround);
     ParticleAtMatrix.erase(ParticleAtMatrix.begin()); 
    }//end else

 
  } 

 
  WriteOutput(p0, ParticleAtGround,0);
}

  std::string PhotonEleCa::getDescription()  {
  	std::stringstream s;
  	s << "PhotonEleCa";
  	return s.str();
  }

