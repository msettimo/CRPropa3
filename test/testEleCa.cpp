#include "crpropa/module/PhotonEleCa.h"

#include "gtest/gtest.h"
#include <fstream>

using namespace std;
using namespace crpropa;
namespace crpropa {

  TEST(testEleCa, step) {
    InitBkgArray("ALL");
    ReadTables("ALL");
    double z1 = 0.001;
    crpropa::Candidate c;
    PhotonEleCa epp1;
    epp1.fdN.resize(170);

    for ( int zi=0; zi<1; zi++) {
      z1= 0.001 + zi*0.001;
      double E = pow(10,20);
      
      c.setRedshift(z1);
      c.current.setId(22);
      
      for (int i=0; i<170; ++i) epp1.fdN[i]=0;
      int ntot=100;
    
      for (int i = 0; i < ntot; i++) {
	if (i%100==0) cout << " zi : " << z1 << " E : " << E << " " << i << " / " << ntot << endl; 
	if (E>1e23) continue;     
	c.current.setEnergy(E);
	epp1.process(c, epp1.GetEthr());
      }
      
     
      double dE = (24-7)/170.;
      char *outname=new char[100];
      sprintf(outname,"flux%2.0f_z%2.3f.txt",log10(E),z1);
      ofstream fout(outname,ios::out);
      if (!fout.is_open()) {
	cout <<  "Unable to open output file! exiting... " << endl;
	return;
      }
      cout << " writing " << outname;
      for (int j=0; j<170; ++j) 
	fout << (7+j*dE) << " " << epp1.fdN[j] << endl; 
      fout.close();
      cout << "  done " << endl;
    }//end if z
  }//end propagate

  int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
  }

} // namespace crpropa
