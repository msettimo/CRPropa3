namespace crpropa {
 
  void WriteOutput(bool spectrum=0) {

    int NsecG=ParticleAtGround.size();
    vector<double> EGround;
    vector<int> typeGround;
    
    for (int ig=0; ig<NsecG;++ig) {
      EGround.push_back(ParticleAtGround.at(ig).GetEnergy());
      typeGround.push_back(ParticleAtGround.at(ig).GetType());
    }
    if (spectrum==1) {
      double emin=7.0;
      double dE=(24.0-7.0)/170.0;
      int dN[170];int ipos=0;
      
      for (int j=0; j<170; ++j) dN[j]=0;
      
      for (int h=0; h<NsecG; ++h) {
	ipos = (int)((log10(EGround.at(h))-emin)/dE);
	dN[ipos]++; 
      }
      for (int j=0; j<170; ++j) 
	if (dN[j]>0) cout << j << " " << dN[j] ; 
    }
    else {
      for (int h=0; h<NsecG; ++h) {
	cout << EGround.at(h) <<  " " << " " << typeGround.at(h) << "  " ;
      }
    }
    cout << endl;
 
  }//end WriteLine
 
}
