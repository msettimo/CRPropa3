{
  gROOT->Reset();
  TFile g("ALLInteractionLengths.root");
  TTree *tree = (TTree*)g->Get("table");
  int nentries=tree->GetEntries();
  ofstream myfile;
  myfile.open("ALLInteractionLengths.dat");
  double Etab, PPle[360], ICSle[360], TPPle[360], DPPle[360];

  tree->SetBranchAddress("Etab",&Etab);
  tree->SetBranchAddress("PPlength",&PPle);
  tree->SetBranchAddress("ICSlength",&ICSle);
  tree->SetBranchAddress("TPPlength",&TPPle);
  tree->SetBranchAddress("DPPlength",&DPPle);

  for (int i=0; i<nentries; ++i) {
    tree->GetEntry(i);
    myfile << Etab << "  " << PPle[0] << " " <<  ICSle[0] << "  " << DPPle[0] << " " << TPPle[0] << endl;
    }
  myfile.close();
}
