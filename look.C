{

  TFile *file = new TFile("mc_fadc.root", "READ");
  TTree *tree = (TTree*)file->Get("tree");

  TH1F  *h1;
  TH1F  *h2;
  TH1F  *h3;

  tree->SetBranchAddress("h1",&h1);
  tree->SetBranchAddress("h2",&h2);
  tree->SetBranchAddress("h3",&h3);

  tree->GetEntry(100);

  TCanvas* canv = new TCanvas("canv","canv",900,600);
  canv->Divide(1,3);
  canv->cd(1);
  h1->Draw("hist");
  canv->cd(2);
  h2->Draw("hist");
  canv->cd(3);
  h3->Draw("hist");

}
