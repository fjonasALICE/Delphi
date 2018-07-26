{

  TFile *f_all = new TFile("all_normalized.root");
  // TFile *f_all = new TFile("WeakBoson_8000GeV_noMPI_normalized.root");
  f_all->ls();

  TH2D *h2 = (TH2D*)f_all->Get("h2_electron_pt_topMotherID");

  vector<TH1D*> vec_all;

  TFile *f_out = new TFile("projections.root","RECREATE");
  TCanvas *c = new TCanvas("c_electron_mothers_pt","electron pt by topMotherID",2000,1000);
  leg = new TLegend(0.5,0.3,0.9,0.9);
  for( int i = 1; i <= h2->GetNbinsX(); i++){
    vec_all.push_back( h2->ProjectionY(h2->GetXaxis()->GetBinLabel(i), i, i) );
    vec_all.at(i-1)->Write();

    vec_all.at(i-1)->SetLineWidth(2);

    if( i == 1){
      vec_all.at(i-1)->SetLineColor(kBlack);
      vec_all.at(i-1)->SetLineStyle(1);}
    else if( i == 2){
      vec_all.at(i-1)->SetLineColor(kBlack);
      vec_all.at(i-1)->SetLineStyle(3);}
    else if( i == 3){
      vec_all.at(i-1)->SetLineColor(kBlack);
      vec_all.at(i-1)->SetLineStyle(7);}
    else if ( i == 4 )
      vec_all.at(i-1)->SetLineColor(kBlue+(i-4)*2);
    else if ( i <= 6 )
      vec_all.at(i-1)->SetLineColor(kGreen+(i-5)*2);
    else if ( i <= 8 )
      vec_all.at(i-1)->SetLineColor(kYellow+(i-7)*2);
    else if ( i <= 10 )
      vec_all.at(i-1)->SetLineColor(kMagenta+(i-9)*2);
    else if ( i == 11 )
      vec_all.at(i-1)->SetLineColor(kRed+(i-11)*2);
    else if ( i == 12 )
      vec_all.at(i-1)->SetLineColor(kOrange+1);
    else if ( i == 13 )
      vec_all.at(i-1)->SetLineColor(kGray+1);
    else if ( i <= 15 ){
      vec_all.at(i-1)->SetLineColor(kCyan+(i-14)*2);
      vec_all.at(i-1)->SetLineStyle(3);}
    else {
      vec_all.at(i-1)->SetLineColor(kMagenta+(i-16)*2);
      vec_all.at(i-1)->SetLineStyle(3);}

    if( i == 1)
      vec_all.at(i-1)->Draw();
    else
      vec_all.at(i-1)->Draw("same");

    vec_all.at(i-1)->SetMarkerSize(0);
    leg->AddEntry(vec_all.at(i-1),Form("%s",h2->GetXaxis()->GetBinLabel(i)),"le");
  }


  leg->Draw();
  gPad->Update();
  c->Write();
 

}
