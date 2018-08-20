void cutFluctuations(TString dirName="bktmin3_bsup23_radfac50_betaZ0.435_noMPI", TString outFileName="cutFluctuations_default_output.root", double vetoFac=70.,const int nFiles = 3000){

  vector <TString> vec_rootInFileName;
  vector <TString> vec_vetoName;
  TString strTemp = "";
  for(int i = 1; i <= nFiles; i++){
    strTemp = dirName;
    strTemp = strTemp.Append(Form("/%d.root",i));
    vec_rootInFileName.push_back(strTemp);
  }
  for(int i = 0; i < vec_rootInFileName.size(); i++){
    printf("%s\n",vec_rootInFileName.at(i).Data());
  }


  TH1::AddDirectory(kFALSE);
  vector <TH1D*> vec_h_isodirectphoton_pt_central;
  vector <TH1D*> vec_h_dPhiJetGamma_central;
  vector <TH1D*> vec_h_xObs_pGoing_central;

  vector <TFile*> vec_rootInFile;
  TH1D *tempHist = 0x0;
  for(int i = 0; i < vec_rootInFileName.size(); i++){
    vec_rootInFile.push_back(TFile::Open(vec_rootInFileName.at(i).Data()));
    if( i == 0) vec_rootInFile.at(i)->ls();

    tempHist = (TH1D*)vec_rootInFile.at(i)->Get("h_isodirectphoton_pt_central");
    vec_h_isodirectphoton_pt_central.push_back((TH1D*)tempHist->Clone("h_isodirectphoton_pt_central"));
    tempHist = (TH1D*)vec_rootInFile.at(i)->Get("h_dPhiJetGamma_central");
    vec_h_dPhiJetGamma_central.push_back((TH1D*)tempHist->Clone("h_dPhiJetGamma_central"));
    tempHist = (TH1D*)vec_rootInFile.at(i)->Get("h_xObs_pGoing_central");
    vec_h_xObs_pGoing_central.push_back((TH1D*)tempHist->Clone("h_xObs_pGoing_central"));

    vec_rootInFile.at(i)->Close();
  }

//----------------------------------------------------------------------
  //----------------------------------------------------------------------
  TH1D *histMeanBin_isodirectphoton_pt_central = (TH1D*)vec_h_isodirectphoton_pt_central.at(0)->Clone("histMeanBin_isodirectphoton_pt_central");
  histMeanBin_isodirectphoton_pt_central->SetTitle("histMeanBin_isodirectphoton_pt_central");
  histMeanBin_isodirectphoton_pt_central->Reset();
  TH1D *histMeanBin_dPhiJetGamma_central = (TH1D*)vec_h_dPhiJetGamma_central.at(0)->Clone("histMeanBin_dPhiJetGamma_central");
  histMeanBin_dPhiJetGamma_central->SetTitle("histMeanBin_dPhiJetGamma_central");
  histMeanBin_dPhiJetGamma_central->Reset();
  TH1D *histMeanBin_xObs_pGoing_central = (TH1D*)vec_h_xObs_pGoing_central.at(0)->Clone("histMeanBin_xObs_pGoing_central");
  histMeanBin_xObs_pGoing_central->SetTitle("histMeanBin_xObs_pGoing_central");
  histMeanBin_xObs_pGoing_central->Reset();
  double tempMeanYield;

  for(int iBin = 0; iBin < vec_rootInFileName.size(); iBin++){
    tempMeanYield = 0.;
    if(vec_h_isodirectphoton_pt_central.at(iBin)->GetBinCenter(iBin) < 15. ||
       vec_h_isodirectphoton_pt_central.at(iBin)->GetBinCenter(iBin) > 30.){
      continue;
    }
    for(int i = 0; i < vec_rootInFileName.size(); i++){
      tempMeanYield += vec_h_isodirectphoton_pt_central.at(i)->GetBinContent(iBin);
    }
    tempMeanYield /= vec_rootInFileName.size();
    histMeanBin_isodirectphoton_pt_central->SetBinContent(iBin,tempMeanYield);
  }

  for(int iBin = 0; iBin < vec_rootInFileName.size(); iBin++){
    tempMeanYield = 0.;
    for(int i = 0; i < vec_rootInFileName.size(); i++){
      tempMeanYield += vec_h_dPhiJetGamma_central.at(i)->GetBinContent(iBin);
    }
    tempMeanYield /= vec_rootInFileName.size();
    histMeanBin_dPhiJetGamma_central->SetBinContent(iBin,tempMeanYield);
  }

  for(int iBin = 0; iBin < vec_rootInFileName.size(); iBin++){
    tempMeanYield = 0.;
    if(vec_h_xObs_pGoing_central.at(iBin)->GetBinCenter(iBin) < 0.0 ||
       vec_h_xObs_pGoing_central.at(iBin)->GetBinCenter(iBin) > 0.03){
      continue;
    }
    for(int i = 0; i < vec_rootInFileName.size(); i++){
      tempMeanYield += vec_h_xObs_pGoing_central.at(i)->GetBinContent(iBin);
    }
    tempMeanYield /= vec_rootInFileName.size();
    histMeanBin_xObs_pGoing_central->SetBinContent(iBin,tempMeanYield);
  }




  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  TCanvas *cSpectrum = new TCanvas("cSpectrum","cSpectrum",1500,1700);
  cSpectrum->Divide(2,3);
  cSpectrum->cd(1);
  gPad->SetLogy(1);
  for(int i = 0; i < vec_rootInFileName.size(); i++){
    if (i==0){
      vec_h_isodirectphoton_pt_central.at(i)->GetXaxis()->SetRangeUser(14.,31.);
      vec_h_isodirectphoton_pt_central.at(i)->Draw("");
    }
    else vec_h_isodirectphoton_pt_central.at(i)->Draw("same");
  }
  histMeanBin_isodirectphoton_pt_central->SetLineWidth(3);
  histMeanBin_isodirectphoton_pt_central->SetLineColor(kRed+1);
  histMeanBin_isodirectphoton_pt_central->Draw("same");

  cSpectrum->cd(2);
  gPad->SetLogy(0);
  vector <TH1D*> vec_ratio_isodirectphoton_pt_central;
  for(int i = 0; i < vec_rootInFileName.size(); i++){
    vec_ratio_isodirectphoton_pt_central.push_back((TH1D*)vec_h_isodirectphoton_pt_central.at(i)->Clone("ratio_isodirectphoton_pt_central"));
    vec_ratio_isodirectphoton_pt_central.at(i)->Divide(histMeanBin_isodirectphoton_pt_central);
    if( i == 0){
      vec_ratio_isodirectphoton_pt_central.at(i)->GetXaxis()->SetRangeUser(14.,31.);
      vec_ratio_isodirectphoton_pt_central.at(i)->GetYaxis()->SetRangeUser(0.01,30.);
      vec_ratio_isodirectphoton_pt_central.at(i)->Draw("");
    }
    else vec_ratio_isodirectphoton_pt_central.at(i)->Draw("same");
  }

  cSpectrum->cd(3);
  gPad->SetLogy(0);
  for(int i = 0; i < vec_rootInFileName.size(); i++){
    if (i==0){
      vec_h_dPhiJetGamma_central.at(i)->Draw("");
    }
    else vec_h_dPhiJetGamma_central.at(i)->Draw("same");
  }
  histMeanBin_dPhiJetGamma_central->SetLineWidth(3);
  histMeanBin_dPhiJetGamma_central->SetLineColor(kRed+1);
  histMeanBin_dPhiJetGamma_central->Draw("same");

  cSpectrum->cd(4);
  gPad->SetLogy(1);
  vector <TH1D*> vec_ratio_dPhiJetGamma_central;
  for(int i = 0; i < vec_rootInFileName.size(); i++){
    vec_ratio_dPhiJetGamma_central.push_back((TH1D*)vec_h_dPhiJetGamma_central.at(i)->Clone("ratio_dPhiJetGamma_central"));
    vec_ratio_dPhiJetGamma_central.at(i)->Divide(histMeanBin_dPhiJetGamma_central);
    if( i == 0){
      vec_ratio_dPhiJetGamma_central.at(i)->GetYaxis()->SetRangeUser(0.1,1000.);
      vec_ratio_dPhiJetGamma_central.at(i)->Draw("");
    }
    else vec_ratio_dPhiJetGamma_central.at(i)->Draw("same");
  }


  cSpectrum->cd(5);
  gPad->SetLogy(0);
  for(int i = 0; i < vec_rootInFileName.size(); i++){
    if (i==0){
      vec_h_xObs_pGoing_central.at(i)->GetXaxis()->SetRangeUser(14.,31.);
      vec_h_xObs_pGoing_central.at(i)->Draw("");
    }
    else vec_h_xObs_pGoing_central.at(i)->Draw("same");
  }
  histMeanBin_xObs_pGoing_central->SetLineWidth(3);
  histMeanBin_xObs_pGoing_central->SetLineColor(kRed+1);
  histMeanBin_xObs_pGoing_central->Draw("same");

  cSpectrum->cd(6);
  gPad->SetLogy(1);
  vector <TH1D*> vec_ratio_xObs_pGoing_central;
  for(int i = 0; i < vec_rootInFileName.size(); i++){
    vec_ratio_xObs_pGoing_central.push_back((TH1D*)vec_h_xObs_pGoing_central.at(i)->Clone("ratio_xObs_pGoing_central"));
    vec_ratio_xObs_pGoing_central.at(i)->Divide(histMeanBin_xObs_pGoing_central);
    if( i == 0){
      vec_ratio_xObs_pGoing_central.at(i)->GetYaxis()->SetRangeUser(0.1,1000.);
      vec_ratio_xObs_pGoing_central.at(i)->Draw("");
    }
    else vec_ratio_xObs_pGoing_central.at(i)->Draw("same");
  }



  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  for(int i = 0; i < vec_rootInFileName.size(); i++){
    double minimum_ratio_dPhiJetGamma_central = vec_ratio_dPhiJetGamma_central.at(i)->GetMinimum();
    double maximum_ratio_dPhiJetGamma_central = vec_ratio_dPhiJetGamma_central.at(i)->GetMaximum();
    if(minimum_ratio_dPhiJetGamma_central < 0.){
      vec_vetoName.push_back(vec_rootInFileName.at(i));
      printf("vetoed file %s with ratio_dPhiJetGamma minimum = %f\n",vec_rootInFileName.at(i).Data(), minimum_ratio_dPhiJetGamma_central);
      continue;
    }
    if(maximum_ratio_dPhiJetGamma_central > vetoFac){
      vec_vetoName.push_back(vec_rootInFileName.at(i));
      printf("vetoed file %s with ratio_dPhiJetGamma maximum = %f\n",vec_rootInFileName.at(i).Data(), maximum_ratio_dPhiJetGamma_central);
      continue;
    }
  }

  for(int i = 0; i < vec_rootInFileName.size(); i++){
    double minimum_ratio_xObs_pGoing_central = vec_ratio_xObs_pGoing_central.at(i)->GetMinimum();
    double maximum_ratio_xObs_pGoing_central = vec_ratio_xObs_pGoing_central.at(i)->GetMaximum();
    if(minimum_ratio_xObs_pGoing_central < 0.){
      vec_vetoName.push_back(vec_rootInFileName.at(i));
      printf("vetoed file %s with ratio_xObs_pGoing minimum = %f\n",vec_rootInFileName.at(i).Data(), minimum_ratio_xObs_pGoing_central);
      continue;
    }
    if(maximum_ratio_xObs_pGoing_central > vetoFac){
      vec_vetoName.push_back(vec_rootInFileName.at(i));
      printf("vetoed file %s with ratio_xObs_pGoing maximum = %f\n",vec_rootInFileName.at(i).Data(), maximum_ratio_xObs_pGoing_central);
      continue;
    }
  }

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  ofstream myfile;
  myfile.open(Form("do_hadd_vetoed.sh",dirName.Data()), ios::trunc | ios::out);
  myfile << Form("hadd -f %s",outFileName.Data());
  for(int i = 0; i < vec_rootInFileName.size(); i++){
    bool isVetoed = false;
    for(int j = 0; j < vec_vetoName.size(); j++){
      if(vec_rootInFileName.at(i) == vec_vetoName.at(j)){
        isVetoed = true;
      }
    }
    if(!isVetoed) myfile << " " << vec_rootInFileName.at(i).Data();
  }
  myfile.close();

}
