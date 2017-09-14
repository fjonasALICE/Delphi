
void normalize_weightSum(const char* rootInFileName){

  char rootOutFileName[1024];
  snprintf( rootOutFileName, sizeof(rootOutFileName), "%s_normalized.root", rootInFileName);

  TFile *infile = new TFile(rootInFileName);
  TFile *target = new TFile(rootOutFileName,"RECREATE");
  infile->ls();

  cout << "Target path: " << target->GetPath() << endl;
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );

  infile->cd();
  TDirectory *current_sourcedir = gDirectory;
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  // loop over all keys in this directory
  TChain *globChain = 0;
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;

  vector <TH1*> vec_histos_pthat_bins;
  // scale pTHat wise bins
  while ( (key = (TKey*)nextkey())) {

    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

    // read object from first source file
    infile->cd();
    TObject *obj = key->ReadObj();

    if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
      TH1 *h1 = (TH1*)obj;
      TString tempString = h1->GetName();
      if(tempString.Contains("weightSum"))
        continue;


      if(tempString.Contains("bin")){

	if(tempString.Contains("bin_00")){
	  vec_histos_pthat_bins.clear();
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_00");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_01"))){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_01");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if(tempString.Contains("bin_02")){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_02");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_03"))){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_03");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if(tempString.Contains("bin_04")){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_04");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_05"))){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_05");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if(tempString.Contains("bin_06")){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_06");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_07"))){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_07");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if(tempString.Contains("bin_08")){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_08");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_09"))){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_09");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if(tempString.Contains("bin_10")){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_10");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_11"))){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_11");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if(tempString.Contains("bin_12")){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_12");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_13"))){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_13");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if(tempString.Contains("bin_14")){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_14");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_15"))){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_15");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}


	target->cd();
	//h1->Write();
      }

      if(!tempString.Contains("bin")){
        target->cd();
        TH1 *final_histo = (TH1*)h1->Clone();
        final_histo->Reset();
        for( int i = 0; i < vec_histos_pthat_bins.size(); i++){
          final_histo->Add(vec_histos_pthat_bins.at(i));
          vec_histos_pthat_bins.at(i)->Write();
        }
        final_histo->Write();
      }

    }

  }



  infile->Close();

  return;
}
