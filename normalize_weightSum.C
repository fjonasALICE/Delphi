void normalize_weightSum(const char* rootInFileName){

  //  char rootOutFileName[1024];
  TString rootOutFileName = TString(rootInFileName);
  rootOutFileName.ReplaceAll(".root","");
  rootOutFileName.Append("_normalized.root");
  //  snprintf( rootOutFileName, sizeof(rootOutFileName), "%s_normalized.root", rootInFileName);

  // for MB production only first bin is considered (i.e. no pthat bins used)
  bool isMB = false;
  if( rootOutFileName.Contains("MB") )
    isMB = true;

  TFile *infile = new TFile(rootInFileName);
  TFile *target = new TFile(rootOutFileName,"RECREATE");
  //  infile->ls();

  cout << "Target path: " << target->GetPath() << endl;
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );

  infile->cd();
  TDirectory *current_sourcedir2 = gDirectory;

  // loop over all keys in this directory
  TChain *globChain = 0;
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TIter nextkey2( current_sourcedir2->GetListOfKeys() );
  TKey *key2, *oldkey2=0;

  vector <TH1*> vec_histos_pthat_bins;
  // scale pTHat wise bins

  while ( (key2 = (TKey*)nextkey2())) {

    //keep only the highest cycle number for each key
    if (oldkey2 && !strcmp(oldkey2->GetName(),key2->GetName())) continue;

    // read object from first source file
    infile->cd();
    TObject *obj2 = key2->ReadObj();
    if ( obj2->IsA()->InheritsFrom( TDirectory::Class() ) ){
      TDirectory *dir = (TDirectory*)obj2;
      dir->cd();
    }

  TKey *key, *oldkey=0;
  TDirectory *current_sourcedir = gDirectory;
  //  gDirectory->ls();
  TIter nextkey( current_sourcedir->GetListOfKeys() );

  while ( (key = (TKey*)nextkey())) {

    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

    // read object from first source file
    infile->cd();
    //infile->ls();
    TObject *obj = key->ReadObj();
    if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
      TH1 *h1 = (TH1*)obj;
      TString tempString = h1->GetName();
      if(tempString.Contains("weightSum"))
        continue;
      // if(tempString.Contains("pTHat"))
      //   continue;

      if(tempString.Contains("bin")){
	
	if(tempString.Contains("bin_00")){
	  vec_histos_pthat_bins.clear();
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_00");
	  if(h_weightSum_bin->GetBinContent(1) > 0.000001)
	    h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_01")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_01");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_02")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_02");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_03")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_03");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_04")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_04");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_05")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_05");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_06")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_06");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_07")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_07");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_08")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_08");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_09")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_09");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_10")) && !isMB ){
	    TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_10");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_11")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_11");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_12")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_12");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_13")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_13");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_14")) && !isMB ) {
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_14");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_15")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_15");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_16")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_16");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_17")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_17");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_18")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_18");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}
	if( (tempString.Contains("bin_19")) && !isMB ){
	  TH1 *h_weightSum_bin = (TH1*)infile->Get("h_weightSum_bin_19");
	  h1->Scale(1./h_weightSum_bin->GetBinContent(1));
	  vec_histos_pthat_bins.push_back(h1);
	}


	target->cd();
	//h1->Write();
      }

      // buggy for pthatbin histos
      if(!tempString.Contains("bin")){
        target->cd();
        TH1 *final_histo = (TH1*)h1->Clone();
        final_histo->Reset();
        for( int i = 0; i < vec_histos_pthat_bins.size(); i++){
          final_histo->Add(vec_histos_pthat_bins.at(i));
	  //          vec_histos_pthat_bins.at(i)->Write();
	  if ( isMB ) break;
        }
        final_histo->Write();
      }

    }

  }
  }




  infile->Close();

  return;
}
