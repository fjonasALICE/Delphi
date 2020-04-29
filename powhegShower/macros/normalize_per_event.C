void normalize_per_event(const char* rootInFileName, bool chooseGammaJetCorr=kFALSE){

  //  char rootOutFileName[1024];
  TString rootOutFileName = TString(rootInFileName);
  rootOutFileName.ReplaceAll(".root","");
  if(chooseGammaJetCorr)
    rootOutFileName.Append("_normalized_GJcorr.root");
  else
    rootOutFileName.Append("_normalized_spectra.root");
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
  TIter nextkey( current_sourcedir2->GetListOfKeys() );
  TKey *key, *oldkey=0;

  TH1D *h_nEvents = 0x0;
  while ( (key = (TKey*)nextkey())) {

    printf("key->GetName() = %s\n", key->GetName());
    
    // read object from first source file
    infile->cd();
    //infile->ls();
    TObject *obj = key->ReadObj();
    if ( obj->IsA()->InheritsFrom( TH1::Class() ) ){
      TH1 *h1 = (TH1*)obj;
      TString tempString = h1->GetName();
      if(tempString.Contains("h_nEvents")){ // only works when h_nEvents is first histo
	h_nEvents = (TH1D*)obj;
      }else{
	h1->Scale(1./h_nEvents->GetBinContent(1));
	target->cd();
	h1->Write();
      }
    }
  }  
  



  infile->Close();
  
  return;
}
