void delete_readmes_in_rootfile(TString rootFileName){

  TFile *f = new TFile(rootFileName,"UPDATE");
  //  f->ls();

  TDirectory *current_sourcedir = gDirectory;
  TH1::AddDirectory(kFALSE);
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key;
  
  int wasAlready = 0;
  while ( (key = (TKey*)nextkey())) {
    TString keyName = key->GetName();
    printf("key->GetName() = %s\n", keyName.Data());
    if(keyName.Contains("README"))
      if(wasAlready++){
	printf("int wasAlready = %d\n", wasAlready);
	printf("to be deleted: %s\n", keyName.Data());
	key->Delete();
      }
  }

  f->Close();
  return;
}
