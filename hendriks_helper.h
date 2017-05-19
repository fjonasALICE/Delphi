// helper functions that make standard pythia stuff

void Set_Pythia_Randomseed(Pythia8::Pythia &p);
void Fill_Non_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h);
void Fill_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h);
void Pass_Parameters_To_Pythia(Pythia8::Pythia &p, int argc, char **argv);
void SoftQCD_HardQCD_Switch(int iBin, double *pTHatBin, char **argv, Pythia8::Pythia &p, int &nEvent);
void Add_Histos_Scale_Write2File( std::vector <TH1D*> &vec_temp_histo, TH1* final_histo, TFile &file, double etaRange);


//----------------------------------------------------------------------
void Fill_Non_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() < 90) h->Fill(event[i].pT());
    }
  }
}

//----------------------------------------------------------------------
void Fill_Decay_Photon_Pt(Pythia8::Event &event, float etaMax, TH1 *h){
  for (int i = 5; i < event.size(); i++) {
    if(event[i].isFinal() && event[i].id() == 22 && TMath::Abs(event[i].eta()) < etaMax ) {
      if(event[i].status() > 90) h->Fill(event[i].pT());
    }
  }
}

//----------------------------------------------------------------------
void Set_Pythia_Randomseed(Pythia8::Pythia &p){

  TRandom3 *rand = new TRandom3(0);
  p.readString("Random:setSeed = on");
  p.readString(Form("Random:seed = %ld", (long)(rand->Rndm()*9e8) ));
  delete rand;
  return;
}

//----------------------------------------------------------------------
void Pass_Parameters_To_Pythia(Pythia8::Pythia &p, int argc, char **argv){

  p.readString(Form("Beams:eCM = %f", strtof(argv[4], NULL) ));
  p.readString(Form("SigmaProcess:renormMultFac = %f", strtof(argv[5],NULL) ));
  p.readString(Form("SigmaProcess:factorMultFac = %f", strtof(argv[6],NULL) ));

  if(argc > 7){
    if( !strcmp(argv[7],"noHadro") ){
      p.readString("HadronLevel:all = off");
    }
    else if( !strcmp(argv[7],"noMPI") ){
      p.readString("PartonLevel:MPI = off");
    }
    else if( !strcmp(argv[7],"noMPInoHadro") ){
      p.readString("PartonLevel:MPI = off");
      p.readString("HadronLevel:all = off");
    }
    else if( !strcmp(argv[7],"noShower") ){
      p.readString("Check:event = off");
      p.readString("PartonLevel:MPI = off");
      p.readString("PartonLevel:FSR = off");
      p.readString("PartonLevel:ISR = off");
      p.readString("PartonLevel:Remnants = off");
    }
  }
  else printf("---no optional argument is given -> full pythia events will be generated---\n");

  return;
}

//----------------------------------------------------------------------
void SoftQCD_HardQCD_Switch(int iBin, double *pTHatBin, char **argv, Pythia8::Pythia &p, int &nEvent){

  if( !strcmp(argv[2],"frag") ||
      !strcmp(argv[2],"decay") ){

    // SoftQCD only in the first pthat bin, HardQCD otherwise
    if (iBin == 0) {
      p.readString("HardQCD:all = off");
      p.readString("SoftQCD:inelastic = on");
      nEvent *= 2;
    } else {
      if(iBin == 1)
	nEvent /= 2;
      p.readString("HardQCD:all = on");
      p.readString("SoftQCD:inelastic = off");
    }

  }
  else if( !strcmp(argv[2],"dir") ){
    p.readString("PromptPhoton:qg2qgamma = on");
    p.readString("PromptPhoton:qqbar2ggamma = on");
  }
  else
    printf("No process switched on! Give \"dir\" or \"frag\" or \"decay\" as second argument");

  p.settings.parm("PhaseSpace:pTHatMin", pTHatBin[iBin]);
  p.settings.parm("PhaseSpace:pTHatMax", pTHatBin[iBin+1]);


  return;
}

//----------------------------------------------------------------------
void Add_Histos_Scale_Write2File( std::vector <TH1D*>& vec, TH1* final_histo, TFile &file, double etaRange){

  file.cd();

  for( int i = 0; i < vec.size(); i++){
    final_histo->Add(vec.at(i));
    vec.at(i)->Scale(1./etaRange, "width");
    vec.at(i)->SetXTitle("p_{T} (GeV/#it{c})");
    vec.at(i)->SetYTitle("#frac{d#sigma}{dp_{T}d#eta}");
    vec.at(i)->Write();
  }

  final_histo->SetXTitle("p_{T} (GeV/#it{c})");
  final_histo->SetYTitle("#frac{d#sigma}{dp_{T}d#eta}");
  final_histo->Scale(1./etaRange, "width");
  final_histo->Write();


  return;
}
