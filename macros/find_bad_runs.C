#include <vector>

void find_bad_runs( const char* dirName1 = "abc/",
		    const char* dirName2 = "abc/",
		    const char* dirName3 = "abc/",
		    const char* dirName4 = "abc/",
                    int firstFileNumber1 = 1, int lastFileNumber1 = 400,
                    int firstFileNumber2 = 1, int lastFileNumber2 = 400,
                    int firstFileNumber3 = 1, int lastFileNumber3 = 400,
                    int firstFileNumber4 = 1, int lastFileNumber4 = 400,
                    const char* rootDirName = "fluctCut_test",
                    const char* histoName = "h1_gamma_pt1_fluctCut30",
                    const char* preString = "",
                    const char* sufString1 = "100000",
                    const char* sufString2 = "100000",
                    const char* sufString3 = "200000",
                    const char* sufString4 = "500000"
                    ){

  vector <TFile*> vec_files;
  for( int i = firstFileNumber1; i <= lastFileNumber1; i++){
    vec_files.push_back(TFile::Open(Form("%s%s%d%s.root",dirName1,preString,i,sufString1)));
  }
  for( int i = firstFileNumber2; i <= lastFileNumber2; i++){
    vec_files.push_back(TFile::Open(Form("%s%s%d%s.root",dirName2,preString,i,sufString2)));
  }
  for( int i = firstFileNumber3; i <= lastFileNumber3; i++){
    vec_files.push_back(TFile::Open(Form("%s%s%d%s.root",dirName3,preString,i,sufString3)));
  }
  for( int i = firstFileNumber4; i <= lastFileNumber4; i++){
    vec_files.push_back(TFile::Open(Form("%s%s%d%s.root",dirName4,preString,i,sufString4)));
  }

  vec_files.at(0)->ls();

  TDirectory *dir;
  vector <TH1D*> vec_h1;
  for (int i = 0; i < vec_files.size(); i++){
    //    dir = (TDirectory*)vec_files.at(i)->Get(rootDirName);
    //    dir->cd();
    vec_files.at(i)->cd();
    //if( i == 0 ) dir->ls();
    vec_h1.push_back( (TH1D*)gDirectory->Get(histoName) );
    cout << "file read in: " << i << endl;
    vec_h1.at(i)->SetName( Form("%d",i+firstFileNumber1) );
  }


  cout << "--------------------------v" << endl;
  cout << "bad runs: " << endl;
  vector <int> sorted_bad_runs1;
  vector <int> sorted_bad_runs2;

  for( int iBin = 1; iBin < vec_h1.at(0)->GetNbinsX(); iBin++){
    double mean = 0.;

    // first check for large weight
    for( int i = 0; i < vec_h1.size(); i++){
      if( vec_h1.at(i)->GetBinContent(iBin) > 0. ){ // sometimes there are negative entries, ignore them
        mean += vec_h1.at(i)->GetBinContent(iBin);
      }
      //      cout << vec_h1.at(i)->GetBinContent(iBin) << endl;
    }
    
    mean /= vec_h1.size();

    for( int i = 0; i < vec_h1.size(); i++){
      if( vec_h1.at(i)->GetBinContent(iBin) > mean*10.){
        // cout << "bad run found: " ;
        // cout << i+firstFileNumber ;
        // cout << " in bin number " ;
        // cout << iBin ;
        // cout << " with mean = " << mean << " and large value = " << vec_h1.at(i)->GetBinContent(iBin) << endl;
        vec_h1.at(i)->SetLineColor(kRed);
        sorted_bad_runs1.push_back(i+firstFileNumber1);
      }
    }

    mean = 0.;
    // second check for large weight with corrected mean
    for( int i = 0; i < vec_h1.size(); i++){
      bool isBadRun = false;
      for(int j = 0; j < sorted_bad_runs1.size(); j++){
        //        cout << i+firstFileNumber1 << " \t " << sorted_bad_runs1.at(j) << endl;
        if( i+firstFileNumber1 == sorted_bad_runs1.at(j) && !isBadRun){
          isBadRun = true;
          break;
        }
      }

      if (!isBadRun){
        //        cout << i+firstFileNumber1 << " is a good run." << endl;
        if( vec_h1.at(i)->GetBinContent(iBin) > 0. ){ // sometimes there are negative entries, ignore them
          mean += vec_h1.at(i)->GetBinContent(iBin);
        }
      }
    }

    mean /= vec_h1.size();

    //cout << "\nStart second iteration with corrected mean" << endl;
    for( int i = 0; i < vec_h1.size(); i++){
      if( vec_h1.at(i)->GetBinContent(iBin) > mean*10. ){
        // cout << "bad run found: " ;
        // cout << i+firstFileNumber1 ;
        // cout << " in bin number " ;
        // cout << iBin ;
        // cout << " with mean = " << mean << " and large value = " << vec_h1.at(i)->GetBinContent(iBin) << endl;
        vec_h1.at(i)->SetLineColor(kRed);
        sorted_bad_runs2.push_back(i+firstFileNumber1);
      }
    }


  } // end of histogram bin loop

  std::sort(sorted_bad_runs2.begin(), sorted_bad_runs2.end());
  // remove double entries
  auto last = std::unique(sorted_bad_runs2.begin(), sorted_bad_runs2.end());
  sorted_bad_runs2.erase(last, sorted_bad_runs2.end());


  // cout << "sorted bad runs after second iteration:" << endl;
  for(int i = 0; i < sorted_bad_runs2.size(); i++)
    cout <<  sorted_bad_runs2.at(i) << endl;


  // for(int i = 0; i < sorted_bad_runs2.size()-1; i++){
  //   if( i == 0 && sorted_bad_runs2.at(0) != firstFileNumber1 )
  //     cout << Form("{{%d..%d},", firstFileNumber1, sorted_bad_runs2.at(i)-1);
  //   if( sorted_bad_runs2.at(i)+1 <= sorted_bad_runs2.at(i+1)-1 ) // solves the problem of two consecutive following runs
  //     cout << Form("{%d..%d},", sorted_bad_runs2.at(i)+1 , sorted_bad_runs2.at(i+1)-1);
  //   if( i == sorted_bad_runs2.size()-2)
  //     cout << Form("{%d..%d}}.root", sorted_bad_runs2.at(i+1)+1, lastFileNumber) << endl;
  // }

  cout << "--------------------------" << endl << endl << endl;

  TCanvas *c = new TCanvas("c","c", 2000, 1000);
  c->cd();
  vec_h1.at(0)->Draw("");
  for( int i = 0; i < vec_h1.size(); i++){
    if( i == 0 ) vec_h1.at(i)->Draw("");
    else vec_h1.at(i)->Draw("same");
  }



   return;

}
