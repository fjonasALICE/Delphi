//macro to add histogram files
//NOTE: This macro is kept for back compatibility only.
//Use instead the executable $ROOTSYS/bin/hadd
//
//This macro will add histograms from a list of root files and write them
//to a target root file. The target file is newly created and must not be
//identical to one of the source files.
//
//Author: Sven A. Schmidt, sven.schmidt@cern.ch
//Date:   13.2.2001

//This code is based on the hadd.C example by Rene Brun and Dirk Geppert,
//which had a problem with directories more than one level deep.
//(see macro hadd_old.C for this previous implementation).
//
//The macro from Sven has been enhanced by
//   Anne-Sylvie Nicollerat <Anne-Sylvie.Nicollerat@cern.ch>
// to automatically add Trees (via a chain of trees).
//
//To use this macro, modify the file names in function hadd.
//
//NB: This macro is provided as a tutorial.
//    Use $ROOTSYS/bin/hadd to merge many histogram files

// uses kisaverage-Bit for histograms!!!
#include <string>
#include <sstream>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include <vector>
#include "Riostream.h"
#include <algorithm>
#include "TString.h"

TList *FileList;
TFile *Target;

void MergeRootfile( TDirectory *target, TList *sourcelist );
double CalcMedian( vector<double> vec );

int main(int argc, char **argv) {
  Target = TFile::Open( "result.root", "RECREATE" );

  FileList = new TList();
  for(int i=1; i<=argc; i++){
    FileList->Add(TFile::Open(argv[i]));
  }

  MergeRootfile( Target, FileList );

  return 0;
}

void MergeRootfile( TDirectory *target, TList *sourcelist ) {


  TH1D *h1_dist_bin1 = new TH1D("h1_dist_bin1","cross section in bin 1 per subjob (i.e. per 1M events)", 2e3, 1.e-2, 10.);
  TH1D *h1_dist_bin2 = new TH1D("h1_dist_bin2","cross section in bin 2 per subjob (i.e. per 1M events)", 2e3, 0.7e-2, 5.);
  TH1D *h1_dist_bin3 = new TH1D("h1_dist_bin3","cross section in bin 3 per subjob (i.e. per 1M events)", 2e3, 0.3e-2, 3.0);
  TH1D *h1_dist_bin4 = new TH1D("h1_dist_bin4","cross section in bin 4 per subjob (i.e. per 1M events)", 2e3, 0.1e-2, 1.5);
  TH1D *h1_dist_bin5 = new TH1D("h1_dist_bin5","cross section in bin 5 per subjob (i.e. per 1M events)", 2e3, 0.04e-2, 1.);

  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );

  TFile *first_source = (TFile*)sourcelist->First();
  first_source->cd( path );
  TDirectory *current_sourcedir = gDirectory;
  //gain time, do not add the objects in the list in memory
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  // loop over all keys in this directory
  TChain *globChain = 0;
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey())) {

    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

    // read object from first source file
    first_source->cd( path );
    TObject *obj = key->ReadObj();

    if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
      // descendant of TH1 -> merge it

      //      cout << "Merging histogram " << obj->GetName() << endl;
      TH1 *h1 = (TH1*)obj;
      vector<TH1*> histvec;
      int nHistos = 0;
      // loop over all source files and add the content of the
      // correspondant histogram to the one pointed to by "h1"
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      while ( nextsource ) {

        // make sure we are at the correct directory level by cd'ing to path
        nextsource->cd( path );
        TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
        if (key2) {
          TH1* h2 = (TH1*)key2->ReadObj();
          histvec.push_back(h2);
          //delete h2;
          ++nHistos;

        }
        nextsource = (TFile*)sourcelist->After( nextsource );
      }

      histvec.push_back(h1); // h1 auch noch mit reinnehmen
      TString tempString = h1->GetName();
      // if( !(tempString.Contains("atlas2017") ||
      // 	    tempString.Contains("atlas2013") ||
      // 	    tempString.Contains("atlas2016")) ) continue;
      //if( !tempString.Contains("rap1") ) continue;

      vector<double> tempVec;
      for(int i = 1; i <= h1->GetNbinsX(); i++){
        for(int j = 0; j <= nHistos; j++){
          tempVec.push_back(histvec.at(j)->GetBinContent(i));
          // if( !tempString.Contains("ptSim") ) continue;
          // if(i == 1) h1_dist_bin1->Fill(histvec.at(j)->GetBinContent(i));
          // if(i == 2) h1_dist_bin2->Fill(histvec.at(j)->GetBinContent(i));
          // if(i == 3) h1_dist_bin3->Fill(histvec.at(j)->GetBinContent(i));
          // if(i == 4) h1_dist_bin4->Fill(histvec.at(j)->GetBinContent(i));
          // if(i == 5) h1_dist_bin5->Fill(histvec.at(j)->GetBinContent(i));
        }
        h1->SetBinContent(i, CalcMedian(tempVec));
        tempVec.clear();
      }
      TFile* f1 = new TFile("dist.root","RECREATE");
      h1_dist_bin1->Write();
      h1_dist_bin2->Write();
      h1_dist_bin3->Write();
      h1_dist_bin4->Write();
      h1_dist_bin5->Write();
      f1->Write();
      f1->Close();

      std::cout << h1->GetName() << "\n" << endl;
      // for(int j = 1; j <= histvec.at(0)->GetNbinsX(); j++){ // bin loop
      //   vector<double> tempVal= {0.};
      //   bool veto = false;

      //   for(int m = 0; m < nHistos; ++m){ //histo loop1
      //     tempVal.push_back(histvec.at(m)->GetBinContent(j));
      //   }

      //   *std::max_element(tempVal.rbegin(), tempVal.rend()) = -1.;

      //   for(int m = 0; m < nHistos; ++m){ //histo loop2
      //     if(tempVal.at(m) > 0.) h1->SetBinContent(j, h1->GetBinContent(j) + histvec.at(m)->GetBinContent(j));
      //   }
      // }
      // h1->Scale(1./(nHistos));

    }
    else if ( obj->IsA()->InheritsFrom( TTree::Class() ) ) {

      // loop over all source files create a chain of Trees "globChain"
      const char* obj_name= obj->GetName();

      globChain = new TChain(obj_name);
      globChain->Add(first_source->GetName());
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      //      const char* file_name = nextsource->GetName();
      // cout << "file name  " << file_name << endl;
      while ( nextsource ) {

        globChain->Add(nextsource->GetName());
        nextsource = (TFile*)sourcelist->After( nextsource );
      }

    } else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
      // it's a subdirectory

      cout << "Found subdirectory " << obj->GetName() << endl;

      // create a new subdir of same name and title in the target file
      target->cd();
      TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

      // newdir is now the starting point of another round of merging
      // newdir still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion
      MergeRootfile( newdir, sourcelist );

    } else {

      // object is of no type that we know or can handle
      cout << "Unknown object type, name: "
           << obj->GetName() << " title: " << obj->GetTitle() << endl;
    }

    // now write the merged histogram (which is "in" obj) to the target file
    // note that this will just store obj in the current directory level,
    // which is not persistent until the complete directory itself is stored
    // by "target->Write()" below
    if ( obj ) {
      target->cd();

      //!!if the object is a tree, it is stored in globChain...
      if(obj->IsA()->InheritsFrom( TTree::Class() ))
        globChain->Merge(target->GetFile(),0,"keep");
      else
        obj->Write( key->GetName() );
    }

  } // while ( ( TKey *key = (TKey*)nextkey() ) )


  // save modifications to target file
  target->SaveSelf(kTRUE);
  TH1::AddDirectory(status);

}

double CalcMedian(vector<double> vec)
{
  double median;
  size_t size = vec.size();

  sort(vec.begin(), vec.end());

  if (size  % 2 == 0)
  {
      median = (vec[size / 2 - 1] + vec[size / 2]) / 2;
  }
  else 
  {
      median = vec[size / 2];
  }

  return median;
}
