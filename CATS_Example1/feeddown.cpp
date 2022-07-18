R__ADD_LIBRARY_PATH("/home/emma/CATS/lib/");
R__LOAD_LIBRARY(libCATSbasic.so);
R__LOAD_LIBRARY(libCATSdev.so);
R__LOAD_LIBRARY(libCATSextended.so);

#include "/home/emma/CATS/include/CATS.h"
#include "/home/emma/CATS/include/CATStools.h"
#include "/home/emma/CATS/include/CATSconstants.h"
#include "/home/emma/CATS/include/DLM_CkModels.h"
#include "/home/emma/CATS/include/DLM_Potentials.h"
#include "/home/emma/CATS/include/DLM_Source.h"
#include "/home/emma/CATS/include/DLM_DecayMatrix.h"
#include "/home/emma/CATS/include/DLM_Ck.h"
#include "/home/emma/CATS/include/DLM_CkDecomposition.h"



//this potential is used to model proton-Lambda interaction

void FeedDown_KDstar_KD(const int& SEED, const int& NumIter){
    DLM_DecayMatrix DECM;
    float Mass_K=493.677;
    float Mass_D=1869.65;
    float Mass_pi0=134.9768;
    float Mass_Dstar=2010.26;

    DECM.SetFileName(TString::Format("/home/emma/CharmingAnalyses/DKDpi/oton_selections/Feeddown/Kdstar_KD_RS%i.root",SEED));
    DECM.SetHistoName("KDstar_KD");
    DECM.SetBins(1024,0,1024);
    DECM.SetNumDaughters1(1);
    DECM.SetMotherMass1(Mass_K);
    DECM.SetNumDaughters2(2);
    DECM.SetDaughterMass2(0,Mass_D);
    DECM.SetDaughterMass2(1,Mass_pi0);
    DECM.SetMotherMass2(Mass_Dstar);
    DECM.SetMeanMomentum(0);
    DECM.SetMomentumSpread(350);
    DECM.Run(SEED,NumIter);
}

TH2F* GetMatix(const char* rootfile, const char* filename) {
//   const char* name =
//       Form("%s/SmearSideband.root", path);
  auto file = TFile::Open(rootfile);
  file->ls();
  TH2F* matrix = (TH2F*)file->FindObjectAny(Form("%s", filename));
 // file->TFile::Close();
  return matrix;
}

void FillCkGraph(DLM_Ck *ck, DLM_CkDecomposition &ckGraph, TGraph *gr) {
        std::cout<<ck->GetNbins()<<" wowowoow"<<std::endl;

  for (unsigned int i = 0; i < ck->GetNbins(); ++i) {
    const float mom = ck->GetBinCenter(0, i);
    std::cout<<mom<<" "<<ckGraph.EvalCk(mom)<<std::endl;
    gr->SetPoint(i, mom, ckGraph.EvalCk(mom));
  }
}

void FeedDown_CF_KDstar_KD_plus(const int& SEED){
    CATS Kitty;
    unsigned NumMomBins = 128;
    double kMin = 0;
    double kMax = 1024;

    double mass1 = 493.677; //kaon
    double mass2 = 2010.26; //dstar

    Kitty.SetMomBins(NumMomBins,kMin,kMax);
    
    // definition of the source:
    CATSparameters source_func(CATSparameters::tSource, 1, true);
    Kitty.SetAnaSource(GaussSource, source_func);
    Kitty.SetAnaSource(0, 0.8115);//in fm the source size
    Kitty.SetQ1Q2(1);//1 for same charge, -1 for opposite charged particle
    Kitty.SetQuantumStatistics(0);
    Kitty.SetRedMass(mass1 * mass2 /(mass1 + mass2));
    Kitty.SetNumChannels(1);
    Kitty.SetChannelWeight(0, 1.);

    Kitty.KillTheCat();

    //GetMomentum gives us the center of the bin, GetCorrFun gives us the corresponding value of C(k)
    printf("1) Result for pp:\n");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        printf("   C(%6.2f)=%.2f  ", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
        if( (uBin+1)%5==0 && uBin!=(NumMomBins-1) ) printf("\n");
    }
    printf("\n\n");

    TH2F* histDecayKindematicsDstar=GetMatix(TString::Format("/home/emma/CharmingAnalyses/DKDpi/oton_selections/Feeddown/matrix/Kdstar_KD_RS%i.root",SEED), "KDstar_KD");

    auto DLM_KDstarplus = new DLM_Ck(1, 0, Kitty);

    DLM_CkDecomposition CkDec_KDstarplus_Smeared("KDstarplus_smear", 0,
                                               *DLM_KDstarplus,
                                               histDecayKindematicsDstar);
    CkDec_KDstarplus_Smeared.Update();
    
    auto grDstarminusSmeared = new TGraph();
    grDstarminusSmeared->SetLineWidth(2);
    grDstarminusSmeared->SetLineColor(2);

    grDstarminusSmeared->SetLineStyle(3);
    grDstarminusSmeared->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");

    FillCkGraph(DLM_KDstarplus, CkDec_KDstarplus_Smeared, grDstarminusSmeared);

    TGraph grGenuine;
    grGenuine.SetName("genuineCF");

    for (int uBin = 0; uBin < NumMomBins; ++uBin) {
        grGenuine.SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));  // genuine p-Phi CF with the parameters obtained in the fit
    }

    auto outputFile = new TFile("/home/emma/CharmingAnalyses/DKDpi/oton_selections/Feeddown/Dstar_KplusDplusOutput.root", "RECREATE");
    outputFile->cd();
    grDstarminusSmeared->Write("smeared");
    grGenuine.Write();
    outputFile->Close();

}

int feeddown(){

    int seed=7;

    FeedDown_KDstar_KD(seed,100000000);
    FeedDown_CF_KDstar_KD_plus(seed);
   
    return 0;
}

