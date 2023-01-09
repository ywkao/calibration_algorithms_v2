#ifndef __HIST_COLLECTION_H__
#define __HIST_COLLECTION_H__

#include <vector>
#include <iostream>
#include <sstream>

#include <TColor.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TString.h>
#include <TStyle.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//class HistCollection : public edm::one::EDAnalyzer<edm::one::SharedResources>
class HistCollection
{
    public:
        HistCollection();
        ~HistCollection();

        void init(TString varName, edm::Service<TFileService> fs);
        void set_string(std::ostringstream &hnamestr, TString token, int idx);
        void fill_histogram(int channleId, double value);
        void fill_hist2D(int channleId, double value, double cm);
        void write_histograms();

        std::vector<TH2D*>     get_vecotor_th2d() {return mH2;};
        std::vector<TProfile*> get_vecotor_tprofile() {return mPr;};
        
    private:
        TString mVarName;
        std::vector<TH1D*> mH1;
        std::vector<TH2D*> mH2;
        std::vector<TProfile*> mPr;
};

HistCollection::HistCollection(){}

HistCollection::~HistCollection()
{
    //printf("end of %s HistCollection\n", mVarName.Data());

    for(int i=0; i<234; ++i) {
        mH1[i]->Delete();
        mH2[i]->Delete();
        mPr[i]->Delete();
    }
}

void HistCollection::init(TString varName, edm::Service<TFileService> fs)
{
    mVarName = varName;

    std::ostringstream histName (std::ostringstream::ate);
    for(int i=0; i<234; ++i) {
        set_string(histName, "h_"+mVarName+"_channel_", i);
        TH1D *m1 = fs->make<TH1D>(histName.str().c_str(), histName.str().c_str(), 175, -25, 150);

        set_string(histName, "h2d_"+mVarName+"_channel_", i);
        TH2D *m2 = fs->make<TH2D>(histName.str().c_str(), ";CM #minus CM_{pedestal};ADC #minus ADC_{pedestal}", 19, -9.5, 9.5, 39, -9.5, 29.5);
        m2 -> SetTitle("");

        set_string(histName, "p2d_"+mVarName+"_channel_", i);
        TProfile *mP = fs->make<TProfile>(histName.str().c_str(), ";CM #minus CM_{pedestal};ADC #minus ADC_{pedestal}", 19, -9.5, 9.5, "S");
        mP -> SetMarkerStyle(20);
        mP -> SetMarkerSize(0.5);
        mP -> SetLineWidth(2);
        mP -> SetStats(0);

        mH1.push_back(m1);
        mH2.push_back(m2);
        mPr.push_back(mP);

        //printf(">>> HistCollection::HistCollection: %s\n", histName.str().c_str());
    }
}

void HistCollection::fill_histogram(int channleId, double value)
{
    if(channleId>233) {
        printf("[ERROR] HistCollection::fill_histogram : channleId out of expectation.\n");
        return;
    }

    mH1[channleId]->Fill(value);
}

void HistCollection::fill_hist2D(int channleId, double value, double cm)
{
    if(channleId>233) {
        printf("[ERROR] HistCollection::fill_hist2D : channleId out of expectation.\n");
        return;
    }

    mH2[channleId]->Fill(cm, value);
    mPr[channleId]->Fill(cm, value);
}

void HistCollection::write_histograms()
{
    for(int i=0; i<234; ++i) { mH1[i]->Write(); }
    for(int i=0; i<234; ++i) { mH2[i]->Write(); }
    for(int i=0; i<234; ++i) { mPr[i]->Write(); }
}

void HistCollection::set_string(std::ostringstream &hnamestr, TString token, int idx)
{
    hnamestr.clear();
    hnamestr.str(token.Data());
    hnamestr<<idx;
}

#endif
