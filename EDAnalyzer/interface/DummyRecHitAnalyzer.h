#ifndef __DummyRecHitAnalyzer_H__
#define __DummyRecHitAnalyzer_H__

// -*- C++ -*-
//
// Package:    EDAnalyzer/DummyRecHitAnalyzer
// Class:      DummyRecHitAnalyzer
//
/**\class DummyRecHitAnalyzer DummyRecHitAnalyzer.cc EDAnalyzer/DummyRecHitAnalyzer/plugins/DummyRecHitAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Yu-Wei Kao
//         Created:  Wed, 02 Nov 2022 10:49:46 GMT
//
//

// system include files
#include <algorithm> // for removing space in string
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

//#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TH2D.h>
//#include <TLatex.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "calibration_algorithms/EDAnalyzer/interface/hgcalhit.h"
#include "calibration_algorithms/EDAnalyzer/interface/RunningCollection.h"

using namespace std;
using namespace edm;

//--------------------------------------------------------------------------------//
// It is not easy to employ TFileService outside a framework module...
//--------------------------------------------------------------------------------//
//#include "calibration_algorithms/EDAnalyzer/interface/HistCollection.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

//using reco::TrackCollection;

class DummyRecHitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit DummyRecHitAnalyzer(const edm::ParameterSet&);
        ~DummyRecHitAnalyzer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        void beginJob() override;
        void analyze(const edm::Event&, const edm::EventSetup&) override;
        void endJob() override;

        virtual void     Init(TTree *tree=0);
        virtual void     Init_my_output_info();

        virtual void     enable_pedestal_subtraction();
        virtual void     enable_cm_subtraction();

        virtual Long64_t LoadTree(Long64_t entry);
        virtual void     Load_metaData(); // pedestal, CM correlation, etc.

        virtual void     fill_histograms();
        virtual void     fill_profiles(int globalChannelId_, double adc_double_);

        virtual void     export_pedestals();
        virtual void     export_cm_parameters();

        virtual void     Show(Long64_t entry = -1);

        //--------------------------------------------------
        // variables initiated from python configure file
        //--------------------------------------------------
        std::vector<int> calibration_flags;
        TString myTag;

        //--------------------------------------------------
        // member data
        //--------------------------------------------------
        //edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file

        RunningCollection myRunStatCollection;
        RunningStatistics myRecorder;

        TString tag_calibration;
        TString tag_channelId;
        int globalChannelId;
        double adc_double;
        double adc_channel_CM;

        std::map<int, double> map_pedestals;
        std::map<int, std::vector<double> > map_cm_parameters;
        bool flag_perform_pedestal_subtraction;
        bool flag_perform_cm_subtraction;

        // output of calibrated RecHits
        TTree    * t_RecHit   ;
        
        // an instance of distributions
        TH1D     * h_adc      ;
        TH1D     * h_adcm     ;
        TH1D     * h_tot      ;
        TH1D     * h_toa      ;
        TH1D     * h_trigtime ;

        // summary of physical quantities
        TProfile * p_adc      ;
        TProfile * p_adcm     ;
        TProfile * p_tot      ;
        TProfile * p_toa      ;
        TProfile * p_trigtime ;
        TProfile * p_status   ;

        // summary for running statistics
        TH1D * h_correlation ;
        TH1D * h_slope       ;
        TH1D * h_intercept   ;

        //--------------------------------------------------
        // for reading tree
        //--------------------------------------------------
        TFile          *f1;
        TTree          *t1;

        TTree          *fChain;   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent; //!current Tree number in a TChain

        // Declaration of leaf types
        Int_t           event;
        Int_t           chip;
        Int_t           half;
        Int_t           channel;
        Int_t           adc;
        Int_t           adcm;
        Int_t           toa;
        Int_t           tot;
        Int_t           totflag;
        Int_t           trigtime;
        Int_t           trigwidth;
        Int_t           corruption;
        Int_t           bxcounter;
        Int_t           eventcounter;
        Int_t           orbitcounter;

        // List of branches
        TBranch        *b_event;
        TBranch        *b_chip;
        TBranch        *b_half;
        TBranch        *b_channel;
        TBranch        *b_adc;
        TBranch        *b_adcm;
        TBranch        *b_toa;
        TBranch        *b_tot;
        TBranch        *b_totflag;
        TBranch        *b_trigtime;
        TBranch        *b_trigwidth;
        TBranch        *b_corruption;
        TBranch        *b_bxcounter;
        TBranch        *b_eventcounter;
        TBranch        *b_orbitcounter;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
        edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DummyRecHitAnalyzer::DummyRecHitAnalyzer(const edm::ParameterSet& iConfig) {
    //: calibration_flags( iConfig.getParameter<std::vector<bool> >( "CalibrationFlags" ) ) {
    //: tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))) {

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif

    //now do what ever initialization is needed
    calibration_flags = iConfig.getParameter<vector<int> >( "CalibrationFlags" );
    myTag = iConfig.getParameter<string>( "DataType" );

    printf("[INFO] DataType: %s\n", myTag.Data());
    // raise error if the value of DataType is unexpected?
}

DummyRecHitAnalyzer::~DummyRecHitAnalyzer() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    //
    // please remove this method altogether if it would be left empty
    f1->Delete();
}

//
// member functions
//

void DummyRecHitAnalyzer::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    if (tree == 0) printf("[ERROR] something goes wrong with input tree\n");

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("event", &event, &b_event);
    fChain->SetBranchAddress("chip", &chip, &b_chip);
    fChain->SetBranchAddress("half", &half, &b_half);
    fChain->SetBranchAddress("channel", &channel, &b_channel);
    fChain->SetBranchAddress("adc", &adc, &b_adc);
    fChain->SetBranchAddress("adcm", &adcm, &b_adcm);
    fChain->SetBranchAddress("toa", &toa, &b_toa);
    fChain->SetBranchAddress("tot", &tot, &b_tot);
    fChain->SetBranchAddress("totflag", &totflag, &b_totflag);
    fChain->SetBranchAddress("trigtime", &trigtime, &b_trigtime);
    fChain->SetBranchAddress("trigwidth", &trigwidth, &b_trigwidth);
    fChain->SetBranchAddress("corruption", &corruption, &b_corruption);
    fChain->SetBranchAddress("bxcounter", &bxcounter, &b_bxcounter);
    fChain->SetBranchAddress("eventcounter", &eventcounter, &b_eventcounter);
    fChain->SetBranchAddress("orbitcounter", &orbitcounter, &b_orbitcounter);

    flag_perform_pedestal_subtraction = false;
    flag_perform_cm_subtraction = false;

    tag_calibration = "";
    tag_channelId = "_channel_7";
}

void DummyRecHitAnalyzer::Init_my_output_info()
{
    usesResource("TFileService");
    edm::Service<TFileService> fs; 

    myRunStatCollection = RunningCollection();

    //--------------------------------------------------
    // output information
    //--------------------------------------------------
    globalChannelId = -1;
    adc_double = 0.;
    adc_channel_CM = 0.;

    t_RecHit = fs->make<TTree>("t_RecHit","");
    t_RecHit -> Branch("event"        , &event        );
    t_RecHit -> Branch("chip"         , &chip         );
    t_RecHit -> Branch("half"         , &half         );
    t_RecHit -> Branch("channel"      , &channel      );
    t_RecHit -> Branch("adc"          , &adc_double   );
    t_RecHit -> Branch("adcm"         , &adcm         );

    // an instance of distributions
    h_adc       = fs->make<TH1D>("h_adc"      + tag_channelId , ";ADC;Entries"      , 175 , -25 , 150 );
    h_adcm      = fs->make<TH1D>("h_adcm"     + tag_channelId , ";ADC-1;Entries"    , 550 , -50 , 500 );
    h_tot       = fs->make<TH1D>("h_tot"      + tag_channelId , ";ToT;Entries"      , 100 , -2  , 2   );
    h_toa       = fs->make<TH1D>("h_toa"      + tag_channelId , ";ToA;Entries"      , 500 , 0   , 500 );
    h_trigtime  = fs->make<TH1D>("h_trigtime" + tag_channelId , ";trigtime;Entries" , 500 , 0   , 500 );

    // summary of physical quantities
    p_adc       = fs->make<TProfile>("p_adc"      , ";channel;ADC"      , 234, 0, 234, "S");
    p_adcm      = fs->make<TProfile>("p_adcm"     , ";channel;ADC-1"    , 234, 0, 234, "S");
    p_tot       = fs->make<TProfile>("p_tot"      , ";channel;ToT"      , 234, 0, 234, "S");
    p_toa       = fs->make<TProfile>("p_toa"      , ";channel;ToA"      , 234, 0, 234, "S");
    p_trigtime  = fs->make<TProfile>("p_trigtime" , ";channel;trigtime" , 234, 0, 234, "S");
    p_status    = fs->make<TProfile>("p_status"   , ";channel;status"   , 234, 0, 234, "S");

    p_adc      -> SetStats(0); p_adc      -> SetLineWidth(2);
    p_adcm     -> SetStats(0); p_adcm     -> SetLineWidth(2);
    p_tot      -> SetStats(0); p_tot      -> SetLineWidth(2);
    p_toa      -> SetStats(0); p_toa      -> SetLineWidth(2);
    p_trigtime -> SetStats(0); p_trigtime -> SetLineWidth(2);
    p_status   -> SetStats(0); p_status   -> SetLineWidth(2);

    // summary of running statistics
    h_correlation = fs->make<TH1D>("h_correlation" , ";channel;Correlation" , 234 , -0.5 , 233.5);
    h_slope       = fs->make<TH1D>("h_slope"       , ";channel;Slope"       , 234 , -0.5 , 233.5);
    h_intercept   = fs->make<TH1D>("h_intercept"   , ";channel;Intercept"   , 234 , -0.5 , 233.5);

    h_correlation -> SetStats(0); h_correlation -> SetMarkerStyle(20); h_correlation -> SetMarkerSize(0.5);
    h_slope       -> SetStats(0); h_slope       -> SetMarkerStyle(20); h_slope       -> SetMarkerSize(0.5);
    h_intercept   -> SetStats(0); h_intercept   -> SetMarkerStyle(20); h_intercept   -> SetMarkerSize(0.5);
}

void DummyRecHitAnalyzer::enable_pedestal_subtraction() { flag_perform_pedestal_subtraction = true; tag_calibration = "_ped_subtracted"; }

void DummyRecHitAnalyzer::enable_cm_subtraction() { flag_perform_cm_subtraction = true; tag_calibration = "_cm_subtracted"; }

Long64_t DummyRecHitAnalyzer::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
    }
    return centry;
}

void DummyRecHitAnalyzer::Load_metaData()
{
    TString csv_file_name = "./meta_conditions/calibration_parameters.csv";
    printf("[INFO] Load calibration parameters: %s\n", csv_file_name.Data());

    std::string line;
    std::ifstream loaded_csv_file(csv_file_name.Data());

    if(loaded_csv_file.is_open()) {
        while(getline(loaded_csv_file, line)) {
            // skip comments
            if(line.find("#")!=std::string::npos) continue;

            // it takes time to remove space... do not create space in the first place
            //std::string::iterator end_pos = std::remove(line.begin(), line.end(), ' ');
            //line.erase(end_pos, line.end());

            std::size_t found_1st_index = line.find(",");
            std::size_t found_2nd_index = line.find(",", found_1st_index+1, 1);
            std::size_t found_3rd_index = line.find(",", found_2nd_index+1, 1);
            std::size_t found_4th_index = line.find(",", found_3rd_index+1, 1);

            int channel_id   = std::stoi( line.substr(0,found_1st_index) );
            double pedestal  = std::stod( line.substr(found_1st_index+1, found_2nd_index) );
            double slope     = std::stod( line.substr(found_2nd_index+1, found_3rd_index) );
            double intercept = std::stod( line.substr(found_3rd_index+1, found_4th_index) );

            std::vector<double> v = {slope, intercept};
            map_pedestals[channel_id] = pedestal;
            map_cm_parameters[channel_id] = v;

            printf("channel_id = %d, pedestal = %.3f, slope = %.3f, intercept = %6.3f\n",
                    channel_id, map_pedestals[channel_id], map_cm_parameters[channel_id][0], map_cm_parameters[channel_id][1] );
        }
        loaded_csv_file.close();
    } else {
        std::cout << "[ERROR] unable to open " << csv_file_name.Data() << std::endl;
    }
}

void DummyRecHitAnalyzer::fill_histograms()
{
    myRecorder.add_entry(adc_channel_CM, adc_double);

    h_adc      -> Fill(adc_double);
    h_adcm     -> Fill(adcm);
    h_tot      -> Fill(tot);
    h_toa      -> Fill(toa);
    h_trigtime -> Fill(trigtime);
}

void DummyRecHitAnalyzer::fill_profiles(int globalChannelId_, double adc_double_)
{
    t_RecHit->Fill();

    myRunStatCollection.add_entry(globalChannelId_, adc_double, adc_channel_CM);

    p_adc      -> Fill(globalChannelId_ , adc_double , 1);
    p_adcm     -> Fill(globalChannelId_ , adcm       , 1);
    p_tot      -> Fill(globalChannelId_ , tot        , 1);
    p_toa      -> Fill(globalChannelId_ , toa        , 1);
    p_trigtime -> Fill(globalChannelId_ , trigtime   , 1);
}

void DummyRecHitAnalyzer::export_pedestals()
{
    if(myTag!="pedestal") return;

    TString csv_file_name = "./meta_conditions/output_EDAnalyzer_pedestals.csv";
    std::ofstream myfile(csv_file_name.Data());
    myfile << "#--------------------------------------------------\n";
    myfile << "# info: " << myTag.Data() << "\n";
    myfile << "# columns: channel, pedestal\n";
    myfile << "#--------------------------------------------------\n";

    // to-do: need to think how to extract pedestal from running parameters
    for(int i=0; i<234; ++i) {
        double mean = p_adc -> GetBinContent(i+1);
        myfile << Form("%d,%.2f\n", i, mean);
    }

    myfile.close();
    printf("[INFO] export pedestal: %s\n", csv_file_name.Data());
}

void DummyRecHitAnalyzer::export_cm_parameters()
{
    TString csv_file_name = "./meta_conditions/output_EDAnalyzer_cm_parameters" + tag_calibration + ".csv";
    std::ofstream myfile(csv_file_name.Data());
    myfile << "#--------------------------------------------------\n";
    myfile << "# info: " << myTag.Data() << "\n";
    myfile << "# columns: channel, slope, intercept, correlation\n";
    myfile << "#--------------------------------------------------\n";

    std::vector<RunningStatistics> mRs = myRunStatCollection.get_vector_running_statistics();

    for(int i=0; i<234; ++i) {
        myfile << Form("%d,%.2f,%.2f,%.2f\n", i, mRs[i].get_slope(), mRs[i].get_intercept(), mRs[i].get_correlation());
    }

    myfile.close();
    printf("[INFO] export CM parameters: %s\n", csv_file_name.Data());
}

void DummyRecHitAnalyzer::Show(Long64_t entry)
{
    std::cout << "event"        << " = " << event        << ", ";
    std::cout << "chip"         << " = " << chip         << ", ";
    std::cout << "half"         << " = " << half         << ", ";
    std::cout << "channel"      << " = " << channel      << ", ";
    std::cout << "adc"          << " = " << adc          << ", ";
    std::cout << "trigtime"     << " = " << trigtime     << std::endl;
    std::cout << "eventcounter" << " = " << eventcounter << std::endl;
    return;

    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}

#endif
