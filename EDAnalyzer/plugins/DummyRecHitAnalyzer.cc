#include "calibration_algorithms/EDAnalyzer/interface/DummyRecHitAnalyzer.h"

// ------------ method called for each event  ------------
void DummyRecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;

    printf(">>> DummyRecHitAnalyzer::analyze: Hello World!\n");

    //for (const auto& track : iEvent.get(tracksToken_)) {
    //  // do something with track parameters, e.g, plot the charge.
    //  // int charge = track.charge();
    //}

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    // if the SetupData is always needed
    auto setup = iSetup.getData(setupToken_);
    // if need the ESHandle to check if the SetupData was there or not
    auto pSetup = iSetup.getHandle(setupToken_);
#endif

    //--------------------------------------------------
    // from standalone c++ code
    //--------------------------------------------------
    Long64_t nentries = fChain->GetEntriesFast();
    int nevent = nentries / 78;

    std::cout << ">>> nentries = " << nentries << std::endl;
    std::cout << ">>> nevent = " << nevent << std::endl;

    // loop over hit
    std::vector<Hit> hits[nevent];
    int counter           = 0;
    int current_half      = 0;
    int recorded_half     = 0;
    double adc_channel_37 = 0.;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry); nbytes += nb;

        // reset cm adc when reaching another ROC half
        current_half = half;
        if(current_half != recorded_half) {
            recorded_half = current_half;
            adc_channel_CM = 0.;
        }

        // get detId & mask bad channel
        auto detid = DetectorId( FromRawData(), chip, half, channel );
        globalChannelId = detid.id(); // chip*78+half*39+channel;
        bool is_bad_channel = globalChannelId==146 || globalChannelId==171;
        if(is_bad_channel) continue;

        // convert adc to double
        adc_double = (double) adc;

        // perform pedestal subtraction
        if(flag_perform_pedestal_subtraction) {
            double pedestal = map_pedestals[globalChannelId];
            adc_double -= pedestal;
        }

        // handle cm information after pedestal subtraction
        bool is_cm_channel = (globalChannelId % 39 == 37 || globalChannelId % 39 == 38);
        if(is_cm_channel) {
            // take average of two cm channels in a half
            adc_channel_CM += adc_double / 2.;
        }

        // record adc of ch37 & fill info of ch37 when processing ch38
        if(globalChannelId % 39 == 37) {
            adc_channel_37 = adc_double;
            continue;

        } else if(globalChannelId % 39 == 38) {
            // CM subtraction for channel 37
            if(flag_perform_cm_subtraction) {
                std::vector<double> parameters = map_cm_parameters[globalChannelId-1];
                double slope = parameters[0];
                double intercept = parameters[1];
                double correction = adc_channel_CM*slope + intercept;
                adc_channel_37 -= correction;
            }

            fill_profiles(globalChannelId-1, adc_channel_37);
        }

        // perform common mode subtraction
        if(flag_perform_cm_subtraction) {
            std::vector<double> parameters = map_cm_parameters[globalChannelId];
            double slope = parameters[0];
            double intercept = parameters[1];
            double correction = adc_channel_CM*slope + intercept;
            adc_double -= correction;
        }

        fill_profiles(globalChannelId, adc_double);

        if(globalChannelId==7) fill_histograms();

        // store information (Arnaud's code)
        Hit hit( event, detid, adc, toa, tot, trigtime ); // Note: need subtraction from pedestal run
        hits[event].push_back( hit );

        continue;

        // print event info
        Show();
        counter+=1;
        if(counter>=1000) break;
    }
}

// ------------ method called once each job just before starting event loop  ------------
void DummyRecHitAnalyzer::beginJob() {
    //--------------------------------------------------
    // load trees from beam data / pedestal run
    //--------------------------------------------------
    TString root_beamRun  = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2022/sps_oct2022/pion_beam_150_320fC/beam_run/run_20221007_191926/beam_run0.root";
    TString root_pedestal = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2022/sps_oct2022/pedestals/pedestal_320fC/pedestal_run/run_20221008_192720/pedestal_run0.root";
    TString input = (myTag=="beam") ? root_beamRun : root_pedestal;
    printf("[INFO] Input rootfile: %s\n", input.Data());

    f1 = new TFile(input, "R");
    t1 = (TTree*) f1->Get("unpacker_data/hgcroc");
    Init(t1); // SetBranchAddress, init variables, etc.
    Init_my_output_info(); // Register histograms, etc. -> to be replaced in DQMEDAnalyzer

    // determine which metaData to load, set boolean and tag of calibrations
    if(calibration_flags[0]) enable_pedestal_subtraction();
    if(calibration_flags[1]) enable_cm_subtraction();
    Load_metaData();
}

// ------------ method called once each job just after ending the event loop  ------------
void DummyRecHitAnalyzer::endJob() {
    // prepare summary for running statistics
    std::vector<RunningStatistics> mRs = myRunStatCollection.get_vector_running_statistics();
    for(int channelId=0; channelId<234; ++channelId) {
        h_correlation -> SetBinContent( channelId+1, mRs[channelId].get_correlation() );
        h_slope       -> SetBinContent( channelId+1, mRs[channelId].get_slope()       );
        h_intercept   -> SetBinContent( channelId+1, mRs[channelId].get_intercept()   );
    }

    if(myTag=="pedestal" && !flag_perform_pedestal_subtraction) { export_pedestals(); }
    if(myTag=="beam" && !flag_perform_cm_subtraction) { export_cm_parameters(); }

    printf("[INFO] This is the end of the job\n");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DummyRecHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);

    //Specify that only 'tracks' is allowed
    //To use, remove the default given above and uncomment below
    //ParameterSetDescription desc;
    //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
    //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DummyRecHitAnalyzer);
