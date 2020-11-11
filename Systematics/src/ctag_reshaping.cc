#include "flashgg/Systematics/interface/ctag_reshaping.h"

retrieve_scale_factor::retrieve_scale_factor()
{
    debug_ = false;
    if(debug_) std::cout << "This is from the retrieve_scale_factor constructor" << std::endl;
}

retrieve_scale_factor::~retrieve_scale_factor()
{
    file->Close();
    if(debug_) std::cout << "This is the end of retrieve_scale_factor destructor" << std::endl;
}

void retrieve_scale_factor::set_ctag_reshape_file(std::string path)
{
    file = TFile::Open(path.c_str());
    printf("[check] sf file: %s\n", path.c_str());

    register_hists();
}

double retrieve_scale_factor::get_scale_factor(TString input, double cvsl, double cvsb)
{
    int token = get_token(input);

    if(token > -1)
    {
        bin_cvsl = h_register[token]->GetXaxis()->FindBin(cvsl);
        bin_cvsb = h_register[token]->GetYaxis()->FindBin(cvsb);
        scale_factor = h_register[token]->GetBinContent(bin_cvsl, bin_cvsb);
        if(debug_) printf("[check] cvsl = %.3f, cvsb = %.3f, sf = %.2f (%s)\n", cvsl, cvsb, scale_factor, input.Data());
        return scale_factor;
    } else {
        printf("[WARNING] the input string does not match the registed SFs! return 1. (from retrieve_scale_factor::get_scale_factor) ");
        return 1.;
    }
}

void retrieve_scale_factor::register_hists()
{
    /* h_register will be registered with histograms in a specified SFs root files */
    TString type_flavour[3] = {"b", "c", "l"};

    TString type_uncertainty[14] = {"Stat", "EleIDSF", "LHEScaleWeight_muF", "LHEScaleWeight_muR", "MuIDSF", "PSWeightFSR", "PSWeightISR",
                                    "PUWeight", "XSec_DYJets", "XSec_ST", "XSec_WJets", "XSec_ttbar", "jer", "jesTotal"};

    TString type_variation[2] = {"Up", "Down"};
    
    int token = -1;
    for(int i=0; i<3; ++i)
    {
        // central value of SFs
        TString string_scale_factor = "SF" + type_flavour[i] + "_" + "hist";
        token = get_token(string_scale_factor);
        if(token > -1) h_register[token] = (TH2D*) file->Get(string_scale_factor);
        if(debug_) printf("token:%2d %s\n", token, string_scale_factor.Data());
        // different systematic sources
        for(int j=0; j<14; ++j)
        {
            // up / donw variation
            for(int k=0; k<2; ++k)
            {
                TString string_scale_factor_sys = "SF" + type_flavour[i] + "_" + "hist" + "_" + type_uncertainty[j] + type_variation[k];
                token = get_token(string_scale_factor_sys);
                if(token > -1) h_register[token] = (TH2D*) file->Get(string_scale_factor_sys);
                if(debug_) printf("token:%2d %s\n", token, string_scale_factor_sys.Data());
            }
        }
    }

}

void retrieve_scale_factor::debug_mode()
{
    debug_ = true;
}

int retrieve_scale_factor::get_token(TString input_name)
{
    if ( input_name == "SFb_hist" ) return 0;
    else if ( input_name == "SFb_hist_StatUp" ) return 1;
    else if ( input_name == "SFb_hist_EleIDSFUp" ) return 2;
    else if ( input_name == "SFb_hist_LHEScaleWeight_muFUp" ) return 3;
    else if ( input_name == "SFb_hist_LHEScaleWeight_muRUp" ) return 4;
    else if ( input_name == "SFb_hist_MuIDSFUp" ) return 5;
    else if ( input_name == "SFb_hist_PSWeightFSRUp" ) return 6;
    else if ( input_name == "SFb_hist_PSWeightISRUp" ) return 7;
    else if ( input_name == "SFb_hist_PUWeightUp" ) return 8;
    else if ( input_name == "SFb_hist_XSec_DYJetsUp" ) return 9;
    else if ( input_name == "SFb_hist_XSec_STUp" ) return 10;
    else if ( input_name == "SFb_hist_XSec_WJetsUp" ) return 11;
    else if ( input_name == "SFb_hist_XSec_ttbarUp" ) return 12;
    else if ( input_name == "SFb_hist_jerUp" ) return 13;
    else if ( input_name == "SFb_hist_jesTotalUp" ) return 14;
    else if ( input_name == "SFb_hist_StatDown" ) return 15;
    else if ( input_name == "SFb_hist_EleIDSFDown" ) return 16;
    else if ( input_name == "SFb_hist_LHEScaleWeight_muFDown" ) return 17;
    else if ( input_name == "SFb_hist_LHEScaleWeight_muRDown" ) return 18;
    else if ( input_name == "SFb_hist_MuIDSFDown" ) return 19;
    else if ( input_name == "SFb_hist_PSWeightFSRDown" ) return 20;
    else if ( input_name == "SFb_hist_PSWeightISRDown" ) return 21;
    else if ( input_name == "SFb_hist_PUWeightDown" ) return 22;
    else if ( input_name == "SFb_hist_XSec_DYJetsDown" ) return 23;
    else if ( input_name == "SFb_hist_XSec_STDown" ) return 24;
    else if ( input_name == "SFb_hist_XSec_WJetsDown" ) return 25;
    else if ( input_name == "SFb_hist_XSec_ttbarDown" ) return 26;
    else if ( input_name == "SFb_hist_jerDown" ) return 27;
    else if ( input_name == "SFb_hist_jesTotalDown" ) return 28;
    
    else if ( input_name == "SFc_hist" ) return 29;
    else if ( input_name == "SFc_hist_StatUp" ) return 30;
    else if ( input_name == "SFc_hist_EleIDSFUp" ) return 31;
    else if ( input_name == "SFc_hist_LHEScaleWeight_muFUp" ) return 32;
    else if ( input_name == "SFc_hist_LHEScaleWeight_muRUp" ) return 33;
    else if ( input_name == "SFc_hist_MuIDSFUp" ) return 34;
    else if ( input_name == "SFc_hist_PSWeightFSRUp" ) return 35;
    else if ( input_name == "SFc_hist_PSWeightISRUp" ) return 36;
    else if ( input_name == "SFc_hist_PUWeightUp" ) return 37;
    else if ( input_name == "SFc_hist_XSec_DYJetsUp" ) return 38;
    else if ( input_name == "SFc_hist_XSec_STUp" ) return 39;
    else if ( input_name == "SFc_hist_XSec_WJetsUp" ) return 40;
    else if ( input_name == "SFc_hist_XSec_ttbarUp" ) return 41;
    else if ( input_name == "SFc_hist_jerUp" ) return 42;
    else if ( input_name == "SFc_hist_jesTotalUp" ) return 43;
    else if ( input_name == "SFc_hist_StatDown" ) return 44;
    else if ( input_name == "SFc_hist_EleIDSFDown" ) return 45;
    else if ( input_name == "SFc_hist_LHEScaleWeight_muFDown" ) return 46;
    else if ( input_name == "SFc_hist_LHEScaleWeight_muRDown" ) return 47;
    else if ( input_name == "SFc_hist_MuIDSFDown" ) return 48;
    else if ( input_name == "SFc_hist_PSWeightFSRDown" ) return 49;
    else if ( input_name == "SFc_hist_PSWeightISRDown" ) return 50;
    else if ( input_name == "SFc_hist_PUWeightDown" ) return 51;
    else if ( input_name == "SFc_hist_XSec_DYJetsDown" ) return 52;
    else if ( input_name == "SFc_hist_XSec_STDown" ) return 53;
    else if ( input_name == "SFc_hist_XSec_WJetsDown" ) return 54;
    else if ( input_name == "SFc_hist_XSec_ttbarDown" ) return 55;
    else if ( input_name == "SFc_hist_jerDown" ) return 56;
    else if ( input_name == "SFc_hist_jesTotalDown" ) return 57;
    
    else if ( input_name == "SFl_hist" ) return 58;
    else if ( input_name == "SFl_hist_StatUp" ) return 59;
    else if ( input_name == "SFl_hist_EleIDSFUp" ) return 60;
    else if ( input_name == "SFl_hist_LHEScaleWeight_muFUp" ) return 61;
    else if ( input_name == "SFl_hist_LHEScaleWeight_muRUp" ) return 62;
    else if ( input_name == "SFl_hist_MuIDSFUp" ) return 63;
    else if ( input_name == "SFl_hist_PSWeightFSRUp" ) return 64;
    else if ( input_name == "SFl_hist_PSWeightISRUp" ) return 65;
    else if ( input_name == "SFl_hist_PUWeightUp" ) return 66;
    else if ( input_name == "SFl_hist_XSec_DYJetsUp" ) return 67;
    else if ( input_name == "SFl_hist_XSec_STUp" ) return 68;
    else if ( input_name == "SFl_hist_XSec_WJetsUp" ) return 69;
    else if ( input_name == "SFl_hist_XSec_ttbarUp" ) return 70;
    else if ( input_name == "SFl_hist_jerUp" ) return 71;
    else if ( input_name == "SFl_hist_jesTotalUp" ) return 72;
    else if ( input_name == "SFl_hist_StatDown" ) return 73;
    else if ( input_name == "SFl_hist_EleIDSFDown" ) return 74;
    else if ( input_name == "SFl_hist_LHEScaleWeight_muFDown" ) return 75;
    else if ( input_name == "SFl_hist_LHEScaleWeight_muRDown" ) return 76;
    else if ( input_name == "SFl_hist_MuIDSFDown" ) return 77;
    else if ( input_name == "SFl_hist_PSWeightFSRDown" ) return 78;
    else if ( input_name == "SFl_hist_PSWeightISRDown" ) return 79;
    else if ( input_name == "SFl_hist_PUWeightDown" ) return 80;
    else if ( input_name == "SFl_hist_XSec_DYJetsDown" ) return 81;
    else if ( input_name == "SFl_hist_XSec_STDown" ) return 82;
    else if ( input_name == "SFl_hist_XSec_WJetsDown" ) return 83;
    else if ( input_name == "SFl_hist_XSec_ttbarDown" ) return 84;
    else if ( input_name == "SFl_hist_jerDown" ) return 85;
    else if ( input_name == "SFl_hist_jesTotalDown" ) return 86;

    else return -1;
}
