#include "flashgg/Systematics/interface/ctag_reshaping.h"

retrieve_scale_factor::retrieve_scale_factor()
{
    file = TFile::Open("../data/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl.root");
    debug_ = false;
}

retrieve_scale_factor::~retrieve_scale_factor()
{
    file->Close();
    if(debug_) std::cout << "This is the end of destructor" << std::endl;
}

void retrieve_scale_factor::set_type_sys_uncertainty(TString input)
{
    type_sys_uncertainty = input;
    h = (TH2D*) file->Get(type_sys_uncertainty);
    if(debug_) printf("type:%s, ", type_sys_uncertainty.Data());
}

void retrieve_scale_factor::set_cvsl_cvsb(double cvsl_, double cvsb_)
{
    cvsl = cvsl_;
    cvsb = cvsb_;
    if(debug_) printf("cvsl = %.3f, cvsb = %.3f, ", cvsl, cvsb);
}

double retrieve_scale_factor::get_scale_factor(TString input, double cvsl, double cvsb)
{
    set_type_sys_uncertainty(input);
    set_cvsl_cvsb(cvsl, cvsb);
    find_bin();
    
    scale_factor = h->GetBinContent(bin_cvsl, bin_cvsb);
    if(debug_) printf("sf = %.2f\n", scale_factor);
    return scale_factor;
}

void retrieve_scale_factor::find_bin()
{
    bin_cvsl = h->GetXaxis()->FindBin(cvsl);
    bin_cvsb = h->GetYaxis()->FindBin(cvsb);
    if(debug_) printf("cvsl = %d\n", bin_cvsl);
    if(debug_) printf("cvsb = %d\n", bin_cvsb);
}

void retrieve_scale_factor::debug_mode()
{
    debug_ = true;
}
