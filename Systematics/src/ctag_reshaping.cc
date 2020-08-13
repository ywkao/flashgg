#include "flashgg/Systematics/interface/ctag_reshaping.h"

retrieve_scale_factor::retrieve_scale_factor(std::string path)
{
    debug_ = false;
    file = TFile::Open(path.c_str());
    printf("[check] retrieve_scale_factor::retrieve_scale_factor\n");
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
    if(debug_) printf("[check] type: %s\n", type_sys_uncertainty.Data());
}

void retrieve_scale_factor::set_cvsl_cvsb(double cvsl_, double cvsb_)
{
    cvsl = cvsl_;
    cvsb = cvsb_;
    if(debug_) printf("[check] cvsl = %.3f, cvsb = %.3f\n", cvsl, cvsb);
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
    if(debug_) printf("[check] bin_cvsl = %d, bin_cvsb = %d\n", bin_cvsl, bin_cvsb);
}

void retrieve_scale_factor::debug_mode()
{
    debug_ = true;
}
