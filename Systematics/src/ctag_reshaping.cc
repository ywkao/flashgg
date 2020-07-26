#include "flashgg/Systematics/interface/ctag_reshaping.h"

retrieve_scale_factor::retrieve_scale_factor()
{
    printf("[check] This is the start of constructor\n");
    file = TFile::Open("flashgg/Systematics/rootfiles/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl.root");
    //file = TFile::Open("../rootfiles/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl.root");
    printf("[check] after open the root file...\n");

    //debug_ = false;
    //nbinx = 50;
    //nbiny = 50;
    //xmin = 0.;
    //xmax = 1.;
    //ymin = 0.;
    //ymax = 1.;
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

void retrieve_scale_factor::setcvsl_temporarycvsb(double cvsl_, double cvsb_)
{
    cvsl = cvsl_;
    cvsb = cvsb_;
    if(debug_) printf("cvsl = %.3f, cvsb = %.3f\n", cvsl, cvsb);
}

double retrieve_scale_factor::get_scale_factor(TString input, double cvsl, double cvsb)
{
    printf("[check] retrieve_scale_factor::get_scale_factor::OK\n");
    set_type_sys_uncertainty(input);
    printf("[check] retrieve_scale_factor::get_scale_factor::set_type_sys_uncertaint is OK\n");
    setcvsl_temporarycvsb(cvsl, cvsb);
    printf("[check] retrieve_scale_factor::get_scale_factor::setcvsl_temporarycvsb is OK\n");
    find_bin();
    printf("[check] retrieve_scale_factor::get_scale_factor::find_bin is OK\n");
    
    scale_factor = h->GetBinContent(bin_cvsl, bin_cvsb);
    if(debug_) printf("sf = %.2f\n", scale_factor);
    return scale_factor;
}

void retrieve_scale_factor::find_bin()
{
    bin_cvsl = h->GetXaxis()->FindBin(cvsl);
    bin_cvsb = h->GetYaxis()->FindBin(cvsb);
    if(debug_) printf("cvsl = %.3f (nbin = %d)\n", cvsl, bin_cvsl);
    if(debug_) printf("cvsb = %.3f (nbin = %d)\n", cvsb, bin_cvsb);

    //int nbinx = h->GetNbinsX();
    //int nbiny = h->GetNbinsY();
    //double xmin = h->GetXaxis()->GetXmin();
    //double xmax = h->GetXaxis()->GetXmax();
    //double ymin = h->GetYaxis()->GetXmin();
    //double ymax = h->GetYaxis()->GetXmax();

    //cvsl_temporary = cvsl;
    //cvsb_temporary = cvsb;

    //step = (ymax - ymin) / (double) nbiny;
    //counter = 0;
    //while(cvsb_temporary >= step)
    //{
    //    cvsb_temporary -= step;
    //    counter += 1;
    //}
    //bin_cvsb = counter + 1;

    //step = (xmax - xmin) / (double) nbinx;
    //counter = 0;
    //while(cvsl_temporary >= step)
    //{
    //    cvsl_temporary -= step;
    //    counter += 1;
    //}
    //bin_cvsl = counter + 1;
}

void retrieve_scale_factor::debug_mode()
{
    debug_ = true;
}
