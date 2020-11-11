#ifndef _CTAG_RESHAPING_H_
#define _CTAG_RESHAPING_H_

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <TFile.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TString.h>


class retrieve_scale_factor{
    public:
        retrieve_scale_factor();
        ~retrieve_scale_factor();

        void set_ctag_reshape_file(std::string path);
        void register_hists();
        int get_token(TString);
        double get_scale_factor(TString, double, double);
        void debug_mode();

    private:

        TFile *file;
        TH2D *h_register[87]; // (1 central + 14 sys. sources * 2 up/donw) *3 flavors = 87
        int bin_cvsl;
        int bin_cvsb;
        double scale_factor;

        bool debug_;
};

#endif
