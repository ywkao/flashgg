#ifndef _THQ_BDT_HELPER_
#define _THQ_BDT_HELPER_

#include <fstream>
#include <string>

#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TLorentzVector.h>

namespace flashgg {

// structure InputVariables{{{
typedef struct{
    float maxIDMVA_;
    float minIDMVA_;
    float max1_btag_;
    float max2_btag_;
    float dipho_delta_R;
    float njets_;
    float nbjets_;
    float ht_;
    float leadptoM_;
    float subleadptoM_;
    float lead_eta_;
    float sublead_eta_;
    float jet1_ptOverM_;
    float jet1_eta_;
    float jet1_btag_;
    float jet2_ptOverM_;
    float jet2_eta_;
    float jet2_btag_;
    float jet3_ptOverM_;
    float jet3_eta_;
    float jet3_btag_;
    float jet4_ptOverM_;
    float jet4_eta_;
    float jet4_btag_;
    float leadPSV_;
    float subleadPSV_;
    float dipho_cosphi_;
    float dipho_rapidity_;
    float met_;
    float dipho_pt_over_mass_;
    float helicity_angle_;
    float chi2_value_;
    float chi2_bjet_ptOverM_;
    float chi2_bjet_eta_;
    float chi2_bjet_btagScores_;
    float chi2_wjet1_ptOverM_;
    float chi2_wjet1_eta_;
    float chi2_wjet1_btagScores_;
    float chi2_wjet2_ptOverM_;
    float chi2_wjet2_eta_;
    float chi2_wjet2_btagScores_;
    float chi2_wjets_deltaR_;
    float chi2_wboson_ptOverM_;
    float chi2_wboson_eta_;
    float chi2_wboson_mass_;
    float chi2_wboson_deltaR_bjet_;
    float chi2_tbw_mass_;
    float chi2_tbw_ptOverM_;
    float chi2_tbw_eta_;
    float chi2_tprime_ptOverM_;
    float chi2_tprime_eta_;
    float tprime_pt_ratio_;
    float helicity_tprime_;
} InputVariables;
//}}}

class THQ_BDT_Helper {
public:
  THQ_BDT_Helper(const string& mva_algo, const std::string& weight_file_name);
  virtual ~THQ_BDT_Helper(){};
  InputVariables varList_;
  double evaluate(const std::string& tag, const InputVariables& varList);
  double convert_tmva_to_prob(double score);

  double chi2_calculator_2x2(double w_mass, double t_mass);
  bool get_the_best_bjj_candidate(std::vector<int> &indices_bjj, std::vector<TLorentzVector> jets, TLorentzVector diphoton, std::vector<float> btag_scores, double &min_chi2_value);

  double Rapidity(double energy, double Pz);
  double Helicity(const TLorentzVector particle_1, const TLorentzVector particle_2);

  void print_details(InputVariables var);

private:
  std::shared_ptr<TMVA::Reader> reader_ = nullptr;
};

// constructor{{{
THQ_BDT_Helper::THQ_BDT_Helper(const string& mva_algo, const std::string& weight_file_name)
{

  reader_ = std::make_shared<TMVA::Reader>( "!Color:!Silent" );

  reader_->AddVariable("maxIDMVA_"                 , &varList_.maxIDMVA_                 );
  reader_->AddVariable("minIDMVA_"                 , &varList_.minIDMVA_                 );
  reader_->AddVariable("max1_btag_"                , &varList_.max1_btag_                );
  reader_->AddVariable("max2_btag_"                , &varList_.max2_btag_                );
  reader_->AddVariable("dipho_delta_R"             , &varList_.dipho_delta_R             );
  reader_->AddVariable("njets_"                    , &varList_.njets_                    );
  reader_->AddVariable("nbjets_"                   , &varList_.nbjets_                   );
  reader_->AddVariable("ht_"                       , &varList_.ht_                       );
  reader_->AddVariable("leadptoM_"                 , &varList_.leadptoM_                 );
  reader_->AddVariable("subleadptoM_"              , &varList_.subleadptoM_              );
  reader_->AddVariable("lead_eta_"                 , &varList_.lead_eta_                 );
  reader_->AddVariable("sublead_eta_"              , &varList_.sublead_eta_              );
  reader_->AddVariable("jet1_ptOverM_"             , &varList_.jet1_ptOverM_             );
  reader_->AddVariable("jet1_eta_"                 , &varList_.jet1_eta_                 );
  reader_->AddVariable("jet1_btag_"                , &varList_.jet1_btag_                );
  reader_->AddVariable("jet2_ptOverM_"             , &varList_.jet2_ptOverM_             );
  reader_->AddVariable("jet2_eta_"                 , &varList_.jet2_eta_                 );
  reader_->AddVariable("jet2_btag_"                , &varList_.jet2_btag_                );
  reader_->AddVariable("jet3_ptOverM_"             , &varList_.jet3_ptOverM_             );
  reader_->AddVariable("jet3_eta_"                 , &varList_.jet3_eta_                 );
  reader_->AddVariable("jet3_btag_"                , &varList_.jet3_btag_                );
  reader_->AddVariable("jet4_ptOverM_"             , &varList_.jet4_ptOverM_             );
  reader_->AddVariable("jet4_eta_"                 , &varList_.jet4_eta_                 );
  reader_->AddVariable("jet4_btag_"                , &varList_.jet4_btag_                );
  reader_->AddVariable("leadPSV_"                  , &varList_.leadPSV_                  );
  reader_->AddVariable("subleadPSV_"               , &varList_.subleadPSV_               );
  reader_->AddVariable("dipho_cosphi_"             , &varList_.dipho_cosphi_             );
  reader_->AddVariable("dipho_rapidity_"           , &varList_.dipho_rapidity_           );
  reader_->AddVariable("met_"                      , &varList_.met_                      );
  reader_->AddVariable("dipho_pt_over_mass_"       , &varList_.dipho_pt_over_mass_       );
  reader_->AddVariable("helicity_angle_"           , &varList_.helicity_angle_           );
  reader_->AddVariable("chi2_value_"               , &varList_.chi2_value_               );
  reader_->AddVariable("chi2_bjet_ptOverM_"        , &varList_.chi2_bjet_ptOverM_        );
  reader_->AddVariable("chi2_bjet_eta_"            , &varList_.chi2_bjet_eta_            );
  reader_->AddVariable("chi2_bjet_btagScores_"     , &varList_.chi2_bjet_btagScores_     );
  reader_->AddVariable("chi2_wjet1_ptOverM_"       , &varList_.chi2_wjet1_ptOverM_       );
  reader_->AddVariable("chi2_wjet1_eta_"           , &varList_.chi2_wjet1_eta_           );
  reader_->AddVariable("chi2_wjet1_btagScores_"    , &varList_.chi2_wjet1_btagScores_    );
  reader_->AddVariable("chi2_wjet2_ptOverM_"       , &varList_.chi2_wjet2_ptOverM_       );
  reader_->AddVariable("chi2_wjet2_eta_"           , &varList_.chi2_wjet2_eta_           );
  reader_->AddVariable("chi2_wjet2_btagScores_"    , &varList_.chi2_wjet2_btagScores_    );
  reader_->AddVariable("chi2_wjets_deltaR_"        , &varList_.chi2_wjets_deltaR_        );
  reader_->AddVariable("chi2_wboson_ptOverM_"      , &varList_.chi2_wboson_ptOverM_      );
  reader_->AddVariable("chi2_wboson_eta_"          , &varList_.chi2_wboson_eta_          );
  reader_->AddVariable("chi2_wboson_mass_"         , &varList_.chi2_wboson_mass_         );
  reader_->AddVariable("chi2_wboson_deltaR_bjet_"  , &varList_.chi2_wboson_deltaR_bjet_  );
  reader_->AddVariable("chi2_tbw_mass_"            , &varList_.chi2_tbw_mass_            );
  reader_->AddVariable("chi2_tbw_ptOverM_"         , &varList_.chi2_tbw_ptOverM_         );
  reader_->AddVariable("chi2_tbw_eta_"             , &varList_.chi2_tbw_eta_             );
  reader_->AddVariable("chi2_tprime_ptOverM_"      , &varList_.chi2_tprime_ptOverM_      );
  reader_->AddVariable("chi2_tprime_eta_"          , &varList_.chi2_tprime_eta_          );
  reader_->AddVariable("tprime_pt_ratio_"          , &varList_.tprime_pt_ratio_          );
  reader_->AddVariable("helicity_tprime_"          , &varList_.helicity_tprime_          );

  reader_->BookMVA(mva_algo.c_str(), weight_file_name);
}
//}}}

double THQ_BDT_Helper::evaluate(const string& mva_algo, const InputVariables& varList) {
  memcpy(&varList_, &varList, sizeof(varList)); // use move syntax here                                                                                                      
  return reader_->EvaluateMVA(mva_algo.c_str());
}

double THQ_BDT_Helper::convert_tmva_to_prob(double score)
{
  // Undo TMVA transformation
  double raw_score = -0.5 * log( (2 / (score + 1)) - 1);
  // Apply logistic (sigmoid) transformation
  double prob = 1 / (1 + exp(-raw_score));
  return prob;
}

//----------------------------------------------------------------------------------------------------
// auxiliary functions
//----------------------------------------------------------------------------------------------------
inline
double THQ_BDT_Helper::Rapidity(double energy, double Pz)
{
    double rapidity = 0.5 * log( (energy+Pz) / (energy-Pz) );
    return rapidity;
};

inline
double THQ_BDT_Helper::Helicity(const TLorentzVector particle_1, const TLorentzVector particle_2)
{
  TLorentzVector p1 = particle_1;
  TLorentzVector parent = particle_1 + particle_2;
  
  TVector3 boost_to_parent = -(parent.BoostVector());
  p1.Boost(boost_to_parent);

  TVector3 v1 = p1.Vect();
  TVector3 vParent = parent.Vect();

  if( isnan(v1.Mag()) ) return -999;

  double cos_theta_1 = (v1.Dot(vParent)) / (v1.Mag() * vParent.Mag());
  return abs(cos_theta_1);  
};

inline
double THQ_BDT_Helper::chi2_calculator_2x2(double w_mass, double t_mass)
{
    // means & covariance matrix elements are taken from covMatrx_Era2017_M1000.json
    TVectorD vec_mean_values(2);

    vec_mean_values(0) = 85.70;
    vec_mean_values(1) = 174.81;

    TVectorD vec_mass(2);
    vec_mass(0) = w_mass - vec_mean_values(0);
    vec_mass(1) = t_mass - vec_mean_values(1);

    TMatrixD matrix(2,2);
    matrix(0,0) = 305.14;
    matrix(0,1) = 282.18;
    matrix(1,0) = 282.18;
    matrix(1,1) = 572.63;

    double chi2_value = matrix.Invert()*vec_mass*vec_mass;
    return chi2_value;
};

inline
bool THQ_BDT_Helper::get_the_best_bjj_candidate(std::vector<int> &indices_bjj, std::vector<TLorentzVector> jets, TLorentzVector diphoton, std::vector<float> btag_scores, double &min_chi2_value)
{
    double pfDeepCSV_btag_loose_wp = 0.1522;

    std::size_t num_jets = jets.size();
    for(std::size_t i = 0; i < num_jets; ++i ){ // b-jet
        if (btag_scores[i] < pfDeepCSV_btag_loose_wp) continue;
        for(std::size_t j = 0; j < num_jets-1; ++j ){ // w-jet1
            if(j==i) continue;
            for(std::size_t k = j+1; k < num_jets; ++k ){ // w-jet2
                if(k==i) continue;

                TLorentzVector bjet = jets[i];
                TLorentzVector w_candidate = jets[j] + jets[k];
                TLorentzVector top_candidate = w_candidate + bjet;
                double w_mass = w_candidate.M();
                double t_mass = top_candidate.M();
                double chi2 = chi2_calculator_2x2(w_mass, t_mass);
                if(chi2 < min_chi2_value){
                    indices_bjj[0] = i;
                    indices_bjj[1] = j;
                    indices_bjj[2] = k;
                    min_chi2_value = chi2;
                }
            }
        }
    }//end of looping jets

    bool has_resonable_reco = indices_bjj[0] >= 0 && indices_bjj[1] >= 0 && indices_bjj[2] >= 0;
    return has_resonable_reco;
};

// print details{{{
inline
void THQ_BDT_Helper::print_details(InputVariables var)
{
    printf("%s: %.10f , " , "maxIDMVA_"                , var.maxIDMVA_                );
    printf("%s: %.10f , " , "minIDMVA_"                , var.minIDMVA_                );
    printf("%s: %.10f , " , "max1_btag_"               , var.max1_btag_               );
    printf("%s: %.10f , " , "max2_btag_"               , var.max2_btag_               );
    printf("%s: %.10f , " , "dipho_delta_R"            , var.dipho_delta_R            );
    printf("%s: %.10f , " , "njets_"                   , var.njets_                   );
    printf("%s: %.10f , " , "nbjets_"                  , var.nbjets_                  );
    printf("%s: %.10f , " , "ht_"                      , var.ht_                      );
    printf("%s: %.10f , " , "leadptoM_"                , var.leadptoM_                );
    printf("%s: %.10f , " , "subleadptoM_"             , var.subleadptoM_             );
    printf("%s: %.10f , " , "lead_eta_"                , var.lead_eta_                );
    printf("%s: %.10f , " , "sublead_eta_"             , var.sublead_eta_             );
    printf("%s: %.10f , " , "jet1_ptOverM_"            , var.jet1_ptOverM_            );
    printf("%s: %.10f , " , "jet1_eta_"                , var.jet1_eta_                );
    printf("%s: %.10f , " , "jet1_btag_"               , var.jet1_btag_               );
    printf("%s: %.10f , " , "jet2_ptOverM_"            , var.jet2_ptOverM_            );
    printf("%s: %.10f , " , "jet2_eta_"                , var.jet2_eta_                );
    printf("%s: %.10f , " , "jet2_btag_"               , var.jet2_btag_               );
    printf("%s: %.10f , " , "jet3_ptOverM_"            , var.jet3_ptOverM_            );
    printf("%s: %.10f , " , "jet3_eta_"                , var.jet3_eta_                );
    printf("%s: %.10f , " , "jet3_btag_"               , var.jet3_btag_               );
    printf("%s: %.10f , " , "jet4_ptOverM_"            , var.jet4_ptOverM_            );
    printf("%s: %.10f , " , "jet4_eta_"                , var.jet4_eta_                );
    printf("%s: %.10f , " , "jet4_btag_"               , var.jet4_btag_               );
    printf("%s: %.10f , " , "leadPSV_"                 , var.leadPSV_                 );
    printf("%s: %.10f , " , "subleadPSV_"              , var.subleadPSV_              );
    printf("%s: %.10f , " , "dipho_cosphi_"            , var.dipho_cosphi_            );
    printf("%s: %.10f , " , "dipho_rapidity_"          , var.dipho_rapidity_          );
    printf("%s: %.10f , " , "met_"                     , var.met_                     );
    printf("%s: %.10f , " , "dipho_pt_over_mass_"      , var.dipho_pt_over_mass_      );
    printf("%s: %.10f , " , "helicity_angle_"          , var.helicity_angle_          );
    printf("%s: %.10f , " , "chi2_value_"              , var.chi2_value_              );
    printf("%s: %.10f , " , "chi2_bjet_ptOverM_"       , var.chi2_bjet_ptOverM_       );
    printf("%s: %.10f , " , "chi2_bjet_eta_"           , var.chi2_bjet_eta_           );
    printf("%s: %.10f , " , "chi2_bjet_btagScores_"    , var.chi2_bjet_btagScores_    );
    printf("%s: %.10f , " , "chi2_wjet1_ptOverM_"      , var.chi2_wjet1_ptOverM_      );
    printf("%s: %.10f , " , "chi2_wjet1_eta_"          , var.chi2_wjet1_eta_          );
    printf("%s: %.10f , " , "chi2_wjet1_btagScores_"   , var.chi2_wjet1_btagScores_   );
    printf("%s: %.10f , " , "chi2_wjet2_ptOverM_"      , var.chi2_wjet2_ptOverM_      );
    printf("%s: %.10f , " , "chi2_wjet2_eta_"          , var.chi2_wjet2_eta_          );
    printf("%s: %.10f , " , "chi2_wjet2_btagScores_"   , var.chi2_wjet2_btagScores_   );
    printf("%s: %.10f , " , "chi2_wjets_deltaR_"       , var.chi2_wjets_deltaR_       );
    printf("%s: %.10f , " , "chi2_wboson_ptOverM_"     , var.chi2_wboson_ptOverM_     );
    printf("%s: %.10f , " , "chi2_wboson_eta_"         , var.chi2_wboson_eta_         );
    printf("%s: %.10f , " , "chi2_wboson_mass_"        , var.chi2_wboson_mass_        );
    printf("%s: %.10f , " , "chi2_wboson_deltaR_bjet_" , var.chi2_wboson_deltaR_bjet_ );
    printf("%s: %.10f , " , "chi2_tbw_mass_"           , var.chi2_tbw_mass_           );
    printf("%s: %.10f , " , "chi2_tbw_ptOverM_"        , var.chi2_tbw_ptOverM_        );
    printf("%s: %.10f , " , "chi2_tbw_eta_"            , var.chi2_tbw_eta_            );
    printf("%s: %.10f , " , "chi2_tprime_ptOverM_"     , var.chi2_tprime_ptOverM_     );
    printf("%s: %.10f , " , "chi2_tprime_eta_"         , var.chi2_tprime_eta_         );
    printf("%s: %.10f , " , "tprime_pt_ratio_"         , var.tprime_pt_ratio_         );
    printf("%s: %.10f , " , "helicity_tprime_"         , var.helicity_tprime_         );
    printf("\n\n");
};
//}}}

} // end of namespace

#endif // _THQ_BDT_HELPER_
