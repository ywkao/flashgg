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
// structure MoreVariables{{{
typedef struct{
    float diphoton_ptOverM_;
    float diphoton_pt_;
    float diphoton_eta_;
    float diphoton_phi_;
    float diphoton_energy_;
    float diphoton_mass_;
    float wboson_ptOverM_;
    float wboson_pt_;
    float wboson_eta_;
    float wboson_phi_;
    float wboson_energy_;
    float wboson_mass_;
    float top_ptOverM_;
    float top_pt_;
    float top_eta_;
    float top_phi_;
    float top_energy_;
    float top_mass_;
    float tprime_ptOverM_;
    float tprime_pt_;
    float tprime_eta_;
    float tprime_phi_;
    float tprime_energy_;
    float tprime_mass_;
    float bjet_ptOverM_;
    float bjet_pt_;
    float bjet_eta_;
    float bjet_phi_;
    float bjet_energy_;
    float bjet_mass_;
    float bjet_btag_;
    float wjet1_ptOverM_;
    float wjet1_pt_;
    float wjet1_eta_;
    float wjet1_phi_;
    float wjet1_energy_;
    float wjet1_mass_;
    float wjet1_btag_;
    float wjet2_ptOverM_;
    float wjet2_pt_;
    float wjet2_eta_;
    float wjet2_phi_;
    float wjet2_energy_;
    float wjet2_mass_;
    float wjet2_btag_;
    
    float jet1_pt_;
    float jet1_eta_;
    float jet1_phi_;
    float jet1_energy_;
    float jet1_mass_;
    float jet1_btag_;
    float jet2_pt_;
    float jet2_eta_;
    float jet2_phi_;
    float jet2_energy_;
    float jet2_mass_;
    float jet2_btag_;
    float jet3_pt_;
    float jet3_eta_;
    float jet3_phi_;
    float jet3_energy_;
    float jet3_mass_;
    float jet3_btag_;
    float jet4_pt_;
    float jet4_eta_;
    float jet4_phi_;
    float jet4_energy_;
    float jet4_mass_;
    float jet4_btag_;
    float jet5_pt_;
    float jet5_eta_;
    float jet5_phi_;
    float jet5_energy_;
    float jet5_mass_;
    float jet5_btag_;
    float jet6_pt_;
    float jet6_eta_;
    float jet6_phi_;
    float jet6_energy_;
    float jet6_mass_;
    float jet6_btag_;
    float jet7_pt_;
    float jet7_eta_;
    float jet7_phi_;
    float jet7_energy_;
    float jet7_mass_;
    float jet7_btag_;
    float jet8_pt_;
    float jet8_eta_;
    float jet8_phi_;
    float jet8_energy_;
    float jet8_mass_;
    float jet8_btag_;
    float jet9_pt_;
    float jet9_eta_;
    float jet9_phi_;
    float jet9_energy_;
    float jet9_mass_;
    float jet9_btag_;
    float jet10_pt_;
    float jet10_eta_;
    float jet10_phi_;
    float jet10_energy_;
    float jet10_mass_;
    float jet10_btag_;
    float jet11_pt_;
    float jet11_eta_;
    float jet11_phi_;
    float jet11_energy_;
    float jet11_mass_;
    float jet11_btag_;
    float jet12_pt_;
    float jet12_eta_;
    float jet12_phi_;
    float jet12_energy_;
    float jet12_mass_;
    float jet12_btag_;
    float jet13_pt_;
    float jet13_eta_;
    float jet13_phi_;
    float jet13_energy_;
    float jet13_mass_;
    float jet13_btag_;
    float jet14_pt_;
    float jet14_eta_;
    float jet14_phi_;
    float jet14_energy_;
    float jet14_mass_;
    float jet14_btag_;
    float jet15_pt_;
    float jet15_eta_;
    float jet15_phi_;
    float jet15_energy_;
    float jet15_mass_;
    float jet15_btag_;

} MoreVariables;
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
  void print_details_cout(InputVariables var);
  void print_details_more(MoreVariables var);
  void print_details_jets(MoreVariables var);

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
    //printf("\n\n");
};
//}}}
// print details cout{{{
inline
void THQ_BDT_Helper::print_details_cout(InputVariables var)
{
    std::cout << "maxIDMVA_: "                << var.maxIDMVA_                << ", ";
    std::cout << "minIDMVA_: "                << var.minIDMVA_                << ", ";
    std::cout << "max1_btag_: "               << var.max1_btag_               << ", ";
    std::cout << "max2_btag_: "               << var.max2_btag_               << ", ";
    std::cout << "dipho_delta_R: "            << var.dipho_delta_R            << ", ";
    std::cout << "njets_: "                   << var.njets_                   << ", ";
    std::cout << "nbjets_: "                  << var.nbjets_                  << ", ";
    std::cout << "ht_: "                      << var.ht_                      << ", ";
    std::cout << "leadptoM_: "                << var.leadptoM_                << ", ";
    std::cout << "subleadptoM_: "             << var.subleadptoM_             << ", ";
    std::cout << "lead_eta_: "                << var.lead_eta_                << ", ";
    std::cout << "sublead_eta_: "             << var.sublead_eta_             << ", ";
    std::cout << "jet1_ptOverM_: "            << var.jet1_ptOverM_            << ", ";
    std::cout << "jet1_eta_: "                << var.jet1_eta_                << ", ";
    std::cout << "jet1_btag_: "               << var.jet1_btag_               << ", ";
    std::cout << "jet2_ptOverM_: "            << var.jet2_ptOverM_            << ", ";
    std::cout << "jet2_eta_: "                << var.jet2_eta_                << ", ";
    std::cout << "jet2_btag_: "               << var.jet2_btag_               << ", ";
    std::cout << "jet3_ptOverM_: "            << var.jet3_ptOverM_            << ", ";
    std::cout << "jet3_eta_: "                << var.jet3_eta_                << ", ";
    std::cout << "jet3_btag_: "               << var.jet3_btag_               << ", ";
    std::cout << "jet4_ptOverM_: "            << var.jet4_ptOverM_            << ", ";
    std::cout << "jet4_eta_: "                << var.jet4_eta_                << ", ";
    std::cout << "jet4_btag_: "               << var.jet4_btag_               << ", ";
    std::cout << "leadPSV_: "                 << var.leadPSV_                 << ", ";
    std::cout << "subleadPSV_: "              << var.subleadPSV_              << ", ";
    std::cout << "dipho_cosphi_: "            << var.dipho_cosphi_            << ", ";
    std::cout << "dipho_rapidity_: "          << var.dipho_rapidity_          << ", ";
    std::cout << "met_: "                     << var.met_                     << ", ";
    std::cout << "dipho_pt_over_mass_: "      << var.dipho_pt_over_mass_      << ", ";
    std::cout << "helicity_angle_: "          << var.helicity_angle_          << ", ";
    std::cout << "chi2_value_: "              << var.chi2_value_              << ", ";
    std::cout << "chi2_bjet_ptOverM_: "       << var.chi2_bjet_ptOverM_       << ", ";
    std::cout << "chi2_bjet_eta_: "           << var.chi2_bjet_eta_           << ", ";
    std::cout << "chi2_bjet_btagScores_: "    << var.chi2_bjet_btagScores_    << ", ";
    std::cout << "chi2_wjet1_ptOverM_: "      << var.chi2_wjet1_ptOverM_      << ", ";
    std::cout << "chi2_wjet1_eta_: "          << var.chi2_wjet1_eta_          << ", ";
    std::cout << "chi2_wjet1_btagScores_: "   << var.chi2_wjet1_btagScores_   << ", ";
    std::cout << "chi2_wjet2_ptOverM_: "      << var.chi2_wjet2_ptOverM_      << ", ";
    std::cout << "chi2_wjet2_eta_: "          << var.chi2_wjet2_eta_          << ", ";
    std::cout << "chi2_wjet2_btagScores_: "   << var.chi2_wjet2_btagScores_   << ", ";
    std::cout << "chi2_wjets_deltaR_: "       << var.chi2_wjets_deltaR_       << ", ";
    std::cout << "chi2_wboson_ptOverM_: "     << var.chi2_wboson_ptOverM_     << ", ";
    std::cout << "chi2_wboson_eta_: "         << var.chi2_wboson_eta_         << ", ";
    std::cout << "chi2_wboson_mass_: "        << var.chi2_wboson_mass_        << ", ";
    std::cout << "chi2_wboson_deltaR_bjet_: " << var.chi2_wboson_deltaR_bjet_ << ", ";
    std::cout << "chi2_tbw_mass_: "           << var.chi2_tbw_mass_           << ", ";
    std::cout << "chi2_tbw_ptOverM_: "        << var.chi2_tbw_ptOverM_        << ", ";
    std::cout << "chi2_tbw_eta_: "            << var.chi2_tbw_eta_            << ", ";
    std::cout << "chi2_tprime_ptOverM_: "     << var.chi2_tprime_ptOverM_     << ", ";
    std::cout << "chi2_tprime_eta_: "         << var.chi2_tprime_eta_         << ", ";
    std::cout << "tprime_pt_ratio_: "         << var.tprime_pt_ratio_         << ", ";
    std::cout << "helicity_tprime_: "         << var.helicity_tprime_         << ", ";
};
//}}}
// print details more{{{
inline
void THQ_BDT_Helper::print_details_more(MoreVariables var)
{
    std::cout << "diphoton_ptOverM_: " << var.diphoton_ptOverM_ << ", ";
    std::cout << "diphoton_pt_: "      << var.diphoton_pt_      << ", ";
    std::cout << "diphoton_eta_: "     << var.diphoton_eta_     << ", ";
    std::cout << "diphoton_phi_: "     << var.diphoton_phi_     << ", ";
    std::cout << "diphoton_energy_: "  << var.diphoton_energy_  << ", ";
    std::cout << "diphoton_mass_: "    << var.diphoton_mass_    << ", ";
    std::cout << "wboson_ptOverM_: "   << var.wboson_ptOverM_   << ", ";
    std::cout << "wboson_pt_: "        << var.wboson_pt_        << ", ";
    std::cout << "wboson_eta_: "       << var.wboson_eta_       << ", ";
    std::cout << "wboson_phi_: "       << var.wboson_phi_       << ", ";
    std::cout << "wboson_energy_: "    << var.wboson_energy_    << ", ";
    std::cout << "wboson_mass_: "      << var.wboson_mass_      << ", ";
    std::cout << "top_ptOverM_: "      << var.top_ptOverM_      << ", ";
    std::cout << "top_pt_: "           << var.top_pt_           << ", ";
    std::cout << "top_eta_: "          << var.top_eta_          << ", ";
    std::cout << "top_phi_: "          << var.top_phi_          << ", ";
    std::cout << "top_energy_: "       << var.top_energy_       << ", ";
    std::cout << "top_mass_: "         << var.top_mass_         << ", ";
    std::cout << "tprime_ptOverM_: "   << var.tprime_ptOverM_   << ", ";
    std::cout << "tprime_pt_: "        << var.tprime_pt_        << ", ";
    std::cout << "tprime_eta_: "       << var.tprime_eta_       << ", ";
    std::cout << "tprime_phi_: "       << var.tprime_phi_       << ", ";
    std::cout << "tprime_energy_: "    << var.tprime_energy_    << ", ";
    std::cout << "tprime_mass_: "      << var.tprime_mass_      << ", ";
    std::cout << "bjet_ptOverM_: "     << var.bjet_ptOverM_     << ", ";
    std::cout << "bjet_pt_: "          << var.bjet_pt_          << ", ";
    std::cout << "bjet_eta_: "         << var.bjet_eta_         << ", ";
    std::cout << "bjet_phi_: "         << var.bjet_phi_         << ", ";
    std::cout << "bjet_energy_: "      << var.bjet_energy_      << ", ";
    std::cout << "bjet_mass_: "        << var.bjet_mass_        << ", ";
    std::cout << "bjet_btag_: "        << var.bjet_btag_        << ", ";
    std::cout << "wjet1_ptOverM_: "    << var.wjet1_ptOverM_    << ", ";
    std::cout << "wjet1_pt_: "         << var.wjet1_pt_         << ", ";
    std::cout << "wjet1_eta_: "        << var.wjet1_eta_        << ", ";
    std::cout << "wjet1_phi_: "        << var.wjet1_phi_        << ", ";
    std::cout << "wjet1_energy_: "     << var.wjet1_energy_     << ", ";
    std::cout << "wjet1_mass_: "       << var.wjet1_mass_       << ", ";
    std::cout << "wjet1_btag_: "       << var.wjet1_btag_       << ", ";
    std::cout << "wjet2_ptOverM_: "    << var.wjet2_ptOverM_    << ", ";
    std::cout << "wjet2_pt_: "         << var.wjet2_pt_         << ", ";
    std::cout << "wjet2_eta_: "        << var.wjet2_eta_        << ", ";
    std::cout << "wjet2_phi_: "        << var.wjet2_phi_        << ", ";
    std::cout << "wjet2_energy_: "     << var.wjet2_energy_     << ", ";
    std::cout << "wjet2_mass_: "       << var.wjet2_mass_       << ", ";
    std::cout << "wjet2_btag_: "       << var.wjet2_btag_       << ", ";
    std::cout << "jet1_pt_: "          << var.jet1_pt_          << ", ";
    std::cout << "jet1_eta_: "         << var.jet1_eta_         << ", ";
    std::cout << "jet1_phi_: "         << var.jet1_phi_         << ", ";
    std::cout << "jet1_energy_: "      << var.jet1_energy_      << ", ";
    std::cout << "jet1_mass_: "        << var.jet1_mass_        << ", ";
    std::cout << "jet1_btag_: "        << var.jet1_btag_        << ", ";
    std::cout << "jet2_pt_: "          << var.jet2_pt_          << ", ";
    std::cout << "jet2_eta_: "         << var.jet2_eta_         << ", ";
    std::cout << "jet2_phi_: "         << var.jet2_phi_         << ", ";
    std::cout << "jet2_energy_: "      << var.jet2_energy_      << ", ";
    std::cout << "jet2_mass_: "        << var.jet2_mass_        << ", ";
    std::cout << "jet2_btag_: "        << var.jet2_btag_        << ", ";
    std::cout << "jet3_pt_: "          << var.jet3_pt_          << ", ";
    std::cout << "jet3_eta_: "         << var.jet3_eta_         << ", ";
    std::cout << "jet3_phi_: "         << var.jet3_phi_         << ", ";
    std::cout << "jet3_energy_: "      << var.jet3_energy_      << ", ";
    std::cout << "jet3_mass_: "        << var.jet3_mass_        << ", ";
    std::cout << "jet3_btag_: "        << var.jet3_btag_        << ", ";
    std::cout << "jet4_pt_: "          << var.jet4_pt_          << ", ";
    std::cout << "jet4_eta_: "         << var.jet4_eta_         << ", ";
    std::cout << "jet4_phi_: "         << var.jet4_phi_         << ", ";
    std::cout << "jet4_energy_: "      << var.jet4_energy_      << ", ";
    std::cout << "jet4_mass_: "        << var.jet4_mass_        << ", ";
    std::cout << "jet4_btag_: "        << var.jet4_btag_        << ", ";
    std::cout << "jet5_pt_: "          << var.jet5_pt_          << ", ";
    std::cout << "jet5_eta_: "         << var.jet5_eta_         << ", ";
    std::cout << "jet5_phi_: "         << var.jet5_phi_         << ", ";
    std::cout << "jet5_energy_: "      << var.jet5_energy_      << ", ";
    std::cout << "jet5_mass_: "        << var.jet5_mass_        << ", ";
    std::cout << "jet5_btag_: "        << var.jet5_btag_        << ", ";
    std::cout << "jet6_pt_: "          << var.jet6_pt_          << ", ";
    std::cout << "jet6_eta_: "         << var.jet6_eta_         << ", ";
    std::cout << "jet6_phi_: "         << var.jet6_phi_         << ", ";
    std::cout << "jet6_energy_: "      << var.jet6_energy_      << ", ";
    std::cout << "jet6_mass_: "        << var.jet6_mass_        << ", ";
    std::cout << "jet6_btag_: "        << var.jet6_btag_        << ", ";
    std::cout << "jet7_pt_: "          << var.jet7_pt_          << ", ";
    std::cout << "jet7_eta_: "         << var.jet7_eta_         << ", ";
    std::cout << "jet7_phi_: "         << var.jet7_phi_         << ", ";
    std::cout << "jet7_energy_: "      << var.jet7_energy_      << ", ";
    std::cout << "jet7_mass_: "        << var.jet7_mass_        << ", ";
    std::cout << "jet7_btag_: "        << var.jet7_btag_        << ", ";
    std::cout << "jet8_pt_: "          << var.jet8_pt_          << ", ";
    std::cout << "jet8_eta_: "         << var.jet8_eta_         << ", ";
    std::cout << "jet8_phi_: "         << var.jet8_phi_         << ", ";
    std::cout << "jet8_energy_: "      << var.jet8_energy_      << ", ";
    std::cout << "jet8_mass_: "        << var.jet8_mass_        << ", ";
    std::cout << "jet8_btag_: "        << var.jet8_btag_        << ", ";
    std::cout << "jet9_pt_: "          << var.jet9_pt_          << ", ";
    std::cout << "jet9_eta_: "         << var.jet9_eta_         << ", ";
    std::cout << "jet9_phi_: "         << var.jet9_phi_         << ", ";
    std::cout << "jet9_energy_: "      << var.jet9_energy_      << ", ";
    std::cout << "jet9_mass_: "        << var.jet9_mass_        << ", ";
    std::cout << "jet9_btag_: "        << var.jet9_btag_        << ", ";
    std::cout << "jet10_pt_: "         << var.jet10_pt_         << ", ";
    std::cout << "jet10_eta_: "        << var.jet10_eta_        << ", ";
    std::cout << "jet10_phi_: "        << var.jet10_phi_        << ", ";
    std::cout << "jet10_energy_: "     << var.jet10_energy_     << ", ";
    std::cout << "jet10_mass_: "       << var.jet10_mass_       << ", ";
    std::cout << "jet10_btag_: "       << var.jet10_btag_       << ", ";
    std::cout << "jet11_pt_: "         << var.jet11_pt_         << ", ";
    std::cout << "jet11_eta_: "        << var.jet11_eta_        << ", ";
    std::cout << "jet11_phi_: "        << var.jet11_phi_        << ", ";
    std::cout << "jet11_energy_: "     << var.jet11_energy_     << ", ";
    std::cout << "jet11_mass_: "       << var.jet11_mass_       << ", ";
    std::cout << "jet11_btag_: "       << var.jet11_btag_       << ", ";
    std::cout << "jet12_pt_: "         << var.jet12_pt_         << ", ";
    std::cout << "jet12_eta_: "        << var.jet12_eta_        << ", ";
    std::cout << "jet12_phi_: "        << var.jet12_phi_        << ", ";
    std::cout << "jet12_energy_: "     << var.jet12_energy_     << ", ";
    std::cout << "jet12_mass_: "       << var.jet12_mass_       << ", ";
    std::cout << "jet12_btag_: "       << var.jet12_btag_       << ", ";
    std::cout << "jet13_pt_: "         << var.jet13_pt_         << ", ";
    std::cout << "jet13_eta_: "        << var.jet13_eta_        << ", ";
    std::cout << "jet13_phi_: "        << var.jet13_phi_        << ", ";
    std::cout << "jet13_energy_: "     << var.jet13_energy_     << ", ";
    std::cout << "jet13_mass_: "       << var.jet13_mass_       << ", ";
    std::cout << "jet13_btag_: "       << var.jet13_btag_       << ", ";
    std::cout << "jet14_pt_: "         << var.jet14_pt_         << ", ";
    std::cout << "jet14_eta_: "        << var.jet14_eta_        << ", ";
    std::cout << "jet14_phi_: "        << var.jet14_phi_        << ", ";
    std::cout << "jet14_energy_: "     << var.jet14_energy_     << ", ";
    std::cout << "jet14_mass_: "       << var.jet14_mass_       << ", ";
    std::cout << "jet14_btag_: "       << var.jet14_btag_       << ", ";
    std::cout << "jet15_pt_: "         << var.jet15_pt_         << ", ";
    std::cout << "jet15_eta_: "        << var.jet15_eta_        << ", ";
    std::cout << "jet15_phi_: "        << var.jet15_phi_        << ", ";
    std::cout << "jet15_energy_: "     << var.jet15_energy_     << ", ";
    std::cout << "jet15_mass_: "       << var.jet15_mass_       << ", ";
    std::cout << "jet15_btag_: "       << var.jet15_btag_       << ", ";
};
//}}}
// print details jets{{{
inline
void THQ_BDT_Helper::print_details_jets(MoreVariables var)
{
    std::cout << "jet1  pt: " << var.jet1_pt_  << ", eta: " << var.jet1_eta_  << ", phi: " << var.jet1_phi_  << ", E: " << var.jet1_energy_  << std::endl;
    std::cout << "jet2  pt: " << var.jet2_pt_  << ", eta: " << var.jet2_eta_  << ", phi: " << var.jet2_phi_  << ", E: " << var.jet2_energy_  << std::endl;
    std::cout << "jet3  pt: " << var.jet3_pt_  << ", eta: " << var.jet3_eta_  << ", phi: " << var.jet3_phi_  << ", E: " << var.jet3_energy_  << std::endl;
    std::cout << "jet4  pt: " << var.jet4_pt_  << ", eta: " << var.jet4_eta_  << ", phi: " << var.jet4_phi_  << ", E: " << var.jet4_energy_  << std::endl;
    std::cout << "jet5  pt: " << var.jet5_pt_  << ", eta: " << var.jet5_eta_  << ", phi: " << var.jet5_phi_  << ", E: " << var.jet5_energy_  << std::endl;
    std::cout << "jet6  pt: " << var.jet6_pt_  << ", eta: " << var.jet6_eta_  << ", phi: " << var.jet6_phi_  << ", E: " << var.jet6_energy_  << std::endl;
    std::cout << "jet7  pt: " << var.jet7_pt_  << ", eta: " << var.jet7_eta_  << ", phi: " << var.jet7_phi_  << ", E: " << var.jet7_energy_  << std::endl;
    std::cout << "jet8  pt: " << var.jet8_pt_  << ", eta: " << var.jet8_eta_  << ", phi: " << var.jet8_phi_  << ", E: " << var.jet8_energy_  << std::endl;
    std::cout << "jet9  pt: " << var.jet9_pt_  << ", eta: " << var.jet9_eta_  << ", phi: " << var.jet9_phi_  << ", E: " << var.jet9_energy_  << std::endl;
    std::cout << "jet10 pt: " << var.jet10_pt_ << ", eta: " << var.jet10_eta_ << ", phi: " << var.jet10_phi_ << ", E: " << var.jet10_energy_ << std::endl;
    std::cout << "jet11 pt: " << var.jet11_pt_ << ", eta: " << var.jet11_eta_ << ", phi: " << var.jet11_phi_ << ", E: " << var.jet11_energy_ << std::endl;
    std::cout << "jet12 pt: " << var.jet12_pt_ << ", eta: " << var.jet12_eta_ << ", phi: " << var.jet12_phi_ << ", E: " << var.jet12_energy_ << std::endl;
    std::cout << "jet13 pt: " << var.jet13_pt_ << ", eta: " << var.jet13_eta_ << ", phi: " << var.jet13_phi_ << ", E: " << var.jet13_energy_ << std::endl;
    std::cout << "jet14 pt: " << var.jet14_pt_ << ", eta: " << var.jet14_eta_ << ", phi: " << var.jet14_phi_ << ", E: " << var.jet14_energy_ << std::endl;
    std::cout << "jet15 pt: " << var.jet15_pt_ << ", eta: " << var.jet15_eta_ << ", phi: " << var.jet15_phi_ << ", E: " << var.jet15_energy_ << std::endl;
};
//}}}

/*
// print details more{{{
inline
void THQ_BDT_Helper::print_details_more(MoreVariables var)
{
    printf("%s: %.10f , " , "diphoton_ptOverM_" , var.diphoton_ptOverM_ );
    printf("%s: %.10f , " , "diphoton_pt_"      , var.diphoton_pt_      );
    printf("%s: %.10f , " , "diphoton_eta_"     , var.diphoton_eta_     );
    printf("%s: %.10f , " , "diphoton_phi_"     , var.diphoton_phi_     );
    printf("%s: %.10f , " , "diphoton_energy_"  , var.diphoton_energy_  );
    printf("%s: %.10f , " , "diphoton_mass_"    , var.diphoton_mass_    );
    printf("%s: %.10f , " , "wboson_ptOverM_"   , var.wboson_ptOverM_   );
    printf("%s: %.10f , " , "wboson_pt_"        , var.wboson_pt_        );
    printf("%s: %.10f , " , "wboson_eta_"       , var.wboson_eta_       );
    printf("%s: %.10f , " , "wboson_phi_"       , var.wboson_phi_       );
    printf("%s: %.10f , " , "wboson_energy_"    , var.wboson_energy_    );
    printf("%s: %.10f , " , "wboson_mass_"      , var.wboson_mass_      );
    printf("%s: %.10f , " , "top_ptOverM_"      , var.top_ptOverM_      );
    printf("%s: %.10f , " , "top_pt_"           , var.top_pt_           );
    printf("%s: %.10f , " , "top_eta_"          , var.top_eta_          );
    printf("%s: %.10f , " , "top_phi_"          , var.top_phi_          );
    printf("%s: %.10f , " , "top_energy_"       , var.top_energy_       );
    printf("%s: %.10f , " , "top_mass_"         , var.top_mass_         );
    printf("%s: %.10f , " , "tprime_ptOverM_"   , var.tprime_ptOverM_   );
    printf("%s: %.10f , " , "tprime_pt_"        , var.tprime_pt_        );
    printf("%s: %.10f , " , "tprime_eta_"       , var.tprime_eta_       );
    printf("%s: %.10f , " , "tprime_phi_"       , var.tprime_phi_       );
    printf("%s: %.10f , " , "tprime_energy_"    , var.tprime_energy_    );
    printf("%s: %.10f , " , "tprime_mass_"      , var.tprime_mass_      );
    printf("%s: %.10f , " , "bjet_ptOverM_"     , var.bjet_ptOverM_     );
    printf("%s: %.10f , " , "bjet_pt_"          , var.bjet_pt_          );
    printf("%s: %.10f , " , "bjet_eta_"         , var.bjet_eta_         );
    printf("%s: %.10f , " , "bjet_phi_"         , var.bjet_phi_         );
    printf("%s: %.10f , " , "bjet_energy_"      , var.bjet_energy_      );
    printf("%s: %.10f , " , "bjet_mass_"        , var.bjet_mass_        );
    printf("%s: %.10f , " , "bjet_btag_"        , var.bjet_btag_        );
    printf("%s: %.10f , " , "wjet1_ptOverM_"    , var.wjet1_ptOverM_    );
    printf("%s: %.10f , " , "wjet1_pt_"         , var.wjet1_pt_         );
    printf("%s: %.10f , " , "wjet1_eta_"        , var.wjet1_eta_        );
    printf("%s: %.10f , " , "wjet1_phi_"        , var.wjet1_phi_        );
    printf("%s: %.10f , " , "wjet1_energy_"     , var.wjet1_energy_     );
    printf("%s: %.10f , " , "wjet1_mass_"       , var.wjet1_mass_       );
    printf("%s: %.10f , " , "wjet1_btag_"       , var.wjet1_btag_       );
    printf("%s: %.10f , " , "wjet2_ptOverM_"    , var.wjet2_ptOverM_    );
    printf("%s: %.10f , " , "wjet2_pt_"         , var.wjet2_pt_         );
    printf("%s: %.10f , " , "wjet2_eta_"        , var.wjet2_eta_        );
    printf("%s: %.10f , " , "wjet2_phi_"        , var.wjet2_phi_        );
    printf("%s: %.10f , " , "wjet2_energy_"     , var.wjet2_energy_     );
    printf("%s: %.10f , " , "wjet2_mass_"       , var.wjet2_mass_       );
    printf("%s: %.10f , " , "wjet2_btag_"       , var.wjet2_btag_       );
    printf("%s: %.10f , " , "jet1_pt_"          , var.jet1_pt_          );
    printf("%s: %.10f , " , "jet1_eta_"         , var.jet1_eta_         );
    printf("%s: %.10f , " , "jet1_phi_"         , var.jet1_phi_         );
    printf("%s: %.10f , " , "jet1_energy_"      , var.jet1_energy_      );
    printf("%s: %.10f , " , "jet1_mass_"        , var.jet1_mass_        );
    printf("%s: %.10f , " , "jet1_btag_"        , var.jet1_btag_        );
    printf("%s: %.10f , " , "jet2_pt_"          , var.jet2_pt_          );
    printf("%s: %.10f , " , "jet2_eta_"         , var.jet2_eta_         );
    printf("%s: %.10f , " , "jet2_phi_"         , var.jet2_phi_         );
    printf("%s: %.10f , " , "jet2_energy_"      , var.jet2_energy_      );
    printf("%s: %.10f , " , "jet2_mass_"        , var.jet2_mass_        );
    printf("%s: %.10f , " , "jet2_btag_"        , var.jet2_btag_        );
    printf("%s: %.10f , " , "jet3_pt_"          , var.jet3_pt_          );
    printf("%s: %.10f , " , "jet3_eta_"         , var.jet3_eta_         );
    printf("%s: %.10f , " , "jet3_phi_"         , var.jet3_phi_         );
    printf("%s: %.10f , " , "jet3_energy_"      , var.jet3_energy_      );
    printf("%s: %.10f , " , "jet3_mass_"        , var.jet3_mass_        );
    printf("%s: %.10f , " , "jet3_btag_"        , var.jet3_btag_        );
    printf("%s: %.10f , " , "jet4_pt_"          , var.jet4_pt_          );
    printf("%s: %.10f , " , "jet4_eta_"         , var.jet4_eta_         );
    printf("%s: %.10f , " , "jet4_phi_"         , var.jet4_phi_         );
    printf("%s: %.10f , " , "jet4_energy_"      , var.jet4_energy_      );
    printf("%s: %.10f , " , "jet4_mass_"        , var.jet4_mass_        );
    printf("%s: %.10f , " , "jet4_btag_"        , var.jet4_btag_        );
    printf("%s: %.10f , " , "jet5_pt_"          , var.jet5_pt_          );
    printf("%s: %.10f , " , "jet5_eta_"         , var.jet5_eta_         );
    printf("%s: %.10f , " , "jet5_phi_"         , var.jet5_phi_         );
    printf("%s: %.10f , " , "jet5_energy_"      , var.jet5_energy_      );
    printf("%s: %.10f , " , "jet5_mass_"        , var.jet5_mass_        );
    printf("%s: %.10f , " , "jet5_btag_"        , var.jet5_btag_        );
    printf("%s: %.10f , " , "jet6_pt_"          , var.jet6_pt_          );
    printf("%s: %.10f , " , "jet6_eta_"         , var.jet6_eta_         );
    printf("%s: %.10f , " , "jet6_phi_"         , var.jet6_phi_         );
    printf("%s: %.10f , " , "jet6_energy_"      , var.jet6_energy_      );
    printf("%s: %.10f , " , "jet6_mass_"        , var.jet6_mass_        );
    printf("%s: %.10f , " , "jet6_btag_"        , var.jet6_btag_        );
    printf("%s: %.10f , " , "jet7_pt_"          , var.jet7_pt_          );
    printf("%s: %.10f , " , "jet7_eta_"         , var.jet7_eta_         );
    printf("%s: %.10f , " , "jet7_phi_"         , var.jet7_phi_         );
    printf("%s: %.10f , " , "jet7_energy_"      , var.jet7_energy_      );
    printf("%s: %.10f , " , "jet7_mass_"        , var.jet7_mass_        );
    printf("%s: %.10f , " , "jet7_btag_"        , var.jet7_btag_        );
    printf("%s: %.10f , " , "jet8_pt_"          , var.jet8_pt_          );
    printf("%s: %.10f , " , "jet8_eta_"         , var.jet8_eta_         );
    printf("%s: %.10f , " , "jet8_phi_"         , var.jet8_phi_         );
    printf("%s: %.10f , " , "jet8_energy_"      , var.jet8_energy_      );
    printf("%s: %.10f , " , "jet8_mass_"        , var.jet8_mass_        );
    printf("%s: %.10f , " , "jet8_btag_"        , var.jet8_btag_        );
    printf("%s: %.10f , " , "jet9_pt_"          , var.jet9_pt_          );
    printf("%s: %.10f , " , "jet9_eta_"         , var.jet9_eta_         );
    printf("%s: %.10f , " , "jet9_phi_"         , var.jet9_phi_         );
    printf("%s: %.10f , " , "jet9_energy_"      , var.jet9_energy_      );
    printf("%s: %.10f , " , "jet9_mass_"        , var.jet9_mass_        );
    printf("%s: %.10f , " , "jet9_btag_"        , var.jet9_btag_        );
    printf("%s: %.10f , " , "jet10_pt_"         , var.jet10_pt_         );
    printf("%s: %.10f , " , "jet10_eta_"        , var.jet10_eta_        );
    printf("%s: %.10f , " , "jet10_phi_"        , var.jet10_phi_        );
    printf("%s: %.10f , " , "jet10_energy_"     , var.jet10_energy_     );
    printf("%s: %.10f , " , "jet10_mass_"       , var.jet10_mass_       );
    printf("%s: %.10f , " , "jet10_btag_"       , var.jet10_btag_       );
    printf("%s: %.10f , " , "jet11_pt_"         , var.jet11_pt_         );
    printf("%s: %.10f , " , "jet11_eta_"        , var.jet11_eta_        );
    printf("%s: %.10f , " , "jet11_phi_"        , var.jet11_phi_        );
    printf("%s: %.10f , " , "jet11_energy_"     , var.jet11_energy_     );
    printf("%s: %.10f , " , "jet11_mass_"       , var.jet11_mass_       );
    printf("%s: %.10f , " , "jet11_btag_"       , var.jet11_btag_       );
    printf("%s: %.10f , " , "jet12_pt_"         , var.jet12_pt_         );
    printf("%s: %.10f , " , "jet12_eta_"        , var.jet12_eta_        );
    printf("%s: %.10f , " , "jet12_phi_"        , var.jet12_phi_        );
    printf("%s: %.10f , " , "jet12_energy_"     , var.jet12_energy_     );
    printf("%s: %.10f , " , "jet12_mass_"       , var.jet12_mass_       );
    printf("%s: %.10f , " , "jet12_btag_"       , var.jet12_btag_       );
    printf("%s: %.10f , " , "jet13_pt_"         , var.jet13_pt_         );
    printf("%s: %.10f , " , "jet13_eta_"        , var.jet13_eta_        );
    printf("%s: %.10f , " , "jet13_phi_"        , var.jet13_phi_        );
    printf("%s: %.10f , " , "jet13_energy_"     , var.jet13_energy_     );
    printf("%s: %.10f , " , "jet13_mass_"       , var.jet13_mass_       );
    printf("%s: %.10f , " , "jet13_btag_"       , var.jet13_btag_       );
    printf("%s: %.10f , " , "jet14_pt_"         , var.jet14_pt_         );
    printf("%s: %.10f , " , "jet14_eta_"        , var.jet14_eta_        );
    printf("%s: %.10f , " , "jet14_phi_"        , var.jet14_phi_        );
    printf("%s: %.10f , " , "jet14_energy_"     , var.jet14_energy_     );
    printf("%s: %.10f , " , "jet14_mass_"       , var.jet14_mass_       );
    printf("%s: %.10f , " , "jet14_btag_"       , var.jet14_btag_       );
    printf("%s: %.10f , " , "jet15_pt_"         , var.jet15_pt_         );
    printf("%s: %.10f , " , "jet15_eta_"        , var.jet15_eta_        );
    printf("%s: %.10f , " , "jet15_phi_"        , var.jet15_phi_        );
    printf("%s: %.10f , " , "jet15_energy_"     , var.jet15_energy_     );
    printf("%s: %.10f , " , "jet15_mass_"       , var.jet15_mass_       );
    printf("%s: %.10f , " , "jet15_btag_"       , var.jet15_btag_       );
    printf("%s: %.10f , " , "jet15_btag_"       , var.jet15_btag_       );
};
//}}}
*/

} // end of namespace

#endif // _THQ_BDT_HELPER_
