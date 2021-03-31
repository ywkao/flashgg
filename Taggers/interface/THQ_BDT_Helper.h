#ifndef _THQ_BDT_HELPER_
#define _THQ_BDT_HELPER_

#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "flashgg/Taggers/interface/BDT_self_defined_objects.h"

namespace flashgg {

//class BDT_Hadronic_tprime {{{
class BDT_Hadronic_tprime {
public:
  BDT_Hadronic_tprime(){};
  BDT_Hadronic_tprime(std::shared_ptr<BDT_Photon> _p1,std::shared_ptr<BDT_Photon> _p2,std::shared_ptr<BDT_Jet> _b,std::shared_ptr<BDT_Jet> _j1,std::shared_ptr<BDT_Jet> _j2)
  {
    leading_photon    = _p1;
    subleading_photon = _p2;

    bjet  = _b;
    wjet1 = _j1;
    wjet2 = _j2;

    *diphoton  =  *(dynamic_cast<BDT_ptvec*>(leading_photon.get())) + *(dynamic_cast<BDT_ptvec*>(subleading_photon.get()));
    *wboson    =  *(dynamic_cast<BDT_ptvec*>(wjet1.get())) + *(dynamic_cast<BDT_ptvec*>(wjet2.get()));
    *top       =  *wboson + *(dynamic_cast<BDT_ptvec*>(bjet.get()));
    *tprime    =  *top + *diphoton;
  }

  std::shared_ptr<BDT_ptvec> diphoton = std::make_shared<BDT_ptvec>(0,0,0,0);
  std::shared_ptr<BDT_ptvec> wboson   = std::make_shared<BDT_ptvec>(0,0,0,0);
  std::shared_ptr<BDT_ptvec> top      = std::make_shared<BDT_ptvec>(0,0,0,0);
  std::shared_ptr<BDT_ptvec> tprime   = std::make_shared<BDT_ptvec>(0,0,0,0);
  std::shared_ptr<BDT_Photon> leading_photon    = nullptr;
  std::shared_ptr<BDT_Photon> subleading_photon = nullptr;
  std::shared_ptr<BDT_Jet> bjet  = nullptr;
  std::shared_ptr<BDT_Jet> wjet1 = nullptr;
  std::shared_ptr<BDT_Jet> wjet2 = nullptr;
};
//}}}
//class THQ_BDT_Helper {{{
class THQ_BDT_Helper {
public:
  THQ_BDT_Helper(std::string weight_file_name){ Init(weight_file_name); };
  ~THQ_BDT_Helper(){ clear(); };

  void addJet(float pt, float eta, float phi, float mass, float deepcsv)  { jets.push_back(std::make_shared<BDT_Jet>(pt,eta,phi,mass,deepcsv));      } ;
  void addPhoton(float pt, float eta, float phi, float mass, float IDMVA) { photons.push_back(std::make_shared<BDT_Photon>(pt,eta,phi,mass,IDMVA)); } ;
  void addPhotonPixelSeed(float a, float b)                               { leadPSV_ = a; subleadPSV_ = b; } ;
  void addNbjets(float n)                                                 { nbjets_  = n;   } ;
  void addHt(float ht)                                                    { ht_      = ht;  } ;
  void addMet(float met)                                                  { met_     = met; } ;

  void Init(std::string weight_file_name);
  void clear();
  void setDebug(bool val){debug = val;};
  double chi2_calculator_2x2(double w_mass, double t_mass);
  double convert_tmva_to_prob(double score);
  double DeltaR(const BDT_ptvec* obj1, const BDT_ptvec* obj2);
  double Rapidity(const BDT_ptvec* obj);
  double Helicity(const BDT_ptvec* particle_1, const BDT_ptvec* particle_2);
  double EvalMVA();

private:
  std::vector<std::shared_ptr<BDT_Jet>> jets;
  std::vector<std::shared_ptr<BDT_Jet>> bjets; // for b-tag scores sorting only
  std::vector<std::shared_ptr<BDT_Photon>> photons;
  std::shared_ptr<TMVA::Reader> reader = nullptr;

  double output = -999;
  double EvalScore(const std::shared_ptr<BDT_Hadronic_tprime>);

  bool debug = false;
  bool has_resonable_reco;
  float min_chi2_value = 99999.;

  float maxIDMVA_                 = -999;
  float minIDMVA_                 = -999;
  float max1_btag_                = -999;
  float max2_btag_                = -999;
  float dipho_delta_R             = -999;
  float njets_                    = -999;
  float nbjets_                   = -999;
  float ht_                       = -999;
  float leadptoM_                 = -999;
  float subleadptoM_              = -999;
  float lead_eta_                 = -999;
  float sublead_eta_              = -999;
  float jet1_ptOverM_             = -999;
  float jet1_eta_                 = -999;
  float jet1_btag_                = -999;
  float jet2_ptOverM_             = -999;
  float jet2_eta_                 = -999;
  float jet2_btag_                = -999;
  float jet3_ptOverM_             = -999;
  float jet3_eta_                 = -999;
  float jet3_btag_                = -999;
  float jet4_ptOverM_             = -999;
  float jet4_eta_                 = -999;
  float jet4_btag_                = -999;
  float leadPSV_                  = -999;
  float subleadPSV_               = -999;
  float dipho_cosphi_             = -999;
  float dipho_rapidity_           = -999;
  float met_                      = -999;
  float dipho_pt_over_mass_       = -999;
  float helicity_angle_           = -999;
  float chi2_value_               = -999;
  float chi2_bjet_ptOverM_        = -999;
  float chi2_bjet_eta_            = -999;
  float chi2_bjet_btagScores_     = -999;
  float chi2_wjet1_ptOverM_       = -999;
  float chi2_wjet1_eta_           = -999;
  float chi2_wjet1_btagScores_    = -999;
  float chi2_wjet2_ptOverM_       = -999;
  float chi2_wjet2_eta_           = -999;
  float chi2_wjet2_btagScores_    = -999;
  float chi2_wjets_deltaR_        = -999;
  float chi2_wboson_ptOverM_      = -999;
  float chi2_wboson_eta_          = -999;
  float chi2_wboson_mass_         = -999;
  float chi2_wboson_deltaR_bjet_  = -999;
  float chi2_tbw_mass_            = -999;
  float chi2_tbw_ptOverM_         = -999;
  float chi2_tbw_eta_             = -999;
  float chi2_tbw_deltaR_dipho_    = -999;
  float chi2_tprime_ptOverM_      = -999;
  float chi2_tprime_eta_          = -999;
  float chi2_tprime_deltaR_tbw_   = -999;
  float chi2_tprime_deltaR_dipho_ = -999;
  float tprime_pt_ratio_          = -999;
};
//}}}

//void THQ_BDT_Helper::Init(std::string weight_file_name){{{
inline
void THQ_BDT_Helper::Init(std::string weight_file_name){

  reader = std::make_shared<TMVA::Reader>( "!Color:!Silent" );

  reader->AddVariable("maxIDMVA_"                 , &maxIDMVA_                 );
  reader->AddVariable("minIDMVA_"                 , &minIDMVA_                 );
  reader->AddVariable("max1_btag_"                , &max1_btag_                );
  reader->AddVariable("max2_btag_"                , &max2_btag_                );
  reader->AddVariable("dipho_delta_R"             , &dipho_delta_R             );
  reader->AddVariable("njets_"                    , &njets_                    );
  reader->AddVariable("nbjets_"                   , &nbjets_                   );
  reader->AddVariable("ht_"                       , &ht_                       );
  reader->AddVariable("leadptoM_"                 , &leadptoM_                 );
  reader->AddVariable("subleadptoM_"              , &subleadptoM_              );
  reader->AddVariable("lead_eta_"                 , &lead_eta_                 );
  reader->AddVariable("sublead_eta_"              , &sublead_eta_              );
  reader->AddVariable("jet1_ptOverM_"             , &jet1_ptOverM_             );
  reader->AddVariable("jet1_eta_"                 , &jet1_eta_                 );
  reader->AddVariable("jet1_btag_"                , &jet1_btag_                );
  reader->AddVariable("jet2_ptOverM_"             , &jet2_ptOverM_             );
  reader->AddVariable("jet2_eta_"                 , &jet2_eta_                 );
  reader->AddVariable("jet2_btag_"                , &jet2_btag_                );
  reader->AddVariable("jet3_ptOverM_"             , &jet3_ptOverM_             );
  reader->AddVariable("jet3_eta_"                 , &jet3_eta_                 );
  reader->AddVariable("jet3_btag_"                , &jet3_btag_                );
  reader->AddVariable("jet4_ptOverM_"             , &jet4_ptOverM_             );
  reader->AddVariable("jet4_eta_"                 , &jet4_eta_                 );
  reader->AddVariable("jet4_btag_"                , &jet4_btag_                );
  reader->AddVariable("leadPSV_"                  , &leadPSV_                  );
  reader->AddVariable("subleadPSV_"               , &subleadPSV_               );
  reader->AddVariable("dipho_cosphi_"             , &dipho_cosphi_             );
  reader->AddVariable("dipho_rapidity_"           , &dipho_rapidity_           );
  reader->AddVariable("met_"                      , &met_                      );
  reader->AddVariable("dipho_pt_over_mass_"       , &dipho_pt_over_mass_       );
  reader->AddVariable("helicity_angle_"           , &helicity_angle_           );
  reader->AddVariable("chi2_value_"               , &chi2_value_               );
  reader->AddVariable("chi2_bjet_ptOverM_"        , &chi2_bjet_ptOverM_        );
  reader->AddVariable("chi2_bjet_eta_"            , &chi2_bjet_eta_            );
  reader->AddVariable("chi2_bjet_btagScores_"     , &chi2_bjet_btagScores_     );
  reader->AddVariable("chi2_wjet1_ptOverM_"       , &chi2_wjet1_ptOverM_       );
  reader->AddVariable("chi2_wjet1_eta_"           , &chi2_wjet1_eta_           );
  reader->AddVariable("chi2_wjet1_btagScores_"    , &chi2_wjet1_btagScores_    );
  reader->AddVariable("chi2_wjet2_ptOverM_"       , &chi2_wjet2_ptOverM_       );
  reader->AddVariable("chi2_wjet2_eta_"           , &chi2_wjet2_eta_           );
  reader->AddVariable("chi2_wjet2_btagScores_"    , &chi2_wjet2_btagScores_    );
  reader->AddVariable("chi2_wjets_deltaR_"        , &chi2_wjets_deltaR_        );
  reader->AddVariable("chi2_wboson_ptOverM_"      , &chi2_wboson_ptOverM_      );
  reader->AddVariable("chi2_wboson_eta_"          , &chi2_wboson_eta_          );
  reader->AddVariable("chi2_wboson_mass_"         , &chi2_wboson_mass_         );
  reader->AddVariable("chi2_wboson_deltaR_bjet_"  , &chi2_wboson_deltaR_bjet_  );
  reader->AddVariable("chi2_tbw_mass_"            , &chi2_tbw_mass_            );
  reader->AddVariable("chi2_tbw_ptOverM_"         , &chi2_tbw_ptOverM_         );
  reader->AddVariable("chi2_tbw_eta_"             , &chi2_tbw_eta_             );
  reader->AddVariable("chi2_tbw_deltaR_dipho_"    , &chi2_tbw_deltaR_dipho_    );
  reader->AddVariable("chi2_tprime_ptOverM_"      , &chi2_tprime_ptOverM_      );
  reader->AddVariable("chi2_tprime_eta_"          , &chi2_tprime_eta_          );
  reader->AddVariable("chi2_tprime_deltaR_tbw_"   , &chi2_tprime_deltaR_tbw_   );
  reader->AddVariable("chi2_tprime_deltaR_dipho_" , &chi2_tprime_deltaR_dipho_ );
  reader->AddVariable("tprime_pt_ratio_"          , &tprime_pt_ratio_          );

  reader->BookMVA("BDT", weight_file_name);
};
//}}}
//void THQ_BDT_Helper::clear(){{{
inline
void THQ_BDT_Helper::clear(){

  jets.clear();
  bjets.clear();
  photons.clear();

  maxIDMVA_                 = -999;
  minIDMVA_                 = -999;
  max1_btag_                = -999;
  max2_btag_                = -999;
  dipho_delta_R             = -999;
  njets_                    = -999;
  nbjets_                   = -999;
  ht_                       = -999;
  leadptoM_                 = -999;
  subleadptoM_              = -999;
  lead_eta_                 = -999;
  sublead_eta_              = -999;
  jet1_ptOverM_             = -999;
  jet1_eta_                 = -999;
  jet1_btag_                = -999;
  jet2_ptOverM_             = -999;
  jet2_eta_                 = -999;
  jet2_btag_                = -999;
  jet3_ptOverM_             = -999;
  jet3_eta_                 = -999;
  jet3_btag_                = -999;
  jet4_ptOverM_             = -999;
  jet4_eta_                 = -999;
  jet4_btag_                = -999;
  leadPSV_                  = -999;
  subleadPSV_               = -999;
  dipho_cosphi_             = -999;
  dipho_rapidity_           = -999;
  met_                      = -999;
  dipho_pt_over_mass_       = -999;
  helicity_angle_           = -999;
  chi2_value_               = -999;
  chi2_bjet_ptOverM_        = -999;
  chi2_bjet_eta_            = -999;
  chi2_bjet_btagScores_     = -999;
  chi2_wjet1_ptOverM_       = -999;
  chi2_wjet1_eta_           = -999;
  chi2_wjet1_btagScores_    = -999;
  chi2_wjet2_ptOverM_       = -999;
  chi2_wjet2_eta_           = -999;
  chi2_wjet2_btagScores_    = -999;
  chi2_wjets_deltaR_        = -999;
  chi2_wboson_ptOverM_      = -999;
  chi2_wboson_eta_          = -999;
  chi2_wboson_mass_         = -999;
  chi2_wboson_deltaR_bjet_  = -999;
  chi2_tbw_mass_            = -999;
  chi2_tbw_ptOverM_         = -999;
  chi2_tbw_eta_             = -999;
  chi2_tbw_deltaR_dipho_    = -999;
  chi2_tprime_ptOverM_      = -999;
  chi2_tprime_eta_          = -999;
  chi2_tprime_deltaR_tbw_   = -999;
  chi2_tprime_deltaR_dipho_ = -999;
  tprime_pt_ratio_          = -999;
};
//}}}
//double THQ_BDT_Helper::chi2_calculator_2x2(double w_mass, double t_mass){{{
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
//}}}
//double THQ_BDT_Helper::convert_tmva_to_prob(double score) {{{
inline
double THQ_BDT_Helper::convert_tmva_to_prob(double score)
{
  // Undo TMVA transformation
  double raw_score = -0.5 * log( (2 / (score + 1)) - 1);

  // Apply logistic (sigmoid) transformation
  double prob = 1 / (1 + exp(-raw_score));
  return prob;
};
//}}}
//double THQ_BDT_Helper::DeltaR(const BDT_ptvec* obj1, const BDT_ptvec* obj2){{{
inline
double THQ_BDT_Helper::DeltaR(const BDT_ptvec* obj1, const BDT_ptvec* obj2)
{
    double delta_R = sqrt( pow((obj1->Eta() - obj2->Eta()), 2) + pow((obj1->Phi() - obj2->Phi()), 2) );
    return delta_R;
};
//}}}
//double THQ_BDT_Helper::Rapidity(const BDT_ptvec* obj){{{
inline
double THQ_BDT_Helper::Rapidity(const BDT_ptvec* obj)
{
    double rapidity = 0.5*log( (obj->E()+obj->Pz()) / (obj->E()-obj->Pz()) );
    return rapidity;
};
//}}}
//double THQ_BDT_Helper::Helicity(const BDT_ptvec* particle_1, const BDT_ptvec* particle_2){{{
inline
double THQ_BDT_Helper::Helicity(const BDT_ptvec* particle_1, const BDT_ptvec* particle_2)
{
  std::shared_ptr<BDT_ptvec> p1_     = std::make_shared<BDT_ptvec>(0,0,0,0);
  std::shared_ptr<BDT_ptvec> parent_ = std::make_shared<BDT_ptvec>(0,0,0,0);

  *p1_ = *particle_1;
  *parent_ = *particle_1 + *particle_2;

  TLorentzVector p1;
  TLorentzVector parent;
  p1.SetPtEtaPhiM(p1_->Pt(), p1_->Eta(), p1_->Phi(), p1_->M());
  parent.SetPtEtaPhiM(parent_->Pt(), parent_->Eta(), parent_->Phi(), parent_->M());
  
  TVector3 boost_to_parent = -(parent.BoostVector());
  p1.Boost(boost_to_parent);

  TVector3 v1 = p1.Vect();
  TVector3 vParent = parent.Vect();
  
  double cos_theta_1 = (v1.Dot(vParent)) / (v1.Mag() * vParent.Mag());

  if(debug){
      printf("[check] v1.Mag() = %.2f\n"        , v1.Mag());
      printf("[check] vParent.Mag() = %.2f\n"   , vParent.Mag());
      printf("[check] v1.Dot(vParent) = %.2f\n" , v1.Dot(vParent));
      printf("[check] cos_theta_1 = %.2f\n"     , cos_theta_1);
  }

  return abs(cos_theta_1);  
};
//}}}
//double THQ_BDT_Helper::EvalMVA() {{{
inline
double THQ_BDT_Helper::EvalMVA(){

  double pfDeepCSV_btag_loose_wp = 0.1522;

  min_chi2_value = 99999.;
  std::vector<int> indices_bjj(3, -1);

  //------------------------------ Minimum chi-2 method (cov.) ------------------------------//
  std::size_t num_jets = jets.size();
  for(std::size_t i = 0; i < num_jets; ++i ){ // b-jet
      if (jets[i]->deepcsv() < pfDeepCSV_btag_loose_wp) continue;
      for(std::size_t j = 0; j < num_jets-1; ++j ){ // w-jet1
          if(j==i) continue;
          for(std::size_t k = j+1; k < num_jets; ++k ){ // w-jet2
              if(k==i) continue;
  
              std::shared_ptr<BDT_ptvec> bjet          = std::make_shared<BDT_ptvec>(0,0,0,0);
              std::shared_ptr<BDT_ptvec> w_candidate   = std::make_shared<BDT_ptvec>(0,0,0,0);
              std::shared_ptr<BDT_ptvec> top_candidate = std::make_shared<BDT_ptvec>(0,0,0,0);

              *bjet = *dynamic_cast<BDT_ptvec*>(jets[i].get());
              *w_candidate = *dynamic_cast<BDT_ptvec*>(jets[j].get()) + *dynamic_cast<BDT_ptvec*>(jets[k].get());
              *top_candidate = *w_candidate + *bjet;

              double w_mass = w_candidate->M();
              double t_mass = top_candidate->M();
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

  //------------------------------ Evaluation ------------------------------//
  has_resonable_reco = indices_bjj[0] >= 0 && indices_bjj[1] >= 0 && indices_bjj[2] >= 0;
  if(has_resonable_reco)
  {
      auto topcand = std::make_shared<BDT_Hadronic_tprime>(photons[0], photons[1], jets[indices_bjj[0]], jets[indices_bjj[1]], jets[indices_bjj[2]]);
      output = EvalScore(topcand);
      if (debug) std::cout << "returning " << output << std::endl;
      return output;
  } else {
      if (debug) std::cout << "returning " << output << std::endl;
      return output;
  }

};
//}}}


inline
double THQ_BDT_Helper::EvalScore(const std::shared_ptr<BDT_Hadronic_tprime> object){
  bool pass_eta_criteria_on_wjets = ( object->wjet1->Eta() < 3. && object->wjet2->Eta() < 3. );

  // sorting for getting max-btag/2nd-max-btag
  bjets = jets;
  std::sort( bjets.begin(), bjets.end(), [](const std::shared_ptr<BDT_Jet> &j1, const std::shared_ptr<BDT_Jet> &j2){return j1->deepcsv() > j2->deepcsv();} );

  // set values before evaluation
  maxIDMVA_            =  photons[0]->IDMVA() > photons[1]->IDMVA() ? photons[0]->IDMVA() : photons[1]->IDMVA();
  minIDMVA_            =  photons[0]->IDMVA() < photons[1]->IDMVA() ? photons[0]->IDMVA() : photons[1]->IDMVA();
  leadptoM_            =  photons[0]->Pt() / object->diphoton->M();
  subleadptoM_         =  photons[1]->Pt() / object->diphoton->M();
  lead_eta_            =  photons[0]->Eta();
  sublead_eta_         =  photons[1]->Eta();
  dipho_pt_over_mass_  =  object->diphoton->Pt() / object->diphoton->M();
  dipho_cosphi_        =  TMath::Cos(photons[0]->Phi() - photons[1]->Phi());
  dipho_rapidity_      =  Rapidity(dynamic_cast<BDT_ptvec*>(object->diphoton.get()));
  dipho_delta_R        =  DeltaR  (dynamic_cast<BDT_ptvec*>(photons[0].get()), dynamic_cast<BDT_ptvec*>(photons[1].get()));
  helicity_angle_      =  Helicity(dynamic_cast<BDT_ptvec*>(photons[0].get()), dynamic_cast<BDT_ptvec*>(photons[1].get()));

  bool reco_condition = has_resonable_reco && pass_eta_criteria_on_wjets;
  njets_         =  jets.size();
  max1_btag_     =  bjets[0]->deepcsv();
  max2_btag_     =  bjets[1]->deepcsv();
  jet1_eta_      =  (njets_ >= 1) ? jets[0]->Eta()     : -999;
  jet2_eta_      =  (njets_ >= 2) ? jets[1]->Eta()     : -999;
  jet3_eta_      =  (njets_ >= 3) ? jets[2]->Eta()     : -999;
  jet4_eta_      =  (njets_ >= 4) ? jets[3]->Eta()     : -999;
  jet1_btag_     =  (njets_ >= 1) ? jets[0]->deepcsv() : -999;
  jet2_btag_     =  (njets_ >= 2) ? jets[1]->deepcsv() : -999;
  jet3_btag_     =  (njets_ >= 3) ? jets[2]->deepcsv() : -999;
  jet4_btag_     =  (njets_ >= 4) ? jets[3]->deepcsv() : -999;
  jet1_ptOverM_  =  (njets_ >= 1 && reco_condition) ? jets[0]->Pt() / object->top->M()  : -999;
  jet2_ptOverM_  =  (njets_ >= 2 && reco_condition) ? jets[1]->Pt() / object->top->M()  : -999;
  jet3_ptOverM_  =  (njets_ >= 3 && reco_condition) ? jets[2]->Pt() / object->top->M()  : -999;
  jet4_ptOverM_  =  (njets_ >= 4 && reco_condition) ? jets[3]->Pt() / object->top->M()  : -999;

  chi2_value_               = min_chi2_value;
  chi2_bjet_btagScores_     = reco_condition ? object->bjet->deepcsv()                           : -999;
  chi2_wjet1_btagScores_    = reco_condition ? object->wjet1->deepcsv()                          : -999;
  chi2_wjet2_btagScores_    = reco_condition ? object->wjet2->deepcsv()                          : -999;
  chi2_bjet_ptOverM_        = reco_condition ? object->bjet->Pt() / object->top->M()             : -999;
  chi2_wjet1_ptOverM_       = reco_condition ? object->wjet1->Pt() / object->wboson->M()         : -999;
  chi2_wjet2_ptOverM_       = reco_condition ? object->wjet2->Pt() / object->wboson->M()         : -999;
  chi2_wboson_ptOverM_      = reco_condition ? object->wboson->Pt() / object->wboson->M()        : -999;
  chi2_tbw_ptOverM_         = reco_condition ? object->top->Pt() / object->top->M()              : -999;
  chi2_tprime_ptOverM_      = reco_condition ? object->tprime->Pt()/object->tprime->M()          : -999;
  chi2_bjet_eta_            = reco_condition ? object->bjet->Eta()                               : -999;
  chi2_wjet1_eta_           = reco_condition ? object->wjet1->Eta()                              : -999;
  chi2_wjet2_eta_           = reco_condition ? object->wjet2->Eta()                              : -999;
  chi2_wboson_eta_          = reco_condition ? object->wboson->Eta()                             : -999;
  chi2_tbw_eta_             = reco_condition ? object->top->Eta()                                : -999;
  chi2_tprime_eta_          = reco_condition ? object->tprime->Eta()                             : -999;
  chi2_wboson_mass_         = reco_condition ? object->wboson->M()                               : -999;
  chi2_tbw_mass_            = reco_condition ? object->top->M()                                  : -999;
  tprime_pt_ratio_          = reco_condition ? (object->top->Pt() + object->diphoton->Pt())/ ht_ : -999;
  chi2_wboson_deltaR_bjet_  = reco_condition ? DeltaR(dynamic_cast<BDT_ptvec*>(object->wboson.get()) , dynamic_cast<BDT_ptvec*>(object->bjet.get()))     : -999;
  chi2_tbw_deltaR_dipho_    = reco_condition ? DeltaR(dynamic_cast<BDT_ptvec*>(object->top.get())    , dynamic_cast<BDT_ptvec*>(object->diphoton.get())) : -999;
  chi2_tprime_deltaR_tbw_   = reco_condition ? DeltaR(dynamic_cast<BDT_ptvec*>(object->tprime.get()) , dynamic_cast<BDT_ptvec*>(object->top.get()))      : -999;
  chi2_tprime_deltaR_dipho_ = reco_condition ? DeltaR(dynamic_cast<BDT_ptvec*>(object->tprime.get()) , dynamic_cast<BDT_ptvec*>(object->diphoton.get())) : -999;
  chi2_wjets_deltaR_        = reco_condition ? DeltaR(dynamic_cast<BDT_ptvec*>(object->wjet1.get())  , dynamic_cast<BDT_ptvec*>(object->wjet2.get()))    : -999;

  double score = convert_tmva_to_prob( reader->EvaluateMVA("BDT") );

  if (debug) {
      printf("%-20s: %.2f\n" , "maxIDMVA_"                 , maxIDMVA_                 );
      printf("%-20s: %.2f\n" , "minIDMVA_"                 , minIDMVA_                 );
      printf("%-20s: %.2f\n" , "max1_btag_"                , max1_btag_                );
      printf("%-20s: %.2f\n" , "max2_btag_"                , max2_btag_                );
      printf("%-20s: %.2f\n" , "dipho_delta_R"             , dipho_delta_R             );
      printf("%-20s: %.2f\n" , "njets_"                    , njets_                    );
      printf("%-20s: %.2f\n" , "nbjets_"                   , nbjets_                   );
      printf("%-20s: %.2f\n" , "ht_"                       , ht_                       );
      printf("%-20s: %.2f\n" , "leadptoM_"                 , leadptoM_                 );
      printf("%-20s: %.2f\n" , "subleadptoM_"              , subleadptoM_              );
      printf("%-20s: %.2f\n" , "lead_eta_"                 , lead_eta_                 );
      printf("%-20s: %.2f\n" , "sublead_eta_"              , sublead_eta_              );
      printf("%-20s: %.2f\n" , "jet1_ptOverM_"             , jet1_ptOverM_             );
      printf("%-20s: %.2f\n" , "jet1_eta_"                 , jet1_eta_                 );
      printf("%-20s: %.2f\n" , "jet1_btag_"                , jet1_btag_                );
      printf("%-20s: %.2f\n" , "jet2_ptOverM_"             , jet2_ptOverM_             );
      printf("%-20s: %.2f\n" , "jet2_eta_"                 , jet2_eta_                 );
      printf("%-20s: %.2f\n" , "jet2_btag_"                , jet2_btag_                );
      printf("%-20s: %.2f\n" , "jet3_ptOverM_"             , jet3_ptOverM_             );
      printf("%-20s: %.2f\n" , "jet3_eta_"                 , jet3_eta_                 );
      printf("%-20s: %.2f\n" , "jet3_btag_"                , jet3_btag_                );
      printf("%-20s: %.2f\n" , "jet4_ptOverM_"             , jet4_ptOverM_             );
      printf("%-20s: %.2f\n" , "jet4_eta_"                 , jet4_eta_                 );
      printf("%-20s: %.2f\n" , "jet4_btag_"                , jet4_btag_                );
      printf("%-20s: %.2f\n" , "leadPSV_"                  , leadPSV_                  );
      printf("%-20s: %.2f\n" , "subleadPSV_"               , subleadPSV_               );
      printf("%-20s: %.2f\n" , "dipho_cosphi_"             , dipho_cosphi_             );
      printf("%-20s: %.2f\n" , "dipho_rapidity_"           , dipho_rapidity_           );
      printf("%-20s: %.2f\n" , "met_"                      , met_                      );
      printf("%-20s: %.2f\n" , "dipho_pt_over_mass_"       , dipho_pt_over_mass_       );
      printf("%-20s: %.2f\n" , "helicity_angle_"           , helicity_angle_           );
      printf("%-20s: %.2f\n" , "chi2_value_"               , chi2_value_               );
      printf("%-20s: %.2f\n" , "chi2_bjet_ptOverM_"        , chi2_bjet_ptOverM_        );
      printf("%-20s: %.2f\n" , "chi2_bjet_eta_"            , chi2_bjet_eta_            );
      printf("%-20s: %.2f\n" , "chi2_bjet_btagScores_"     , chi2_bjet_btagScores_     );
      printf("%-20s: %.2f\n" , "chi2_wjet1_ptOverM_"       , chi2_wjet1_ptOverM_       );
      printf("%-20s: %.2f\n" , "chi2_wjet1_eta_"           , chi2_wjet1_eta_           );
      printf("%-20s: %.2f\n" , "chi2_wjet1_btagScores_"    , chi2_wjet1_btagScores_    );
      printf("%-20s: %.2f\n" , "chi2_wjet2_ptOverM_"       , chi2_wjet2_ptOverM_       );
      printf("%-20s: %.2f\n" , "chi2_wjet2_eta_"           , chi2_wjet2_eta_           );
      printf("%-20s: %.2f\n" , "chi2_wjet2_btagScores_"    , chi2_wjet2_btagScores_    );
      printf("%-20s: %.2f\n" , "chi2_wjets_deltaR_"        , chi2_wjets_deltaR_        );
      printf("%-20s: %.2f\n" , "chi2_wboson_ptOverM_"      , chi2_wboson_ptOverM_      );
      printf("%-20s: %.2f\n" , "chi2_wboson_eta_"          , chi2_wboson_eta_          );
      printf("%-20s: %.2f\n" , "chi2_wboson_mass_"         , chi2_wboson_mass_         );
      printf("%-20s: %.2f\n" , "chi2_wboson_deltaR_bjet_"  , chi2_wboson_deltaR_bjet_  );
      printf("%-20s: %.2f\n" , "chi2_tbw_mass_"            , chi2_tbw_mass_            );
      printf("%-20s: %.2f\n" , "chi2_tbw_ptOverM_"         , chi2_tbw_ptOverM_         );
      printf("%-20s: %.2f\n" , "chi2_tbw_eta_"             , chi2_tbw_eta_             );
      printf("%-20s: %.2f\n" , "chi2_tbw_deltaR_dipho_"    , chi2_tbw_deltaR_dipho_    );
      printf("%-20s: %.2f\n" , "chi2_tprime_ptOverM_"      , chi2_tprime_ptOverM_      );
      printf("%-20s: %.2f\n" , "chi2_tprime_eta_"          , chi2_tprime_eta_          );
      printf("%-20s: %.2f\n" , "chi2_tprime_deltaR_tbw_"   , chi2_tprime_deltaR_tbw_   );
      printf("%-20s: %.2f\n" , "chi2_tprime_deltaR_dipho_" , chi2_tprime_deltaR_dipho_ );
      printf("%-20s: %.2f\n" , "tprime_pt_ratio_"          , tprime_pt_ratio_          );

      std::cout << score << std::endl;
  }

  return score;
};

}

#endif // _THQ_BDT_HELPER_
