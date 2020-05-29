#ifdef _TOP_RECO_HELPER_H_
#define _TOP_RECO_HELPER_H_

//#chi-2 related
//--------------------------------------------------//
//get_bjet_indices{{{
double pfDeepCSVJetTags_loose  = 0.1522;
vector<int> get_bjet_indices(vector<TLorentzVector> Jets, vector<double> btag_scores)
{
    vector<int> indices;
    //# consider all b-jet candidates
    //for(size_t i=0; i!=Jets.size(); ++i)
    //{
    //    double btag_score = btag_scores[i];
    //    if(btag_score >= pfDeepCSVJetTags_loose) indices.push_back(i);
    //}

    //# consider the jet with the highest b-tag score only
    int index_bjet = std::max_element(btag_scores.begin(), btag_scores.end()) - btag_scores.begin();
    indices = {index_bjet};

    return indices;
}
//}}}
//chi-2 with 3x3 covariance matrix{{{
std::vector<int> get_bjjq_indices_chi2_3x3(std::vector<TLorentzVector> Jets, int index_bjet, TLorentzVector diphoton, double &chi2_min)
{
    // 1. pick up 3 jets
    // 2. q-jj three combinations 
    int index_tqh_qjet;
    std::vector<int> index_wjets(2);
    TLorentzVector wjets[2], tqh_qjet;
    TLorentzVector w_boson, sm_top, fcnc_top; 
    TLorentzVector bjet = Jets[index_bjet];
    std::size_t num_jets = Jets.size();
    for(std::size_t i=0; i!=num_jets; ++i){
        if(i==index_bjet) continue;//bypass bjet
        for(std::size_t j=i+1; j!=num_jets; ++j){
            if(j==index_bjet) continue;//bypass bjet
            for(std::size_t k=j+1; k!=num_jets; ++k){
                if(k==index_bjet) continue;//bypass bjet
                //--- combinations{{{
                TLorentzVector jets_chosen[3];
                jets_chosen[0] = Jets[i];
                jets_chosen[1] = Jets[j];
                jets_chosen[2] = Jets[k];

                TLorentzVector tqh_q_chosen[3];
                tqh_q_chosen[0] = jets_chosen[0];
                tqh_q_chosen[1] = jets_chosen[1];
                tqh_q_chosen[2] = jets_chosen[2];

                TLorentzVector w_candidate[3];
                w_candidate[0] = jets_chosen[1] + jets_chosen[2];
                w_candidate[1] = jets_chosen[0] + jets_chosen[2];
                w_candidate[2] = jets_chosen[0] + jets_chosen[1];

                TLorentzVector top_candidate[3];
                TLorentzVector fcnc_top_candidate[3];
                //}}}
                //--- calculation{{{
                std::vector<double> chi2;
                double w_mass[3], t_mass[3], fcnc_top_mass[3];
                for(int x = 0; x<3; ++x){
                    w_mass[x] =  w_candidate[x].M();

                    top_candidate[x] = w_candidate[x] + bjet;
                    t_mass[x] =  top_candidate[x].M();

                    fcnc_top_candidate[x] = tqh_q_chosen[x] + diphoton;
                    fcnc_top_mass[x] = fcnc_top_candidate[x].M();

                    chi2.push_back( Chi2_calculator_improved(w_mass[x], t_mass[x], fcnc_top_mass[x]) );
                }
                //}}}
                //--- sorting{{{
                int smallest_chi2_index = std::min_element(chi2.begin(),chi2.end()) - chi2.begin();
                double smallest_chi2 = *std::min_element(chi2.begin(),chi2.end());

                if(smallest_chi2 < chi2_min){
                    if(smallest_chi2_index == 0){
                        tqh_qjet = tqh_q_chosen[0];
                        wjets[0] = jets_chosen[1];
                        wjets[1] = jets_chosen[2];
                        index_tqh_qjet = i;
                        index_wjets[0] = j;
                        index_wjets[1] = k;
                    } else if(smallest_chi2_index == 1){
                        tqh_qjet = tqh_q_chosen[1];
                        wjets[0] = jets_chosen[0];
                        wjets[1] = jets_chosen[2];
                        index_tqh_qjet = j;
                        index_wjets[0] = i;
                        index_wjets[1] = k;
                    } else{
                        tqh_qjet = tqh_q_chosen[2];
                        wjets[0] = jets_chosen[0];
                        wjets[1] = jets_chosen[1];
                        index_tqh_qjet = k;
                        index_wjets[0] = i;
                        index_wjets[1] = j;
                    }

                    chi2_min = smallest_chi2;
                }
                //}}}
            }
        }
    }//end of looping jets
    w_boson  = wjets[0] + wjets[1];
    sm_top   = w_boson + bjet;
    fcnc_top = tqh_qjet + diphoton;

    std::vector<int> indices = {index_bjet, index_wjets[0], index_wjets[1], index_tqh_qjet};
    return indices;
}


std::vector<int> get_bjjq_indices_min_chi2_3x3(std::vector<TLorentzVector> Jets, std::vector<int> indices_bjet, TLorentzVector diphoton)
{
    vector<double> vec_chi2;
    vector<vector<int>> vec_indices_jet;
    for(std::size_t i=0; i!=indices_bjet.size(); ++i){
        double chi2_min = 99999;
        std::vector<int> index_jet_chi2 = get_bjjq_indices_chi2_3x3(Jets, indices_bjet[i], diphoton, chi2_min); // indices of b j j q
        vec_chi2.push_back(chi2_min);
        vec_indices_jet.push_back(index_jet_chi2);
        //printf("[check-imp] [%d] (b, j, j, q) = ", i);
        //printf("(%d, ", index_jet_chi2[0]);
        //printf("%d, " , index_jet_chi2[1]);
        //printf("%d, " , index_jet_chi2[2]);
        //printf("%d), ", index_jet_chi2[3]);
        //printf("chi2 = %7.3f\n", chi2_min);
    }

    int min_index =  std::min_element(vec_chi2.begin(), vec_chi2.end()) - vec_chi2.begin();
    double min    = *std::min_element(vec_chi2.begin(), vec_chi2.end());
    //printf("[curiosity] begin, end, diff = %d, %d, %d\n", vec_chi2.begin(), vec_chi2.end(), vec_chi2.end() - vec_chi2.begin());
    //printf("[curiosity] *begin, *end     = %f, %f\n", *vec_chi2.begin(), *(vec_chi2.end()-1));
    //printf("[check-imp] min: [%d] (b, j, j, q) = ", min_index);
    //printf("(%d, ", vec_indices_jet[min_index][0]);
    //printf("%d, " , vec_indices_jet[min_index][1]);
    //printf("%d, " , vec_indices_jet[min_index][2]);
    //printf("%d), ", vec_indices_jet[min_index][3]);
    //printf("chi2 = %7.3f\n", min);
    return vec_indices_jet[min_index];
}


//}}}
// chi-2 with 2x2 covariance matrix{{{
int get_q_index_min_chi2(std::vector<TLorentzVector> Jets, std::vector<int> indices_bjj, TLorentzVector diphoton)
{
    std::vector<int> indices;
    std::vector<double> top_fcnh_chi2;
    for(std::size_t i=0; i!=Jets.size(); ++i){
        if(i==indices_bjj[0] || i==indices_bjj[1] || i==indices_bjj[2]) continue; //skip the selected jets for bjj
        TLorentzVector top_fcnh_tmp = diphoton + Jets[i];
        double chi2 = (top_fcnh_tmp.M() - top_quark_mass) * (top_fcnh_tmp.M() - top_quark_mass);
        indices.push_back(i);
        top_fcnh_chi2.push_back(chi2);
        //printf("[check-ywk] q = ");
        //printf("%d, " , i);
        //printf("chi2 = %7.3f\n", chi2);
    }

    int min_index =  std::min_element(top_fcnh_chi2.begin(), top_fcnh_chi2.end()) - top_fcnh_chi2.begin();
    double min    = *std::min_element(top_fcnh_chi2.begin(), top_fcnh_chi2.end());
    //printf("[check-ywk] min: q = ");
    //printf("%d, " , indices[min_index]);
    //printf("chi2 = %7.3f\n", top_fcnh_chi2[min_index]);

    return indices[min_index];
}


std::vector<int> get_indices_chi2(std::vector<TLorentzVector> Jets, int index_bjet, double &chi2_min, bool is_chi2_modified)
{
    std::vector<int> indices_selected_jet(3, -999);
    indices_selected_jet[0] = index_bjet;
    TLorentzVector bjet = Jets[index_bjet];

    std::size_t num_jets = Jets.size();
    for(std::size_t i=0; i<num_jets; ++i){
        if(i==index_bjet) continue;//bypass bjet
        for(std::size_t j=i+1; j<num_jets; ++j){
            if(j==index_bjet) continue;//bypass bjet
            TLorentzVector w_candidate = Jets[i] + Jets[j];
            double w_mass = w_candidate.M();
            TLorentzVector top_candidate = w_candidate + bjet;
            double t_mass = top_candidate.M();
            double chi2 = !is_chi2_modified ? Chi2_calculator_simple(w_mass, t_mass) : Chi2_calculator_modified(w_mass, t_mass);
            if(chi2 < chi2_min){
                indices_selected_jet[1] = i;
                indices_selected_jet[2] = j;
                chi2_min = chi2;
            }
        }
    }//end of looping jets

    return indices_selected_jet; //indices of b j j
}


vector<int> get_bjj_indices_min_chi2(std::vector<TLorentzVector> Jets, std::vector<int> indices_bjet, bool is_chi2_modified)
{
    vector<double> vec_chi2;
    vector<vector<int>> vec_indices_jet;
    for(std::size_t i=0; i!=indices_bjet.size(); ++i){
        double chi2_min = 99999;
        std::vector<int> index_jet_chi2 = get_indices_chi2(Jets, indices_bjet[i], chi2_min, is_chi2_modified); // indices of b j j
        vec_chi2.push_back(chi2_min);
        vec_indices_jet.push_back(index_jet_chi2);
    }

    //printf("[check-ywk] ScanChain_ttHHadronic.h::get_bjj_indices_min_chi2::indices_bjet.size() = %d\n", indices_bjet.size());
    //printf("[check-ywk] ScanChain_ttHHadronic.h::get_bjj_indices_min_chi2::vec_chi2.size()     = %d\n", vec_chi2.size());

    int min_index =  std::min_element(vec_chi2.begin(), vec_chi2.end()) - vec_chi2.begin();
    double min    = *std::min_element(vec_chi2.begin(), vec_chi2.end());
    //printf("[curiosity] begin, end, diff = %d, %d, %d\n", vec_chi2.begin(), vec_chi2.end(), vec_chi2.end() - vec_chi2.begin());
    //printf("[curiosity] *begin, *end     = %f, %f\n", *vec_chi2.begin(), *(vec_chi2.end()-1));
    //printf("[check-ywk] min: [%d] (b, j, j) = ", min_index);
    //printf("(%d, ", vec_indices_jet[min_index][0]);
    //printf("%d, " , vec_indices_jet[min_index][1]);
    //printf("%d), ", vec_indices_jet[min_index][2]);
    //printf("chi2 = %7.3f\n\n", min);
    return vec_indices_jet[min_index];
}


vector<int> get_bjjq_indices_min_chi2(std::vector<TLorentzVector> Jets, std::vector<int> indices_bjet, TLorentzVector diphoton, bool is_chi2_modified)
{
    std::vector<int> indices = get_bjj_indices_min_chi2(Jets, indices_bjet, is_chi2_modified);
    int index_q = get_q_index_min_chi2(Jets, indices, diphoton);
    indices.push_back(index_q);
    return indices;
}
//}}}
//--------------------------------------------------//
#endif
