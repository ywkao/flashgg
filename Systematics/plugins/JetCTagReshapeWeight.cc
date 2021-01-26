#include <TSystem.h>
#include "flashgg/Systematics/interface/ObjectSystMethodBinnedByFunctor.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "flashgg/Systematics/interface/ctag_reshaping.h"


namespace flashgg {

    class JetCTagReshapeWeight: public ObjectSystMethodBinnedByFunctor<flashgg::Jet, int>
    {

    public:
        typedef StringCutObjectSelector<Jet, true> selector_type;

        JetCTagReshapeWeight( const edm::ParameterSet &conf, edm::ConsumesCollector && iC, const GlobalVariablesComputer * gv );
        float makeWeight( const flashgg::Jet &y, int syst_shift ) override;
        std::string shiftLabel( int syst_shift ) const override;

    private:
        selector_type overall_range_;
        bool debug_;
        std::string cTag_;
        bool isApply_;
        std::string cTagReshapeSystOption_;
        std::string cTagReshapeFile_;
        retrieve_scale_factor sf_retriever;
    };

    JetCTagReshapeWeight::JetCTagReshapeWeight( const edm::ParameterSet &conf, edm::ConsumesCollector && iC, const GlobalVariablesComputer * gv ) : 
        ObjectSystMethodBinnedByFunctor( conf, std::forward<edm::ConsumesCollector>(iC), gv  ),
        overall_range_( conf.getParameter<std::string>( "OverallRange" ) ),        
        debug_( conf.getUntrackedParameter<bool>( "Debug", false ) ),
        cTag_( conf.getParameter<std::string>("cTag") ),
        isApply_( conf.getParameter<bool>("isApply") )
    {
        this->setMakesWeight( true );
        cTagReshapeSystOption_ = conf.getParameter<std::string>( "cTagReshapeSystOption"); 
        cTagReshapeFile_ = conf.getParameter<edm::FileInPath>("cTagCalibrationFile").fullPath();
        sf_retriever.set_ctag_reshape_file(cTagReshapeFile_);
    }

    std::string JetCTagReshapeWeight::shiftLabel( int syst_value ) const
    {
        std::string result;
        if( syst_value == 0 ) {
            result = Form( "%sCentral", label().c_str() );
        } else if( syst_value > 0 ) {
            result = Form( "%sUp%.2dsigma", label().c_str(), syst_value );
        } else {
            result = Form( "%sDown%.2dsigma", label().c_str(), -1 * syst_value );
        }
        return result;
    }

    float JetCTagReshapeWeight::makeWeight( const flashgg::Jet &obj, int syst_shift ) 
    {
        if( this->debug_ ) {
            sf_retriever.debug_mode();
            std::cout<<"In JetCTagReshapeProducer and syst_shift="<<  syst_shift <<std::endl;
        }

        float theWeight = 1.;

        if(!isApply_) {
            return theWeight;
        } else {
            if( overall_range_( obj ) ) {

                float central = 1., errup = 1., errdown = 1.;

                //obtaining scale factors
                float JetPt = obj.pt();
                float JetEta = fabs(obj.eta());
                int JetFlav = obj.hadronFlavour();

                //float cvsl = obj.bDiscriminator("pfDeepCSVJetTags:probc") / (obj.bDiscriminator("pfDeepCSVJetTags:probc") + obj.bDiscriminator("pfDeepCSVJetTags:probudsg")) ;
                //float cvsb = obj.bDiscriminator("pfDeepCSVJetTags:probc") / (obj.bDiscriminator("pfDeepCSVJetTags:probc") + obj.bDiscriminator("pfDeepCSVJetTags:probb") + obj.bDiscriminator("pfDeepCSVJetTags:probbb")) ;


                if (this->debug_) {
                    std::cout << "In JetCTagReshapeProducer, jet discriminator values: " << std::endl;
                    std::cout << "mini_pfDeepFlavourJetTags:probc: " << obj.bDiscriminator("mini_pfDeepFlavourJetTags:probc") << std::endl;
                    std::cout << "mini_pfDeepFlavourJetTags:probuds: " << obj.bDiscriminator("mini_pfDeepFlavourJetTags:probuds") << std::endl;
                    std::cout << "mini_pfDeepFlavourJetTags:probb: " << obj.bDiscriminator("mini_pfDeepFlavourJetTags:probb") << std::endl;
                    std::cout << "mini_pfDeepFlavourJetTags:probbb: " << obj.bDiscriminator("mini_pfDeepFlavourJetTags:probbb") << std::endl;
                    std::cout << "mini_pfDeepFlavourJetTags:problepb: " << obj.bDiscriminator("mini_pfDeepFlavourJetTags:problepb") << std::endl;
                    std::cout << "mini_pfDeepFlavourJetTags:probg: " << obj.bDiscriminator("mini_pfDeepFlavourJetTags:probg") << std::endl;

                }

                float cvsl;
                float cvsb;

                if (obj.bDiscriminator("mini_pfDeepFlavourJetTags:probc") == -1000 || obj.bDiscriminator("mini_pfDeepFlavourJetTags:probuds") == -1000 || obj.bDiscriminator("mini_pfDeepFlavourJetTags:probg") == -10000 || obj.bDiscriminator("mini_pfDeepFlavourJetTags:probb") == -1000 || obj.bDiscriminator("mini_pfDeepFlavourJetTags:probbb") == -1000 || obj.bDiscriminator("mini_pfDeepFlavourJetTags:problepb") == -1000) {
                    cvsl = 0.0;
                    cvsb = 0.0;                    
                }

                else {
                    cvsl = obj.bDiscriminator("mini_pfDeepFlavourJetTags:probc") / (obj.bDiscriminator("mini_pfDeepFlavourJetTags:probc") + obj.bDiscriminator("mini_pfDeepFlavourJetTags:probuds") + obj.bDiscriminator("mini_pfDeepFlavourJetTags:probg")) ;
                    cvsb = obj.bDiscriminator("mini_pfDeepFlavourJetTags:probc") / (obj.bDiscriminator("mini_pfDeepFlavourJetTags:probc") + obj.bDiscriminator("mini_pfDeepFlavourJetTags:probb") + obj.bDiscriminator("mini_pfDeepFlavourJetTags:probbb") + obj.bDiscriminator("mini_pfDeepFlavourJetTags:problepb")) ;
                }

                if( this->debug_ ) {
                    std::cout << " In JetCTagReshapeWeight before calib reader: " << shiftLabel( syst_shift ) << ": Object has pt= " << obj.pt() << " eta=" << obj.eta() << " flavour=" << obj.hadronFlavour()
                              << " values for scale factors : "<< JetFlav << ", c-tagger = " << cTag_
                              << " CvsL values: " << cvsl
                              << " CvsB values: " << cvsb <<endl;
                }

                //get scale factors from calib reader
                double jet_scalefactor = 1.0;
                double jet_scalefactor_up =  1.0;
                double jet_scalefactor_do =  1.0;
                
                if( cvsl < 0.0 ) cvsl = -0.05;
                if( cvsl > 1.0 ) cvsl = 1.0;
                if( cvsb < 0.0 ) cvsb = -0.05;
                if( cvsb > 1.0 ) cvsb = 1.0;

                TString type_flavour;
                if(JetFlav == 5){ // b jets
                    type_flavour = "b";
                } else if(JetFlav == 4){ // c jets
                    type_flavour = "c";
                } else{ //light jets
                    type_flavour = "l";
                }

                TString string_scale_factor     = "SF" + type_flavour + "_" + "hist";
                TString string_scale_factor_sys = "SF" + type_flavour + "_" + "hist" + "_" + cTagReshapeSystOption_.c_str();

                jet_scalefactor = sf_retriever.get_scale_factor(string_scale_factor, cvsl, cvsb);
                jet_scalefactor_up = sf_retriever.get_scale_factor(string_scale_factor_sys + "Up", cvsl, cvsb); 
                jet_scalefactor_do = sf_retriever.get_scale_factor(string_scale_factor_sys + "Down", cvsl, cvsb); 
                if( this->debug_ )  { std::cout << " In JetCTagReshapeWeight Systematics : "<< cTagReshapeSystOption_ << " : "
                                                << jet_scalefactor << " " << jet_scalefactor_up << " " << jet_scalefactor_do <<std::endl; }
                
                
                if( this->debug_ ) {
                    std::cout << " In JetCTagReshapeWeight after obtaining SF: " << shiftLabel( syst_shift ) << " : Object has pt= " << obj.pt() << " eta=" << obj.eta() << " flavour=" << obj.hadronFlavour()
                              << " scale factors : " << jet_scalefactor << " " << jet_scalefactor_up << " " << jet_scalefactor_do << std::endl;
                }
                
                if ( syst_shift < 0 ){
                    jet_scalefactor_do = abs(syst_shift) * (jet_scalefactor_do - jet_scalefactor) + jet_scalefactor;
                }
                if ( syst_shift > 0 ){
                    jet_scalefactor_up = abs(syst_shift) * (jet_scalefactor_up - jet_scalefactor) + jet_scalefactor;
                }
                
                if( this->debug_ ) {
                    std::cout << " In JetCTagReshapeWeight : " << shiftLabel( syst_shift ) << " : Object has pt= " << obj.pt() << " eta=" << obj.eta() << " flavour=" << obj.hadronFlavour()
                              << " scale factors : "<< jet_scalefactor <<" "<< jet_scalefactor_up <<" "<< jet_scalefactor_do << std::endl;
                }
                
                central = jet_scalefactor ;
                errdown = jet_scalefactor_do ;
                errup   = jet_scalefactor_up ;
                            
                theWeight = central;
                if ( syst_shift < 0 ) theWeight = errdown;
                if ( syst_shift > 0 ) theWeight = errup;
                
                if( this->debug_ ) {
                    std::cout << " In JetCTagReshapeWeight : " << shiftLabel( syst_shift ) <<  " : Object has pt= " << obj.pt() << " eta=" << obj.eta() << " flavour=" << obj.hadronFlavour()
                              << " and we apply a weight of " << theWeight << std::endl;
                }
            }
            return theWeight;
        }
    } // end of float JetCTagReshapeWeight::makeWeight()

}

DEFINE_EDM_PLUGIN( FlashggSystematicJetMethodsFactory,
                   flashgg::JetCTagReshapeWeight,
                   "FlashggJetCTagReshapeWeight" );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
