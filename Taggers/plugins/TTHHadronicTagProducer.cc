#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/TTHHadronicTag.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include <vector>
#include <algorithm>
#include <string>
#include <utility>

#include "TMVA/Reader.h"
#include "TRandom.h"

TRandom* myRandHadronic = new TRandom(0);

using namespace std;
using namespace edm;

namespace flashgg {

    class TTHHadronicTagProducer : public EDProducer
    {

    public:
        typedef math::XYZPoint Point;

        TTHHadronicTagProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;

        std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > tokenJets_;
        EDGetTokenT<View<DiPhotonCandidate> > diPhotonToken_;
        std::vector<edm::InputTag> inputTagJets_;
        EDGetTokenT<View<Electron> > electronToken_;
        EDGetTokenT<View<flashgg::Muon> > muonToken_;
        EDGetTokenT<View<reco::Vertex> > vertexToken_;
        EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;
        EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
        //EDGetTokenT<View<reco::GenParticleCollection> > genParticleCollectionToken_;
        EDGetTokenT<View<flashgg::Met> > METToken_;
        EDGetTokenT<int> stage0catToken_, stage1catToken_, njetsToken_;
        EDGetTokenT<float> pTHToken_,pTVToken_;
        EDGetTokenT<double> rhoTag_;
        string systLabel_;


        typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;
        bool useTTHHadronicMVA_;
        double tthHadMVAThreshold_;
        //---thresholds---
        //---photons
        double MVAThreshold_;
        double MVATTHHMVAThreshold_;
        double PhoMVAThreshold_;
        double leadPhoPtThreshold_;
        bool   leadPhoUseVariableTh_;
        double leadPhoOverMassThreshold_;
        double leadPhoOverMassTTHHMVAThreshold_;
        double subleadPhoPtThreshold_;
        bool   subleadPhoUseVariableTh_;
        double subleadPhoOverMassThreshold_;

        double diphotonMassThreshold_;

        //---jets
        double jetPtThreshold_;
        double jetEtaThreshold_;
        double dRJetPhoLeadCut_;
        double dRJetPhoSubleadCut_;
        vector<double> bDiscriminator_;
        double jetsNumberThreshold_;
        double bjetsNumberThreshold_;
        double bjetsLooseNumberThreshold_;
        double jetsNumberTTHHMVAThreshold_;
        double bjetsNumberTTHHMVAThreshold_;
        double bjetsLooseNumberTTHHMVAThreshold_;
        string bTag_;
        //leptons
        double leptonPtThreshold_;
        double muonEtaThreshold_;
        vector<double> nonTrigMVAThresholds_;
        vector<double> nonTrigMVAEtaCuts_;
        double electronIsoThreshold_;
        double elMiniIsoEBThreshold_;
        double elMiniIsoEEThreshold_;
        double electronNumOfHitsThreshold_;
        vector<double>  electronEtaThresholds_;
        double muPFIsoSumRelThreshold_;
        double muMiniIsoSumRelThreshold_;
        double TransverseImpactParam_EB;
        double LongitudinalImpactParam_EB;
        double TransverseImpactParam_EE;
        double LongitudinalImpactParam_EE;

        bool useStdLeptonID_;
        bool useElectronMVARecipe_;
        bool useElectronLooseID_;

        unique_ptr<TMVA::Reader>TThMva_;
        FileInPath tthMVAweightfile_;
        string _MVAMethod;

        int jetcount_;
        float nJets_;
        int njets_btagloose_;
        int njets_btagmedium_;
        int njets_btagtight_;
        double idmva1_;
        double idmva2_;
        float leadJetPt_;
        float subLeadJetPt_;
        float sumJetPt_;
        float maxBTagVal_;
        float secondMaxBTagVal_;
        float tthMvaVal_;

    };

    TTHHadronicTagProducer::TTHHadronicTagProducer( const ParameterSet &iConfig ) :
        diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
        inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
        electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
        muonToken_( consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
        vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
        mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag>( "MVAResultTag" ) ) ),
        genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
        //genParticleCollectionToken_( consumes<View<reco::GenParticleCollection> >( iConfig.getParameter<InputTag> ( "GenParticleCollectionTag" ) ) ),
        METToken_( consumes<View<flashgg::Met> >( iConfig.getParameter<InputTag>( "MetTag" ) ) ),
        rhoTag_( consumes<double>( iConfig.getParameter<InputTag>( "rhoTag" ) ) ),
        systLabel_( iConfig.getParameter<string> ( "SystLabel" ) ),
        _MVAMethod( iConfig.getParameter<string> ( "MVAMethod" ) )
    {
        ParameterSet HTXSps = iConfig.getParameterSet( "HTXSTags" );
        stage0catToken_ = consumes<int>( HTXSps.getParameter<InputTag>("stage0cat") );
        stage1catToken_ = consumes<int>( HTXSps.getParameter<InputTag>("stage1cat") );
        njetsToken_ = consumes<int>( HTXSps.getParameter<InputTag>("njets") );
        pTHToken_ = consumes<float>( HTXSps.getParameter<InputTag>("pTH") );
        pTVToken_ = consumes<float>( HTXSps.getParameter<InputTag>("pTV") );

        MVAThreshold_ = iConfig.getParameter<double>( "MVAThreshold");
        MVATTHHMVAThreshold_ = iConfig.getParameter<double>( "MVATTHHMVAThreshold");
        PhoMVAThreshold_ = iConfig.getParameter<double>( "PhoMVAThreshold");

        leptonPtThreshold_ = iConfig.getParameter<double>( "leptonPtThreshold");
        muonEtaThreshold_ = iConfig.getParameter<double>( "muonEtaThreshold");
        leadPhoPtThreshold_ = iConfig.getParameter<double>( "leadPhoPtThreshold");
        leadPhoUseVariableTh_ = iConfig.getParameter<bool>( "leadPhoUseVariableThreshold");
        leadPhoOverMassThreshold_ = iConfig.getParameter<double>( "leadPhoOverMassThreshold");
        leadPhoOverMassTTHHMVAThreshold_ = iConfig.getParameter<double>( "leadPhoOverMassTTHHMVAThreshold");
        subleadPhoPtThreshold_ = iConfig.getParameter<double>( "subleadPhoPtThreshold");
        subleadPhoUseVariableTh_ = iConfig.getParameter<bool>( "subleadPhoUseVariableThreshold");
        subleadPhoOverMassThreshold_ = iConfig.getParameter<double>( "subleadPhoOverMassThreshold");
        jetPtThreshold_ = iConfig.getParameter<double>( "jetPtThreshold");
        jetEtaThreshold_ = iConfig.getParameter<double>( "jetEtaThreshold");
        dRJetPhoLeadCut_ = iConfig.getParameter<double>( "dRJetPhoLeadCut");
        dRJetPhoSubleadCut_ = iConfig.getParameter<double>( "dRJetPhoSubleadCut");
        bDiscriminator_ = iConfig.getParameter<vector<double > >( "bDiscriminator");
        jetsNumberThreshold_ = iConfig.getParameter<int>( "jetsNumberThreshold");
        bjetsNumberThreshold_ = iConfig.getParameter<int>( "bjetsNumberThreshold");
        bTag_ = iConfig.getParameter<string> ( "bTag");

        useTTHHadronicMVA_ = iConfig.getParameter<bool>( "useTTHHadronicMVA");
        tthHadMVAThreshold_ =  iConfig.getParameter<double>( "tthHadMVAThreshold");
        bjetsLooseNumberThreshold_ = iConfig.getParameter<int>( "bjetsLooseNumberThreshold");
        jetsNumberTTHHMVAThreshold_ = iConfig.getParameter<int>( "jetsNumberTTHHMVAThreshold");
        bjetsNumberTTHHMVAThreshold_ = iConfig.getParameter<int>( "bjetsNumberTTHHMVAThreshold");
        bjetsLooseNumberTTHHMVAThreshold_ = iConfig.getParameter<int>( "bjetsLooseNumberTTHHMVAThreshold");

        muPFIsoSumRelThreshold_ = iConfig.getParameter<double>( "muPFIsoSumRelThreshold");
        muMiniIsoSumRelThreshold_ = iConfig.getParameter<double>( "muMiniIsoSumRelThreshold");
        nonTrigMVAThresholds_ =  iConfig.getParameter<vector<double > >( "nonTrigMVAThresholds");
        nonTrigMVAEtaCuts_ =  iConfig.getParameter<vector<double > >( "nonTrigMVAEtaCuts");
        electronIsoThreshold_ = iConfig.getParameter<double>( "electronIsoThreshold");
        elMiniIsoEBThreshold_ = iConfig.getParameter<double>( "elMiniIsoEBThreshold");
        elMiniIsoEEThreshold_ = iConfig.getParameter<double>( "elMiniIsoEEThreshold");
        electronNumOfHitsThreshold_ = iConfig.getParameter<double>( "electronNumOfHitsThreshold");
        TransverseImpactParam_EB = iConfig.getParameter<double>( "TransverseImpactParamEB");
        LongitudinalImpactParam_EB = iConfig.getParameter<double>( "LongitudinalImpactParamEB");
        TransverseImpactParam_EE = iConfig.getParameter<double>( "TransverseImpactParamEE");
        LongitudinalImpactParam_EE = iConfig.getParameter<double>( "LongitudinalImpactParamEE");
        electronEtaThresholds_ = iConfig.getParameter<vector<double > >( "electronEtaThresholds");
        useStdLeptonID_=iConfig.getParameter<bool>("useStdLeptonID");
        useElectronMVARecipe_=iConfig.getParameter<bool>("useElectronMVARecipe");
        useElectronLooseID_=iConfig.getParameter<bool>("useElectronLooseID");
        

        diphotonMassThreshold_ = iConfig.getParameter<double>( "diphotonMassThreshold");

        tthMVAweightfile_ = iConfig.getParameter<edm::FileInPath>( "tthMVAweightfile" ); 

        nJets_ = 0;
        leadJetPt_ = 0.;
        subLeadJetPt_ = 0.;
        sumJetPt_ = 0.;
        maxBTagVal_ = -999.;
        secondMaxBTagVal_ = -999.;

        if (_MVAMethod != ""){
            TThMva_.reset( new TMVA::Reader( "!Color:Silent" ) );
            TThMva_->AddVariable( "nJets", &nJets_);
            TThMva_->AddVariable( "maxBTagVal",&maxBTagVal_);
            TThMva_->AddVariable( "secondMaxBTagVal", &secondMaxBTagVal_);
            TThMva_->AddVariable( "leadJetPt", &leadJetPt_);
        
            TThMva_->BookMVA( _MVAMethod.c_str() , tthMVAweightfile_.fullPath() );
        
        }       

        for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
            auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
            tokenJets_.push_back(token);
        }

        produces<vector<TTHHadronicTag> >();
        produces<vector<TagTruthBase> >();
    }

    void TTHHadronicTagProducer::produce( Event &evt, const EventSetup & )
    {

        Handle<int> stage0cat, stage1cat, njets;
        Handle<float> pTH, pTV;
        evt.getByToken(stage0catToken_, stage0cat);
        evt.getByToken(stage1catToken_,stage1cat);
        evt.getByToken(njetsToken_,njets);
        evt.getByToken(pTHToken_,pTH);
        evt.getByToken(pTVToken_,pTV);

        //Handle<View<flashgg::Jet> > theJets;
        //evt.getByToken( thejetToken_, theJets );
        // const PtrVector<flashgg::Jet>& jetPointers = theJets->ptrVector();
        JetCollectionVector Jets( inputTagJets_.size() );
        for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
            evt.getByToken( tokenJets_[j], Jets[j] );
        }

        Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
        evt.getByToken( diPhotonToken_, diPhotons );
        // const PtrVector<flashgg::DiPhotonCandidate>& diPhotonPointers = diPhotons->ptrVector();

        Handle<View<flashgg::Muon> > theMuons;
        evt.getByToken( muonToken_, theMuons );

        Handle<View<flashgg::Electron> > theElectrons;
        evt.getByToken( electronToken_, theElectrons );

        Handle<View<flashgg::Met> > theMet_;
        evt.getByToken( METToken_, theMet_ );

        edm::Handle<double>  rho;
        evt.getByToken(rhoTag_,rho);
        double rho_    = *rho;

        Handle<View<reco::Vertex> > vertices;
        evt.getByToken( vertexToken_, vertices );

        Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;
        evt.getByToken( mvaResultToken_, mvaResults );
        // const PtrVector<flashgg::DiPhotonMVAResult>& mvaResultPointers = mvaResults->ptrVector();

        Handle<View<reco::GenParticle> > genParticles;

        //Handle<View<reco::GenParticle> > genParticleCollection;
        //edm::Handle<reco::GenParticleCollection> genpsHandle;
        //evt.getByToken(genParticlesToken, genpsHandle);
        //const vector<reco::GenParticle>* genps_coll = genParticleCollection.product();

        std::unique_ptr<vector<TTHHadronicTag> > tthhtags( new vector<TTHHadronicTag> );
        std::unique_ptr<vector<TagTruthBase> > truths( new vector<TagTruthBase> );

        Point higgsVtx;
        if( ! evt.isRealData() ) {
            evt.getByToken( genParticleToken_, genParticles );
            for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                if( pdgid == 25 || pdgid == 22 ) {
                    higgsVtx = genParticles->ptrAt( genLoop )->vertex();
                    break;
                }
            }
        }

        edm::RefProd<vector<TagTruthBase> > rTagTruth = evt.getRefBeforePut<vector<TagTruthBase> >();
        unsigned int idx = 0;

        std::vector<edm::Ptr<flashgg::Muon> > goodMuons;
        if( !useStdLeptonID_) {
            goodMuons = selectAllMuonsSum16( theMuons->ptrs(), vertices->ptrs(), muonEtaThreshold_ , leptonPtThreshold_, muMiniIsoSumRelThreshold_ );
        } else {
            goodMuons = selectAllMuons( theMuons->ptrs(), vertices->ptrs(), muonEtaThreshold_ , leptonPtThreshold_, muPFIsoSumRelThreshold_ );
        }
        
       
        std::vector<edm::Ptr<Electron> > goodElectrons ;

        //if( !useStdLeptonID_) {
             // goodElectrons = selectAllElectrons( theElectrons->ptrs(), vertices->ptrs(), leptonPtThreshold_, 
             //                                     TransverseImpactParam_, LongitudinalImpactParam_, nonTrigMVAThresholds_, nonTrigMVAEtaCuts_,
             //                                     electronIsoThreshold_, electronNumOfHitsThreshold_, electronEtaThresholds_ );
            //goodElectrons = selectAllElectronsSum16( theElectrons->ptrs(), vertices->ptrs(), 
            //                                       leptonPtThreshold_, electronEtaThresholds_,
            //                                       true, true, elMiniIsoEBThreshold_, elMiniIsoEEThreshold_);
            //goodElectrons = selectAllElectronsSum16( theElectrons->ptrs(), vertices->ptrs(), leptonPtThreshold_, electronEtaThresholds_,
            //                                         true, true, elMiniIsoEBThreshold_, elMiniIsoEEThreshold_,
            //                                         TransverseImpactParam_EB, LongitudinalImpactParam_EB, TransverseImpactParam_EE, LongitudinalImpactParam_EE,
            //                                         rho_, evt.isRealData() );
        //} else {
            goodElectrons = selectStdAllElectrons(theElectrons->ptrs(), vertices->ptrs(), leptonPtThreshold_, electronEtaThresholds_,
                                                  useElectronMVARecipe_, useElectronLooseID_,
                                                  rho_, evt.isRealData() );
            //}
        
        for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ ) {

            if( goodElectrons.size() > 0 ||  goodMuons.size() > 0 )  continue; 

            jetcount_ = 0;
            nJets_ = 0;
            njets_btagloose_ = 0;
            njets_btagmedium_ = 0;
            njets_btagtight_ = 0;
            idmva1_ = -999.;
            idmva2_ = -999.;
            leadJetPt_ = 0.;
            subLeadJetPt_ = 0.;
            sumJetPt_ = 0.;
            maxBTagVal_ = -999.;
            secondMaxBTagVal_ = -999.;
            tthMvaVal_ = -999.;

            unsigned int jetCollectionIndex = diPhotons->ptrAt( diphoIndex )->jetCollectionIndex();

            std::vector<edm::Ptr<flashgg::Jet> > JetVect;
            JetVect.clear();
            std::vector<edm::Ptr<flashgg::Jet> > BJetVect;
            BJetVect.clear();
            std::vector<edm::Ptr<flashgg::Jet> > BJetTTHHMVAVect;
            BJetTTHHMVAVect.clear();
            std::vector<float> JetBTagVal;
            JetBTagVal.clear();

            edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );

            if (dipho->mass() < diphotonMassThreshold_) { continue; }

            idmva1_ = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
            idmva2_ = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );

            if( idmva1_ <= PhoMVAThreshold_ || idmva2_ <= PhoMVAThreshold_ ) { continue; }

            edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );

            if( !dipho->leadingPhoton()->passElectronVeto() || !dipho->subLeadingPhoton()->passElectronVeto() ) { continue; }

            double leadPhoPtCut = leadPhoPtThreshold_;
            double subleadPhoPtCut = subleadPhoPtThreshold_;
            if( leadPhoUseVariableTh_ )
            { 
                leadPhoPtCut = leadPhoOverMassThreshold_ * dipho->mass(); 
                if(useTTHHadronicMVA_){
                    leadPhoPtCut = leadPhoOverMassTTHHMVAThreshold_ * dipho->mass();
                }
            }
            if( subleadPhoUseVariableTh_ )
            { subleadPhoPtCut = subleadPhoOverMassThreshold_ * dipho->mass(); }
            double diphoMVAcut = MVAThreshold_;
            if(useTTHHadronicMVA_){
                    diphoMVAcut = MVATTHHMVAThreshold_;
            }

            if( dipho->leadingPhoton()->pt() < leadPhoPtCut || dipho->subLeadingPhoton()->pt() < subleadPhoPtCut ) { continue; }
            if( mvares->mvaValue() < diphoMVAcut ) { continue; }

            for( unsigned int jetIndex = 0; jetIndex < Jets[jetCollectionIndex]->size() ; jetIndex++ ) {
                edm::Ptr<flashgg::Jet> thejet = Jets[jetCollectionIndex]->ptrAt( jetIndex );
                if( fabs( thejet->eta() ) > jetEtaThreshold_ ) { continue; }
                if(!thejet->passesJetID  ( flashgg::Loose ) ) { continue; }
                float dRPhoLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ;
                float dRPhoSubLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->subLeadingPhoton()->superCluster()->eta(),
                                                dipho->subLeadingPhoton()->superCluster()->phi() );

                if( dRPhoLeadJet < dRJetPhoLeadCut_ || dRPhoSubLeadJet < dRJetPhoSubleadCut_ ) { continue; }
                if( thejet->pt() < jetPtThreshold_ ) { continue; }

                jetcount_++;
                nJets_ = jetcount_;
                JetVect.push_back( thejet );
                
                float jetPt = thejet->pt();
                if(jetPt > leadJetPt_){
                    if(leadJetPt_ > subLeadJetPt_) { subLeadJetPt_ = leadJetPt_; }
                    leadJetPt_ = jetPt;
                } else if(jetPt > subLeadJetPt_){
                    subLeadJetPt_ = jetPt;
                }
                sumJetPt_ += jetPt;
                
                float bDiscriminatorValue = -2.;
                bDiscriminatorValue = thejet->bDiscriminator( bTag_ );

                if(bDiscriminatorValue > maxBTagVal_){
                    BJetTTHHMVAVect.insert( BJetTTHHMVAVect.begin(), thejet );
                    if(BJetTTHHMVAVect.size() >= 3){ BJetTTHHMVAVect.pop_back(); }
                    if(maxBTagVal_ > secondMaxBTagVal_) { secondMaxBTagVal_ = maxBTagVal_; }
                    maxBTagVal_ = bDiscriminatorValue;

                } else if(bDiscriminatorValue > secondMaxBTagVal_){
                    secondMaxBTagVal_ = bDiscriminatorValue;
                    if(BJetTTHHMVAVect.size() >= 2){BJetTTHHMVAVect.pop_back();} 
                    BJetTTHHMVAVect.push_back( thejet );
                }
                
                JetBTagVal.push_back( bDiscriminatorValue );
                if( bDiscriminatorValue > bDiscriminator_[0] ) njets_btagloose_++;
                if( bDiscriminatorValue > bDiscriminator_[1] ){
                
    
                    njets_btagmedium_++;
                    //JetVect.pop_back();
                    BJetVect.push_back( thejet );
                }
                if( bDiscriminatorValue > bDiscriminator_[2] ) njets_btagtight_++;
            }

            if(useTTHHadronicMVA_){
                BJetVect.clear();
                BJetVect = BJetTTHHMVAVect;

                 if(njets_btagloose_ >= bjetsLooseNumberTTHHMVAThreshold_ && njets_btagmedium_ >= bjetsNumberTTHHMVAThreshold_ && jetcount_ >= jetsNumberTTHHMVAThreshold_ && _MVAMethod != ""){

                     tthMvaVal_ = TThMva_->EvaluateMVA( _MVAMethod.c_str() );
                     //                mvares_tth = TThMva_->EvaluateMVA( "BDT" );

                     /*
                    cout << "input variables : " << endl;
                    cout << "nJets_ = " << nJets_ << endl;
                    cout << "maxBTagVal_ = " << maxBTagVal_ << endl;
                    cout << "secondMaxBTagVal_ = " << secondMaxBTagVal_ << endl;
                    cout << "leadJetPt_ = " << leadJetPt_ << endl;

                    cout << "mva result :" << endl;
                    cout << "tthMvaVal_ = " << tthMvaVal_  << endl;
                     */
                 }
            }

            bool isTTHHadronicTagged = false;
            
            if( !useTTHHadronicMVA_ && njets_btagloose_ >= bjetsLooseNumberThreshold_ && njets_btagmedium_ >= bjetsNumberThreshold_ && jetcount_ >= jetsNumberThreshold_ ) {
                
                isTTHHadronicTagged = true;
            
            } else if ( useTTHHadronicMVA_ && njets_btagloose_ >= bjetsLooseNumberTTHHMVAThreshold_ && njets_btagmedium_ >= bjetsNumberTTHHMVAThreshold_ && jetcount_ >= jetsNumberTTHHMVAThreshold_ && tthMvaVal_ > tthHadMVAThreshold_ ) {
                
                isTTHHadronicTagged = true;
            }

            if( isTTHHadronicTagged ) {
                TTHHadronicTag tthhtags_obj( dipho, mvares, JetVect, BJetVect );
                tthhtags_obj.setNjet( jetcount_ );
                tthhtags_obj.setNBLoose( njets_btagloose_ );
                tthhtags_obj.setNBMedium( njets_btagmedium_ );
                tthhtags_obj.setNBTight( njets_btagtight_ );
                tthhtags_obj.setDiPhotonIndex( diphoIndex );
                tthhtags_obj.setLeadJetPt( leadJetPt_ );
                tthhtags_obj.setSubLeadJetPt( subLeadJetPt_ );
                tthhtags_obj.setSumJetPt( sumJetPt_ );
                tthhtags_obj.setMaxBTagVal( maxBTagVal_ );
                tthhtags_obj.setSecondMaxBTagVal( secondMaxBTagVal_ );
                tthhtags_obj.setSystLabel( systLabel_ );
                tthhtags_obj.setMVAres(tthMvaVal_);
                tthhtags_obj.setDiphoMVARes(mvares->mvaValue());
                cout << "random number: " << myRandHadronic->Rndm() << endl;   
                tthhtags_obj.setRand(myRandHadronic->Rndm());

                if( theMet_ -> size() != 1 )
                  std::cout << "WARNING number of MET is not equal to 1" << std::endl;
                Ptr<flashgg::Met> Met = theMet_->ptrAt( 0 );
                tthhtags_obj.setMetPt((float)Met->pt());
                tthhtags_obj.setMetPhi((float)Met->phi());

                // Gen photon info
                /*
                if( ! evt.isRealData() ) {
                    evt.getByToken( genParticleToken_, genParticles );
                    for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                        int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                        if (abs(pdgid) != 22)   continue;
                        double pt = genParticles->ptrAt( genLoop )->p4().pt();
                        bool isPromptFinalState = genParticles->ptrAt( genLoop )->isPromptFinalState();
                        if (!isPromptFinalState)    continue;
                        if (pt < 10)            continue;
                        double deltaR_pho1 = deltaR( genParticles->ptrAt( genLoop )->p4().eta(), genParticles->ptrAt( genLoop )->p4().phi(), dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ; 
                        double deltaR_pho2 = deltaR( genParticles->ptrAt( genLoop )->p4().eta(), genParticles->ptrAt( genLoop )->p4().phi(), dipho->subLeadingPhoton()->superCluster()->eta(), dipho->subLeadingPhoton()->superCluster()->phi() ) ;
                        if (deltaR_pho1 < gen_photon_deltaR_thresh) {
                          
                        }
                        else if (deltaR_pho2 < gen_photon_deltaR_thresh) {

                        }
 
                    }
                } */
                // Gen lepton info
                if( ! evt.isRealData() ) {
                    evt.getByToken( genParticleToken_, genParticles );
                    int nGoodEls(0), nGoodMus(0), nGoodElsFromTau(0), nGoodMusFromTau(0), nGoodTaus(0);
                    cout << "Number of gen particles: " << genParticles->size() << endl;
                    for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                        int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                        double pt = genParticles->ptrAt( genLoop )->p4().pt();
                        int status = genParticles->ptrAt( genLoop )->status();
                        bool isPromptFinalState = genParticles->ptrAt( genLoop )->isPromptFinalState();
                        bool isPromptDecayed = genParticles->ptrAt( genLoop )->isPromptDecayed();
                        bool isDirectPromptTauDecayProductFinalState = genParticles->ptrAt( genLoop )->isDirectPromptTauDecayProductFinalState();
                        if (pt < 20) continue;
                        if (abs(pdgid) == 11 || abs(pdgid) == 13 || abs(pdgid) == 15) {
                            cout << "Found a gen lepton/tau with pT > 20" << endl;
                            cout << "pdgid: " << pdgid << endl;
                            cout << "pt: " << pt << endl;
                            cout << "status: " << status << endl;
                            cout << "isPromptFinalState: " << isPromptFinalState << endl;
                            cout << "isPromptDecayed: " << isPromptDecayed << endl;
                            cout << "isDirectPromptTauDecayProductFinalState: " << isDirectPromptTauDecayProductFinalState << endl;
                        }
                        if (abs(pdgid) == 11 && status == 1 && (isPromptFinalState || isDirectPromptTauDecayProductFinalState)) {
                            nGoodEls++;
                            if (isDirectPromptTauDecayProductFinalState)
                                nGoodElsFromTau++;
                        } 
                        else if (abs(pdgid) == 13 && status == 1 && (isPromptFinalState || isDirectPromptTauDecayProductFinalState)) {
                            nGoodMus++;
                            if (isDirectPromptTauDecayProductFinalState)
                                nGoodMusFromTau++;
                        }
                        if (abs(pdgid) == 15 && status == 2 && isPromptDecayed) {
                            nGoodTaus++;
                        }
                    } // end gen loop
                    bool lead_photon_is_electron(false), sublead_photon_is_electron(false);
                    double lead_photon_eta = dipho->leadingPhoton()->eta();
                    double lead_photon_phi = dipho->leadingPhoton()->phi();
                    double sublead_photon_eta = dipho->subLeadingPhoton()->eta();
                    double sublead_photon_phi = dipho->subLeadingPhoton()->phi();
                    //double sublead_photon_eta = dipho->subLeadingPhoton()->superCluster()->eta(); // should use subLeadingPhoton()->eta() instead. Not sure why they use supercluster in jet dR matching ...
                    //double sublead_photon_phi = dipho->subLeadingPhoton()->superCluster()->phi();
                    for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                        int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                        if (abs(pdgid) != 11) continue;
                        if (genParticles->ptrAt( genLoop )->p4().pt() < 15) continue;
                        if (!genParticles->ptrAt( genLoop )->isPromptFinalState()) continue;
                        double electron_eta = genParticles->ptrAt( genLoop )->p4().eta();
                        double electron_phi = genParticles->ptrAt( genLoop )->p4().phi();

                        cout << "Found an electron with pT: " << genParticles->ptrAt( genLoop )->p4().pt() << endl;


                        double deltaR_lead = sqrt( pow(electron_eta - lead_photon_eta, 2) + pow(electron_phi - lead_photon_phi, 2));
                        double deltaR_sublead = sqrt( pow(electron_eta - sublead_photon_eta, 2) + pow(electron_phi - sublead_photon_phi, 2));

                        cout << "Delta R with leading photon: " << deltaR_lead << endl;
                        cout << "Delta R with subleading photon: " << deltaR_sublead << endl;

                        const double deltaR_thresh = 0.1;
                        if (deltaR_lead < deltaR_thresh)
                            lead_photon_is_electron = true;
                        if (deltaR_sublead < deltaR_thresh)
                            sublead_photon_is_electron = true;

                    } // end gen loop

                    int lead_photon_type;
                    if (dipho->leadingPhoton()->genMatchType() != 1) { // fake
                        if (!lead_photon_is_electron)
                            lead_photon_type = 3;   // fake
                        else
                            lead_photon_type = 2;   // fake from  electron 
                    }
                    else
                        lead_photon_type = 1; // prompt
                    int sublead_photon_type;
                    if (dipho->subLeadingPhoton()->genMatchType() != 1) { // fake
                        if (!sublead_photon_is_electron)
                            sublead_photon_type = 3;   // fake
                        else
                            sublead_photon_type = 2;   // fake from  electron 
                    }
                    else
                        sublead_photon_type = 1; // prompt

                    // Do our own gen-matching
                    double lead_photon_closest_match = dipho->leadingPhoton()->genMatchType() == 1 ? sqrt( pow(lead_photon_eta - dipho->leadingPhoton()->matchedGenPhoton()->eta(), 2) + pow(lead_photon_phi - dipho->leadingPhoton()->matchedGenPhoton()->phi(), 2)) : 999;
                    double sublead_photon_closest_match = dipho->subLeadingPhoton()->genMatchType() == 1 ? sqrt( pow(sublead_photon_eta - dipho->subLeadingPhoton()->matchedGenPhoton()->eta(), 2) + pow(sublead_photon_phi - dipho->subLeadingPhoton()->matchedGenPhoton()->phi(), 2)) : 999;

                    double lead_photon_closest_match_pt = dipho->leadingPhoton()->genMatchType() == 1 ? dipho->leadingPhoton()->pt() : -1;
                    double sublead_photon_closest_match_pt = dipho->subLeadingPhoton()->genMatchType() == 1 ? dipho->subLeadingPhoton()->pt() : -1;
                    for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                        int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                        if (abs(pdgid) != 22) continue;
                        if (genParticles->ptrAt( genLoop )->p4().pt() < 10) continue;
                        if (!genParticles->ptrAt( genLoop )->isPromptFinalState()) continue;
                        double gen_photon_candidate_eta = genParticles->ptrAt( genLoop )->p4().eta();
                        double gen_photon_candidate_phi = genParticles->ptrAt( genLoop )->p4().phi();

                        double deltaR_lead = sqrt( pow(gen_photon_candidate_eta - lead_photon_eta, 2) + pow(gen_photon_candidate_phi - lead_photon_phi, 2));
                        double deltaR_sublead = sqrt( pow(gen_photon_candidate_eta - sublead_photon_eta, 2) + pow(gen_photon_candidate_phi - sublead_photon_phi, 2));

                        if (deltaR_lead < lead_photon_closest_match) {   
                            lead_photon_closest_match = deltaR_lead;
                            lead_photon_closest_match_pt = genParticles->ptrAt( genLoop )->p4().pt();
                        }
                        if (deltaR_sublead < sublead_photon_closest_match) {
                            sublead_photon_closest_match = deltaR_sublead;
                            sublead_photon_closest_match_pt = genParticles->ptrAt( genLoop )->p4().pt();
                        }
                    } // end gen loop
                    string lead_fake = dipho->leadingPhoton()->genMatchType() == 1 ? "prompt" : "fake";
                    //cout << "Leading photon is " << lead_fake << ". Closest match:" << lead_photon_closest_match << endl;

                    string sublead_fake = dipho->subLeadingPhoton()->genMatchType() == 1 ? "prompt" : "fake";
                    //cout << "Subleading photon is " << sublead_fake << ". Closest match:" << sublead_photon_closest_match << endl;


                    tthhtags_obj.setnGoodEls(nGoodEls);
                    tthhtags_obj.setnGoodElsFromTau(nGoodElsFromTau);
                    tthhtags_obj.setnGoodMus(nGoodMus);
                    tthhtags_obj.setnGoodMusFromTau(nGoodMusFromTau);
                    tthhtags_obj.setnGoodTaus(nGoodTaus);
                    tthhtags_obj.setLeadPhotonType(lead_photon_type);
                    tthhtags_obj.setSubleadPhotonType(sublead_photon_type); 
                    tthhtags_obj.setLeadPhotonClosestDeltaR(lead_photon_closest_match);
                    tthhtags_obj.setSubleadPhotonClosestDeltaR(sublead_photon_closest_match);
                    tthhtags_obj.setLeadPhotonClosestPt(lead_photon_closest_match_pt);
                    tthhtags_obj.setSubleadPhotonClosestPt(sublead_photon_closest_match_pt);
                }
                else {
                    tthhtags_obj.setnGoodEls(-1);
                    tthhtags_obj.setnGoodElsFromTau(-1);
                    tthhtags_obj.setnGoodMus(-1);
                    tthhtags_obj.setnGoodMusFromTau(-1);
                    tthhtags_obj.setnGoodTaus(-1);
                    tthhtags_obj.setLeadPhotonType(-1);
                    tthhtags_obj.setSubleadPhotonType(-1);
                    tthhtags_obj.setLeadPhotonClosestDeltaR(-1);
                    tthhtags_obj.setSubleadPhotonClosestDeltaR(-1);
                    tthhtags_obj.setLeadPhotonClosestPt(-1);
                    tthhtags_obj.setSubleadPhotonClosestPt(-1);
                }


                if(!useTTHHadronicMVA_){
                    for( unsigned num = 0; num < JetVect.size(); num++ ) {
                        tthhtags_obj.includeWeightsByLabel( *JetVect[num] , "JetBTagCutWeight");
                        //tthhtags_obj.includeWeightsByLabel( *JetVect[num] , "JetBTagReshapeWeight");
                    }
                } else {
                    for( unsigned num = 0; num < JetVect.size(); num++ ) {
                        tthhtags_obj.includeWeightsByLabel( *JetVect[num] , "JetBTagReshapeWeight");
                    }                    
                }
                tthhtags_obj.includeWeights( *dipho );
                tthhtags->push_back( tthhtags_obj );
                if( ! evt.isRealData() ) {
                    TagTruthBase truth_obj;
                    truth_obj.setGenPV( higgsVtx );
                    if ( stage0cat.isValid() ) {
                        truth_obj.setHTXSInfo( *( stage0cat.product() ),
                                               *( stage1cat.product() ),
                                               *( njets.product() ),
                                               *( pTH.product() ),
                                               *( pTV.product() ) );
                    } else {
                        truth_obj.setHTXSInfo( 0, 0, 0, 0., 0. );
                    }
                    truths->push_back( truth_obj );
                    tthhtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, idx++ ) ) );
                }
                // count++;
            }
        }
        evt.put( std::move( tthhtags ) );
        evt.put( std::move( truths ) );
        // cout << "tagged events = " << count << endl;
    }
}
typedef flashgg::TTHHadronicTagProducer FlashggTTHHadronicTagProducer;
DEFINE_FWK_MODULE( FlashggTTHHadronicTagProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

