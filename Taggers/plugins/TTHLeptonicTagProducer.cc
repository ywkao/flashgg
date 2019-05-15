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
#include "flashgg/DataFormats/interface/TTHLeptonicTag.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "flashgg/Taggers/interface/LeptonSelection2018.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"

#include "flashgg/Taggers/interface/BDT_resolvedTopTagger.h"

#include <vector>
#include <algorithm>
#include <string>
#include <utility>
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMVA/Reader.h"
#include "TRandom.h"

TRandom* myRandLeptonic = new TRandom();

using namespace std;
using namespace edm;


namespace flashgg {
    class TTHLeptonicTagProducer : public EDProducer
    {

    public:
        typedef math::XYZPoint Point;

        TTHLeptonicTagProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;
        int  chooseCategory( float tthmvavalue , bool debug_ );

        const  reco::GenParticle* motherID(const reco::GenParticle* gp);
        bool PassFrixione(Handle<View<reco::GenParticle> > genParticles, const reco::GenParticle* gp, int nBinsForFrix, double cone_frix);
        vector<int> IsPromptAfterOverlapRemove(Handle<View<reco::GenParticle> > genParticles, const edm::Ptr<reco::GenParticle> genPho);
        int  GenPhoIndex(Handle<View<reco::GenParticle> > genParticles, const flashgg::Photon* pho, int usedIndex);
        double NearestDr(Handle<View<reco::GenParticle> > genParticles, const reco::GenParticle* gp);

        struct Sorter {
            bool operator()( const std::pair<unsigned int, float>  pair1, const std::pair<unsigned int,float>  pair2 )
            {
                return ( pair1.second > pair2.second );
            };
        };
        

        std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > tokenJets_;
        EDGetTokenT<View<DiPhotonCandidate> > diPhotonToken_;
        //EDGetTokenT<View<Jet> > thejetToken_;
        std::vector<edm::InputTag> inputTagJets_;
        EDGetTokenT<View<Electron> > electronToken_;
        EDGetTokenT<View<flashgg::Muon> > muonToken_;
        EDGetTokenT<View<flashgg::Met> > METToken_;
        EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;
        EDGetTokenT<View<Photon> > photonToken_;
        EDGetTokenT<View<reco::Vertex> > vertexToken_;
        EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
        EDGetTokenT<int> stage0catToken_, stage1catToken_, njetsToken_;
        EDGetTokenT<HTXS::HiggsClassification> newHTXSToken_;
        EDGetTokenT<float> pTHToken_,pTVToken_;
        string systLabel_;
        FileInPath topTaggerXMLfile_;

        typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

        unique_ptr<TMVA::Reader> DiphotonMva_;
        FileInPath MVAweightfile_;

        //Thresholds

        double MuonEtaCut_;
        double MuonPtCut_;
        double MuonIsoCut_;
        double MuonPhotonDrCut_;
        int    MinNLep_;
        int    MaxNLep_;

        double ElePtCut_;
        std::vector<double> EleEtaCuts_;
        double ElePhotonDrCut_;
        double ElePhotonZMassCut_;
        double DeltaRTrkEle_;

        double LeptonsZMassCut_;

        double DiLeptonJetThreshold_;
        double DiLeptonbJetThreshold_;
        double DiLeptonMVAThreshold_;

        double leadPhoOverMassThreshold_;
        double subleadPhoOverMassThreshold_;
        vector<double> MVAThreshold_;
        double deltaRJetLeadPhoThreshold_;
        double deltaRJetSubLeadPhoThreshold_;
        double jetsNumberThreshold_;
        double bjetsNumberThreshold_;
        double jetPtThreshold_;
        double jetEtaThreshold_;
        double deltaRJetLepton_;
        double leadingJetPtThreshold_;
        vector<double> bDiscriminator_;
        string bTag_;
        double PhoMVAThreshold_;

        bool UseCutBasedDiphoId_;
        bool debug_;
        bool SplitDiLeptEv_;
        vector<double> CutBasedDiphoId_;

        float leadeta_;
        float subleadeta_;
        float leadptom_;
        float subleadptom_;
        float leadIDMVA_;
        float subleadIDMVA_;
        float deltaphi_;
        float leadPSV_;
        float subleadPSV_;
        float nJets_;
        float nJets_bTagMedium_;
        float jet_pt1_;
        float jet_pt2_;
        float jet_pt3_;
        float jet_eta1_;
        float jet_eta2_;
        float jet_eta3_;
        float bTag1_;
        float bTag2_;
        float MetPt_;
        float MetPhi_;
        float lepton_leadPt_;
        float lepton_leadEta_;
        float lepton_leadPhi_;
        BDT_resolvedTopTagger *topTagger;
    };

    const reco::GenParticle* TTHLeptonicTagProducer::motherID(const reco::GenParticle* gp)
    {
        const reco::GenParticle* mom_lead = gp;
        //cout << "in id: " << gp->pdgId() << endl;
        while( mom_lead->numberOfMothers() > 0 ) {
             for(uint j = 0; j < mom_lead->numberOfMothers(); ++j) {
                 mom_lead = dynamic_cast<const reco::GenParticle*>( mom_lead->mother(j) );
                 //cout << j << "th id: " << mom_lead->pdgId() << ", gpid: " << gp->pdgId() << endl;
                 if( mom_lead->pdgId() != gp->pdgId() )
                     return mom_lead; 
                     //break;
                }
             }

        return mom_lead;
    }

    bool TTHLeptonicTagProducer::PassFrixione(Handle<View<reco::GenParticle> > genParticles, const reco::GenParticle* gp, int nBinsForFrix, double cone_frix)
    {

        bool passFrix = true;

        double pho_et = gp->p4().Et();
        double pho_eta = gp->p4().eta();
        double pho_phi = gp->p4().phi();
        const double initConeFrix = 1E-10;
        double ets[nBinsForFrix] = {};
   
        for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ )
           {
           int status = genParticles->ptrAt( genLoop )->status();
           if (!(status >= 21 && status <=24)) continue;
           int pdgid = genParticles->ptrAt( genLoop )->pdgId();
           if (!(abs(pdgid) == 11 || abs(pdgid) == 13 || abs(pdgid) == 15 || abs(pdgid) < 10 || abs(pdgid) == 21) ) continue; 
           double et = genParticles->ptrAt( genLoop )->p4().Et();
           double eta = genParticles->ptrAt( genLoop )->p4().eta();
           double phi = genParticles->ptrAt( genLoop )->p4().phi();
           double dR = deltaR(pho_eta,pho_phi,eta,phi); 
           for (int j = 0; j < nBinsForFrix; j++) { 
               double cone_var = (initConeFrix + j*cone_frix/nBinsForFrix);
               if (dR < cone_var && cone_var < cone_frix) ets[j] += et/(1-cos(cone_var) )*(1-cos(cone_frix) );
               }
           }

           for (int k = 0; k < nBinsForFrix; k++) {
               if (ets[k] > pho_et) {
                  passFrix = false; break;
                  }
           }

        return passFrix;
    }


    double TTHLeptonicTagProducer::NearestDr(Handle<View<reco::GenParticle> > genParticles, const reco::GenParticle* gp)
    {

        double nearDr = 999;

        double pho_eta = gp->p4().eta();
        double pho_phi = gp->p4().phi();

        for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ )
           {
           int status = genParticles->ptrAt( genLoop )->status();
           if (!(status >= 21 && status <=24)) continue;
           int pdgid = genParticles->ptrAt( genLoop )->pdgId();
           if (!(abs(pdgid) == 11 || abs(pdgid) == 13 || abs(pdgid) == 15 || abs(pdgid) < 10 || abs(pdgid) == 21) ) continue;
           double eta = genParticles->ptrAt( genLoop )->p4().eta();
           double phi = genParticles->ptrAt( genLoop )->p4().phi();
           double dR = deltaR(pho_eta,pho_phi,eta,phi);
           if (dR < nearDr) {
              nearDr = dR;
              }
           }

        return nearDr;
    }


    vector<int> TTHLeptonicTagProducer::IsPromptAfterOverlapRemove(Handle<View<reco::GenParticle> > genParticles, const edm::Ptr<reco::GenParticle> genPho)
    {

        vector<int> flags;
        bool isPromptAfterOverlapRemove = false;

        int simplemotherID = -1;
        int simplemotherStatus = -1;
        if (genPho->numberOfMothers() > 0) {
            simplemotherID = genPho->mother(0)->pdgId();
            simplemotherStatus = genPho->mother(0)->status();
        }
        const reco::GenParticle* mom = motherID(&(*genPho));
        const reco::GenParticle* mommom = motherID(mom);
        bool isMad = false;
        if (simplemotherID == 22 && simplemotherStatus >= 21 && simplemotherStatus <= 24) isMad = true;
        bool isFromWb = false;
        if (abs(mom->pdgId()) == 24 || abs(mommom->pdgId()) == 24 || (abs(mom->pdgId()) == 5 && abs(mommom->pdgId()) == 6) ) isFromWb = true;
        bool isFromQuark = false;
        if (abs(mom->pdgId()) == 21 || abs(mommom->pdgId()) == 21 || (abs(mom->pdgId()) <= 6 && abs(mommom->pdgId()) != 6 && abs(mommom->pdgId()) != 24 ) ) isFromQuark = true;
        bool isFromProton = false;
        if (abs(mom->pdgId()) == 2212) isFromProton = true;
        bool failFrix = !PassFrixione(genParticles, &(*genPho), 100, 0.05);
        bool isPythia = false;
        if (!isMad && genPho->isPromptFinalState() && ( isFromWb || (isFromQuark && failFrix) || isFromProton) ) isPythia = true;
        if (isMad || isPythia) isPromptAfterOverlapRemove = true;
        flags.push_back(isPromptAfterOverlapRemove ? 1 : 0);
        flags.push_back(isMad ? 1 : 0);
        flags.push_back(isPythia ? 1 : 0);
        flags.push_back((!failFrix) ? 1 : 0);
        flags.push_back(simplemotherID);
        flags.push_back(simplemotherStatus);
        flags.push_back(abs(mom->pdgId()));
        flags.push_back(abs(mommom->pdgId()));
        return flags;

    }


    int TTHLeptonicTagProducer::GenPhoIndex(Handle<View<reco::GenParticle> > genParticles, const flashgg::Photon* pho, int usedIndex)
    {
        double maxDr = 0.2;
        double ptDiffMax = 99e15;
        int index = -1;
                
        for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
            int pdgid = genParticles->ptrAt( genLoop )->pdgId();
            if (pho->genMatchType() != 1) continue;
            if (int(genLoop) == usedIndex) continue;
            if (abs(pdgid) != 22) continue;
            if (genParticles->ptrAt( genLoop )->p4().pt() < 10) continue;
            if (!genParticles->ptrAt( genLoop )->isPromptFinalState()) continue;

            double gen_photon_candidate_pt = genParticles->ptrAt( genLoop )->p4().pt();
            double gen_photon_candidate_eta = genParticles->ptrAt( genLoop )->p4().eta();
            double gen_photon_candidate_phi = genParticles->ptrAt( genLoop )->p4().phi();
            double deltaR_ = deltaR(gen_photon_candidate_eta, gen_photon_candidate_phi, pho->p4().eta(), pho->p4().phi());

            if (deltaR_ > maxDr) continue;

            double ptdiff = abs(gen_photon_candidate_pt - pho->p4().pt());
            if (ptdiff < ptDiffMax) {
                ptDiffMax = ptdiff;
                index = genLoop;
            }
            
        } // end gen loop

        return index;
    }



    TTHLeptonicTagProducer::TTHLeptonicTagProducer( const ParameterSet &iConfig ) :
        diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
        //thejetToken_( consumes<View<flashgg::Jet> >( iConfig.getParameter<InputTag>( "JetTag" ) ) ),
        inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
        electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
        muonToken_( consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
        METToken_( consumes<View<flashgg::Met> >( iConfig.getParameter<InputTag>( "MetTag" ) ) ),
        mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag> ( "MVAResultTag" ) ) ),
        vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
        genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
        systLabel_( iConfig.getParameter<string> ( "SystLabel" ) )
    {
        leadPhoOverMassThreshold_ = iConfig.getParameter<double>( "leadPhoOverMassThreshold");
        subleadPhoOverMassThreshold_ = iConfig.getParameter<double>( "subleadPhoOverMassThreshold");
        MVAThreshold_ = iConfig.getParameter<std::vector<double>>( "MVAThreshold");
        PhoMVAThreshold_ = iConfig.getParameter<double>( "PhoMVAThreshold");
        jetsNumberThreshold_ = iConfig.getParameter<double>( "jetsNumberThreshold");
        bjetsNumberThreshold_ = iConfig.getParameter<double>( "bjetsNumberThreshold");
        jetPtThreshold_ = iConfig.getParameter<double>( "jetPtThreshold");
        jetEtaThreshold_ = iConfig.getParameter<double>( "jetEtaThreshold");
        deltaRJetLepton_ = iConfig.getParameter<double>( "deltaRJetLepton");
        leadingJetPtThreshold_ = iConfig.getParameter<double>("leadingJetPtThreshold");

        MuonEtaCut_ = iConfig.getParameter<double>( "MuonEtaCut");
        MuonPtCut_ = iConfig.getParameter<double>( "MuonPtCut");
        MuonIsoCut_ = iConfig.getParameter<double>( "MuonIsoCut");
        MuonPhotonDrCut_ = iConfig.getParameter<double>( "MuonPhotonDrCut");

        MinNLep_ = iConfig.getParameter<int>( "MinNLep");
        MaxNLep_ = iConfig.getParameter<int>( "MaxNLep");
 
        EleEtaCuts_ = iConfig.getParameter<std::vector<double>>( "EleEtaCuts");
        ElePtCut_ = iConfig.getParameter<double>( "ElePtCut");
        ElePhotonDrCut_ = iConfig.getParameter<double>( "ElePhotonDrCut");
        ElePhotonZMassCut_ = iConfig.getParameter<double>( "ElePhotonZMassCut");
        DeltaRTrkEle_ = iConfig.getParameter<double>( "DeltaRTrkEle");

        LeptonsZMassCut_ = iConfig.getParameter<double>( "LeptonsZMassCut");

        DiLeptonJetThreshold_ = iConfig.getParameter<double>( "DiLeptonJetThreshold");
        DiLeptonbJetThreshold_ = iConfig.getParameter<double>( "DiLeptonbJetThreshold");
        DiLeptonMVAThreshold_ = iConfig.getParameter<double>( "DiLeptonMVAThreshold");

        deltaRJetLeadPhoThreshold_ = iConfig.getParameter<double>( "deltaRJetLeadPhoThreshold");
        deltaRJetSubLeadPhoThreshold_ = iConfig.getParameter<double>( "deltaRJetSubLeadPhoThreshold");

        bDiscriminator_ = iConfig.getParameter<vector<double > >( "bDiscriminator");
        bTag_ = iConfig.getParameter<string>( "bTag");

        UseCutBasedDiphoId_ = iConfig.getParameter<bool>( "UseCutBasedDiphoId" );
        debug_ = iConfig.getParameter<bool>( "debug" );
        SplitDiLeptEv_ = iConfig.getParameter<bool>( "SplitDiLeptEv" );
        CutBasedDiphoId_ = iConfig.getParameter<std::vector<double>>( "CutBasedDiphoId" );

        ParameterSet HTXSps = iConfig.getParameterSet( "HTXSTags" );
        stage0catToken_ = consumes<int>( HTXSps.getParameter<InputTag>("stage0cat") );
        stage1catToken_ = consumes<int>( HTXSps.getParameter<InputTag>("stage1cat") );
        njetsToken_ = consumes<int>( HTXSps.getParameter<InputTag>("njets") );
        pTHToken_ = consumes<float>( HTXSps.getParameter<InputTag>("pTH") );
        pTVToken_ = consumes<float>( HTXSps.getParameter<InputTag>("pTV") );
        newHTXSToken_ = consumes<HTXS::HiggsClassification>( HTXSps.getParameter<InputTag>("ClassificationObj") );

        MVAweightfile_ = iConfig.getParameter<edm::FileInPath>( "MVAweightfile" );
        topTaggerXMLfile_ = iConfig.getParameter<edm::FileInPath>( "topTaggerXMLfile" );

        DiphotonMva_.reset( new TMVA::Reader( "!Color:Silent" ) );
        DiphotonMva_->AddVariable( "dipho_leadEta", &leadeta_ );
        DiphotonMva_->AddVariable( "dipho_subleadEta", &subleadeta_ );
        DiphotonMva_->AddVariable( "dipho_lead_ptoM", &leadptom_ );
        DiphotonMva_->AddVariable( "dipho_sublead_ptoM", &subleadptom_ );
        DiphotonMva_->AddVariable( "dipho_leadIDMVA", &leadIDMVA_ );
        DiphotonMva_->AddVariable( "dipho_subleadIDMVA", &subleadIDMVA_ );
        DiphotonMva_->AddVariable( "dipho_deltaphi", &deltaphi_ );
        DiphotonMva_->AddVariable( "dipho_lead_PSV", &leadPSV_ );
        DiphotonMva_->AddVariable( "dipho_sublead_PSV", &subleadPSV_ );
        DiphotonMva_->AddVariable( "nJets", &nJets_ );
        DiphotonMva_->AddVariable( "nJets_bTagMedium", &nJets_bTagMedium_ );
        DiphotonMva_->AddVariable( "jet1_pt", &jet_pt1_ );
        DiphotonMva_->AddVariable( "jet2_pt", &jet_pt2_ );
        DiphotonMva_->AddVariable( "jet3_pt", &jet_pt3_ );
        DiphotonMva_->AddVariable( "jet1_eta", &jet_eta1_ );
        DiphotonMva_->AddVariable( "jet2_eta", &jet_eta2_ );
        DiphotonMva_->AddVariable( "jet3_eta", &jet_eta3_ );
        DiphotonMva_->AddVariable( "bTag1", &bTag1_ );
        DiphotonMva_->AddVariable( "bTag2", &bTag2_ );
        DiphotonMva_->AddVariable( "MetPt", &MetPt_ );
        DiphotonMva_->AddVariable( "lepton_leadPt", &lepton_leadPt_ );
        DiphotonMva_->AddVariable( "lepton_leadEta", &lepton_leadEta_ );

        DiphotonMva_->BookMVA( "BDT", MVAweightfile_.fullPath() );

        for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
            auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
            tokenJets_.push_back(token);
        }

        produces<vector<TTHLeptonicTag> >();
        produces<vector<TagTruthBase> >();

        topTagger = new BDT_resolvedTopTagger(topTaggerXMLfile_.fullPath());
    }

    int TTHLeptonicTagProducer::chooseCategory( float tthmvavalue , bool debug_)
    {
        // should return 0 if mva above all the numbers, 1 if below the first, ..., boundaries.size()-N if below the Nth, ...
        int n;
        for( n = 0 ; n < ( int )MVAThreshold_.size() ; n++ ) {
            if( ( double )tthmvavalue > MVAThreshold_[MVAThreshold_.size() - n - 1] ) { return n; }
        }

      if(debug_)
        {   cout << "Checking class, thresholds: ";
            for(unsigned int i=0; i<MVAThreshold_.size(); ++i)
                cout << MVAThreshold_[i] << " ";
        }

        return -1; // Does not pass, object will not be produced
    }

    void TTHLeptonicTagProducer::produce( Event &evt, const EventSetup & )
    {
        Handle<int> stage0cat, stage1cat, njets;
        Handle<float> pTH, pTV;
        evt.getByToken(stage0catToken_, stage0cat);
        evt.getByToken(stage1catToken_,stage1cat);
        evt.getByToken(njetsToken_,njets);
        evt.getByToken(pTHToken_,pTH);
        evt.getByToken(pTVToken_,pTV);
        Handle<HTXS::HiggsClassification> htxsClassification;
        evt.getByToken(newHTXSToken_,htxsClassification);


        //Handle<View<flashgg::Jet> > theJets;
        //evt.getByToken( thejetToken_, theJets );
        //const PtrVector<flashgg::Jet>& jetPointers = theJets->ptrVector();
        JetCollectionVector Jets( inputTagJets_.size() );
        for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
            evt.getByToken( tokenJets_[j], Jets[j] );
        }

        Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
        evt.getByToken( diPhotonToken_, diPhotons );

        Handle<View<flashgg::Muon> > theMuons;
        evt.getByToken( muonToken_, theMuons );

        Handle<View<flashgg::Electron> > theElectrons;
        evt.getByToken( electronToken_, theElectrons );

        Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;
        evt.getByToken( mvaResultToken_, mvaResults );

        Handle<View<reco::GenParticle> > genParticles;

        Handle<View<reco::Vertex> > vertices;
        evt.getByToken( vertexToken_, vertices );

        Handle<View<flashgg::Met> > theMet_;
        evt.getByToken( METToken_, theMet_ );


        std::unique_ptr<vector<TTHLeptonicTag> > tthltags( new vector<TTHLeptonicTag> );
        std::unique_ptr<vector<TagTruthBase> > truths( new vector<TagTruthBase> );
        edm::RefProd<vector<TagTruthBase> > rTagTruth = evt.getRefBeforePut<vector<TagTruthBase> >();
        unsigned int idx = 0;

        Point higgsVtx;

        if( ! evt.isRealData() )
        {
            evt.getByToken( genParticleToken_, genParticles );
            for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ )
            {
                int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                if( pdgid == 25 || pdgid == 22 )
                {
                    higgsVtx = genParticles->ptrAt( genLoop )->vertex();
                    break;
                }
            }
        }

        assert( diPhotons->size() == mvaResults->size() );

        double idmva1 = 0.;
        double idmva2 = 0.;

        for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ )
        {
            unsigned int jetCollectionIndex = diPhotons->ptrAt( diphoIndex )->jetCollectionIndex();

            edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );
            edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );

            if( dipho->leadingPhoton()->pt() < ( dipho->mass() )*leadPhoOverMassThreshold_ ) { continue; }
            if( dipho->subLeadingPhoton()->pt() < ( dipho->mass() )*subleadPhoOverMassThreshold_ ) { continue; }
            idmva1 = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
            idmva2 = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );

            if(debug_)
                cout << "Photon pair with PhoIdMVA values: " << idmva1 << " " << idmva2 << endl;
 
            if( idmva1 < PhoMVAThreshold_ || idmva2 < PhoMVAThreshold_ ) { continue; }

            bool passDiphotonSelection = true;
            if(UseCutBasedDiphoId_)
            {
                assert(CutBasedDiphoId_.size()==6);
                if(dipho->leadingPhoton()->pt()/dipho->mass() < CutBasedDiphoId_[0]) passDiphotonSelection = false;
                if(dipho->subLeadingPhoton()->pt()/dipho->mass() < CutBasedDiphoId_[1]) passDiphotonSelection = false;
                if(dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() ) < CutBasedDiphoId_[2]) passDiphotonSelection = false;
                if(dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() ) < CutBasedDiphoId_[3]) passDiphotonSelection = false;
                if(abs (dipho->leadingPhoton()->eta() - dipho->subLeadingPhoton()->eta()) > CutBasedDiphoId_[4]) passDiphotonSelection = false;
                if(deltaPhi(dipho->leadingPhoton()->phi(), dipho->subLeadingPhoton()->phi() ) > CutBasedDiphoId_[5] ) passDiphotonSelection = false;
            }

            if(!passDiphotonSelection) continue;

            if(debug_)
                cout << "Passed photon selection, checking leptons: " << idmva1 << " " << idmva2 << endl;

            std::vector<edm::Ptr<flashgg::Muon> >     Muons;
            std::vector<edm::Ptr<flashgg::Muon> >     MuonsLoose;
            std::vector<edm::Ptr<flashgg::Muon> >     MuonsMedium;
            std::vector<edm::Ptr<flashgg::Muon> >     MuonsTight;

            std::vector<edm::Ptr<flashgg::Electron> > Electrons;
            std::vector<edm::Ptr<flashgg::Electron> > ElectronsLoose;
            std::vector<edm::Ptr<flashgg::Electron> > ElectronsMedium;
            std::vector<edm::Ptr<flashgg::Electron> > ElectronsTight;

            std::vector<double> lepPt;
            std::vector<double> lepEta;
            std::vector<double> lepPhi;
            std::vector<double> lepE;
            std::vector<int>    lepType;

            if(theMuons->size()>0) {
                Muons = selectMuons(theMuons->ptrs(), dipho, vertices->ptrs(), MuonPtCut_, MuonEtaCut_, MuonIsoCut_, MuonPhotonDrCut_, debug_);
                MuonsTight = selectMuons(theMuons->ptrs(), dipho, vertices->ptrs(), MuonPtCut_, MuonEtaCut_, MuonIsoCut_, MuonPhotonDrCut_, debug_, 3);
            }
            if(theElectrons->size()>0) {
                Electrons = selectElectrons(theElectrons->ptrs(), dipho, ElePtCut_, EleEtaCuts_, ElePhotonDrCut_, ElePhotonZMassCut_, DeltaRTrkEle_, debug_);
                ElectronsTight = selectElectrons(theElectrons->ptrs(), dipho, ElePtCut_, EleEtaCuts_, ElePhotonDrCut_, ElePhotonZMassCut_, DeltaRTrkEle_, debug_, 3);
            } 


            //If 2 same flavour leptons are found remove the pairs with mass compatible with a Z boson

            if(Muons.size()>=2)
            {
                std::vector<edm::Ptr<flashgg::Muon>> Muons_0;
                Muons_0 = Muons;
                std::vector<int> badIndexes;

                for(unsigned int i=0; i<Muons_0.size(); ++i)
                {
                    for(unsigned int j=i+1; j<Muons_0.size(); ++j)
                    {
                        TLorentzVector l1, l2;
                        l1.SetPtEtaPhiE(Muons_0[i]->pt(), Muons_0[i]->eta(), Muons_0[i]->phi(), Muons_0[i]->energy());
                        l2.SetPtEtaPhiE(Muons_0[j]->pt(), Muons_0[j]->eta(), Muons_0[j]->phi(), Muons_0[j]->energy());

                        if(fabs((l1+l2).M() - 91.187) < LeptonsZMassCut_)
                        {
                            badIndexes.push_back(i);
                            badIndexes.push_back(j);
                        }
                    }
                }

                if(badIndexes.size()!=0)
                {
                    Muons.clear();
                    for(unsigned int i=0; i<Muons_0.size(); ++i)
                    {
                       bool isBad = false;
                       for(unsigned int j=0; j<badIndexes.size(); ++j)
                       {
                          if(badIndexes[j]==(int)i)
                               isBad = true;
                      }
                      if(!isBad) Muons.push_back(Muons_0[i]);
                    }
                }
            }        

            if(MuonsTight.size()>=2)
            {
                std::vector<edm::Ptr<flashgg::Muon>> MuonsTight_0;
                MuonsTight_0 = MuonsTight;
                std::vector<int> badIndexes;

                for(unsigned int i=0; i<MuonsTight_0.size(); ++i)
                {
                    for(unsigned int j=i+1; j<MuonsTight_0.size(); ++j)
                    {
                        TLorentzVector l1, l2;
                        l1.SetPtEtaPhiE(MuonsTight_0[i]->pt(), MuonsTight_0[i]->eta(), MuonsTight_0[i]->phi(), MuonsTight_0[i]->energy());
                        l2.SetPtEtaPhiE(MuonsTight_0[j]->pt(), MuonsTight_0[j]->eta(), MuonsTight_0[j]->phi(), MuonsTight_0[j]->energy());

                        if(fabs((l1+l2).M() - 91.187) < LeptonsZMassCut_)
                        {
                            badIndexes.push_back(i);
                            badIndexes.push_back(j);
                        }
                    }
                }

                if(badIndexes.size()!=0)
                {
                    MuonsTight.clear();
                    for(unsigned int i=0; i<MuonsTight_0.size(); ++i)
                    {
                       bool isBad = false;
                       for(unsigned int j=0; j<badIndexes.size(); ++j)
                       {
                          if(badIndexes[j]==(int)i)
                               isBad = true;
                      }
                      if(!isBad) MuonsTight.push_back(MuonsTight_0[i]);
                    }
                }
            }


            if(Electrons.size()>=2)
            {
                std::vector<int> badIndexes;
                std::vector<edm::Ptr<flashgg::Electron> > Electrons_0;
                Electrons_0 = Electrons;
                for(unsigned int i=0; i<Electrons_0.size(); ++i)
                {
                    for(unsigned int j=i+1; j<Electrons_0.size(); ++j)
                    {
                        TLorentzVector l1, l2;
                        l1.SetPtEtaPhiE(Electrons_0[i]->pt(), Electrons_0[i]->eta(), Electrons_0[i]->phi(), Electrons_0[i]->energy());
                        l2.SetPtEtaPhiE(Electrons_0[j]->pt(), Electrons_0[j]->eta(), Electrons_0[j]->phi(), Electrons_0[j]->energy());

                        if(fabs((l1+l2).M() - 91.187) < LeptonsZMassCut_)
                        {
                            badIndexes.push_back(i);
                            badIndexes.push_back(j);
                        }
                    }
                }
                if(badIndexes.size()!=0)
                {
                    Electrons.clear();

                    for(unsigned int i=0; i<Electrons_0.size(); ++i)
                    {
                         bool isBad = false;
                         for(unsigned int j=0; j<badIndexes.size(); ++j)
                         {
                             if(badIndexes[j]==(int)i)
                                 isBad = true;
                         }
                         if(!isBad) Electrons.push_back(Electrons_0[i]);
                    }
                 }
             }        

            if(ElectronsTight.size()>=2)
            {
                std::vector<int> badIndexes;
                std::vector<edm::Ptr<flashgg::Electron> > ElectronsTight_0;
                ElectronsTight_0 = ElectronsTight;
                for(unsigned int i=0; i<ElectronsTight_0.size(); ++i)
                {
                    for(unsigned int j=i+1; j<ElectronsTight_0.size(); ++j)
                    {
                        TLorentzVector l1, l2;
                        l1.SetPtEtaPhiE(ElectronsTight_0[i]->pt(), ElectronsTight_0[i]->eta(), ElectronsTight_0[i]->phi(), ElectronsTight_0[i]->energy());
                        l2.SetPtEtaPhiE(ElectronsTight_0[j]->pt(), ElectronsTight_0[j]->eta(), ElectronsTight_0[j]->phi(), ElectronsTight_0[j]->energy());

                        if(fabs((l1+l2).M() - 91.187) < LeptonsZMassCut_)
                        {
                            badIndexes.push_back(i);
                            badIndexes.push_back(j);
                        }
                    }
                }
                if(badIndexes.size()!=0)
                {
                    ElectronsTight.clear();

                    for(unsigned int i=0; i<ElectronsTight_0.size(); ++i)
                    {
                         bool isBad = false;
                         for(unsigned int j=0; j<badIndexes.size(); ++j)
                         {
                             if(badIndexes[j]==(int)i)
                                 isBad = true;
                         }
                         if(!isBad) ElectronsTight.push_back(ElectronsTight_0[i]);
                    }
                 }
             }

            if( (Muons.size() + Electrons.size()) < (unsigned) MinNLep_ || (Muons.size() + Electrons.size()) > (unsigned) MaxNLep_) continue;

            // Fill lepton vectors            
            // ===================
            
            std::vector<std::pair< unsigned int, float > > sorter;

            if(debug_) cout<<" nMuons="<<Muons.size()<<" nElectrons="<<Electrons.size()<< endl;

            for(unsigned int i=0;i<Muons.size();i++){
                float pt=Muons[i]->pt();
                int index=100;
                index+=i;
                std::pair<unsigned int, float>pairToSort = std::make_pair(index, pt);
                sorter.push_back( pairToSort );
                if(debug_) cout<<" muon "<< i <<" pt="<<pt<< endl;
            }
            for(unsigned int i=0;i<Electrons.size();i++){
                float pt=Electrons[i]->pt();
                int index=200;
                index+=i;
                std::pair<unsigned int, float>pairToSort = std::make_pair(index, pt);
                sorter.push_back( pairToSort );
                if(debug_) cout<<" elec "<< i <<" pt="<<pt<< endl;
            }
            // sort map by pt

            std::sort( sorter.begin(), sorter.end(), Sorter() );
                
            // fill vectors
            for (unsigned int i=0;i<sorter.size();i++){

                if(debug_) cout<<" Filling lepton vector index:"<<sorter[i].first<<" pt:" <<sorter[i].second << endl;

                lepPt.push_back(sorter[i].second);                

                int type=0;
                type=(sorter[i].first)/100;
                int n=(sorter[i].first)%100;

                lepType.push_back(type);

                if(debug_) cout<<" type="<<type<<" n="<<n<< endl;
                    
                if(type==1){
                    if(debug_) cout<<"MUON LEPPTCHECK "<<   sorter[i].second<<" "<<Muons[n]->pt()<< endl;
                    lepEta.push_back(Muons[n]->eta());
                    lepPhi.push_back(Muons[n]->phi());
                    lepE.push_back(Muons[n]->energy());
                    
                }else if(type==2){
                    if(debug_) cout<<"ELEC LEPPTCHECK "<<   sorter[i].second<<" "<<Electrons[n]->pt()<< endl;
                    lepEta.push_back(Electrons[n]->eta());
                    lepPhi.push_back(Electrons[n]->phi());
                    lepE.push_back(Electrons[n]->energy());
                }                
            }
            
            int njet_ = 0;
            int njets_btagloose_ = 0;
            int njets_btagmedium_ = 0;
            int njets_btagtight_ = 0;
            std::vector<edm::Ptr<flashgg::Jet>> tagJets;
            std::vector<edm::Ptr<flashgg::Jet>> tagBJets;
            std::vector<float> bTags;

            for( unsigned int jetIndex = 0; jetIndex < Jets[jetCollectionIndex]->size() ; jetIndex++ )
            {
                edm::Ptr<flashgg::Jet> thejet = Jets[jetCollectionIndex]->ptrAt( jetIndex );

                if( fabs( thejet->eta() ) > jetEtaThreshold_ ) { continue; }
                if(!thejet->passesJetID  ( flashgg::Loose ) ) { continue; }
                if( thejet->pt() < jetPtThreshold_ ) { continue; }

                float dRPhoLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->leadingPhoton()->eta(), dipho->leadingPhoton()->phi() ) ;
                float dRPhoSubLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->subLeadingPhoton()->eta(), dipho->subLeadingPhoton()->phi() );

                if( dRPhoLeadJet < deltaRJetLeadPhoThreshold_ || dRPhoSubLeadJet < deltaRJetSubLeadPhoThreshold_ ) { continue; }

                bool passDrLeptons = 1;                

                for( unsigned int muonIndex = 0; muonIndex < Muons.size(); muonIndex++ )
                {
                    Ptr<flashgg::Muon> muon = Muons[muonIndex];

                    float dRLept = deltaR( thejet->eta(), thejet->phi(), muon->eta(),  muon->phi() ) ;
                    if( dRLept > deltaRJetLepton_) { continue; }

                    passDrLeptons = 0;
                    break;
                }

                if(passDrLeptons)
                {   for( unsigned int eleIndex = 0; eleIndex < Electrons.size(); eleIndex++ )
                    {
                        Ptr<flashgg::Electron> ele = Electrons[eleIndex];
                    
                        float dRLept = deltaR( thejet->eta(), thejet->phi(), ele->eta(),  ele->phi() ) ;
                        if( dRLept > deltaRJetLepton_) { continue; }

                        passDrLeptons = 0;
                        break;
                    }
                }

                if(passDrLeptons)
                {
                    njet_++;
                    tagJets.push_back( thejet );
                    float bDiscriminatorValue = -2.;
                    if(bTag_ == "pfDeepCSV") bDiscriminatorValue = thejet->bDiscriminator("pfDeepCSVJetTags:probb")+thejet->bDiscriminator("pfDeepCSVJetTags:probbb") ;
                    else  bDiscriminatorValue = thejet->bDiscriminator( bTag_ );

                    bDiscriminatorValue >= 0. ? bTags.push_back(bDiscriminatorValue) : bTags.push_back(-1.);

                    if( bDiscriminatorValue > bDiscriminator_[0] ) njets_btagloose_++;
                    if( bDiscriminatorValue > bDiscriminator_[1] ) njets_btagmedium_++;
                    if( bDiscriminatorValue > bDiscriminator_[2] ) njets_btagtight_++;

                    if( bDiscriminatorValue > bDiscriminator_[1] )
                        tagBJets.push_back( thejet );

                    bool eval_top_tagger = true;

                    if (eval_top_tagger) {
                      float cvsl = thejet->bDiscriminator("pfDeepCSVJetTags:probc") + thejet->bDiscriminator("pfDeepCSVJetTags:probudsg") ;
                      float cvsb = thejet->bDiscriminator("pfDeepCSVJetTags:probc") + thejet->bDiscriminator("pfDeepCSVJetTags:probb")+thejet->bDiscriminator("pfDeepCSVJetTags:probbb") ;
                      float ptD = thejet->userFloat("ptD") ;
                      float axis1 = thejet->userFloat("axis1") ;
                      int mult = thejet->userFloat("totalMult") ;

                      topTagger->addJet(thejet->pt(), thejet->eta(), thejet->phi(), thejet->mass(), bDiscriminatorValue, cvsl, cvsb, ptD, axis1, mult);           
                    }
                }
            }

            if(njet_ < jetsNumberThreshold_ || njets_btagmedium_ < bjetsNumberThreshold_) continue;
 
            if(debug_)
                cout << "Jets after selections " << njet_ << ", bJets " << njets_btagmedium_ << endl;

            std::sort(bTags.begin(),bTags.end(),std::greater<float>());

            // Set variables to compute MVA value

            leadeta_ = dipho->leadingPhoton()->eta();
            subleadeta_ = dipho->subLeadingPhoton()->eta();
            leadptom_ = dipho->leadingPhoton()->pt()/dipho->mass();
            subleadptom_ = dipho->subLeadingPhoton()->pt()/dipho->mass();
            leadIDMVA_ = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
            subleadIDMVA_ = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
            deltaphi_ = deltaPhi( dipho->leadingPhoton()->phi(), dipho->subLeadingPhoton()->phi() );
            leadPSV_ = dipho->leadingPhoton()->hasPixelSeed();
            subleadPSV_ = dipho->subLeadingPhoton()->hasPixelSeed();
            nJets_ = njet_;
            nJets_bTagMedium_ = njets_btagmedium_;

            if(tagJets.size()==0)
            {
                jet_pt1_ = -1;
                jet_pt2_ = -1.;
                jet_pt3_ = -1.;
                jet_eta1_ = -5;
                jet_eta2_ = -5.;
                jet_eta3_ = -5.;

                bTag1_ = -1;
                bTag2_ = -1;
            
            }else if(tagJets.size()==1)
            {
                jet_pt1_ = tagJets[0]->pt();
                jet_pt2_ = -1.;
                jet_pt3_ = -1.;
                jet_eta1_ = tagJets[0]->eta();
                jet_eta2_ = -5.;
                jet_eta3_ = -5.;

                bTag1_ = bTags[0];
                bTag2_ = -1.;
            }
            else if(tagJets.size()==2)
            {
                jet_pt1_ = tagJets[0]->pt();
                jet_pt2_ = tagJets[1]->pt();
                jet_pt3_ = -1.;
                jet_eta1_ = tagJets[0]->eta();
                jet_eta2_ = tagJets[1]->eta();
                jet_eta3_ = -5.;

                bTag1_ = bTags[0];
                bTag2_ = bTags[1];
            }
            else
            {
                jet_pt1_ = tagJets[0]->pt();
                jet_pt2_ = tagJets[1]->pt();
                jet_pt3_ = tagJets[2]->pt();
                jet_eta1_ = tagJets[0]->eta();
                jet_eta2_ = tagJets[1]->eta();
                jet_eta3_ = tagJets[2]->eta();

                bTag1_ = bTags[0];
                bTag2_ = bTags[1];
            }

            if( theMet_ -> size() != 1 )
                std::cout << "WARNING number of MET is not equal to 1" << std::endl;
            MetPt_ = theMet_->ptrAt( 0 ) -> getCorPt();
            MetPhi_ = theMet_->ptrAt( 0 ) -> phi();

            vector<float> mvaEval = topTagger->EvalMVA();
            cout << "eval top tagger: " << mvaEval.size() << endl;
            for (unsigned topt = 0; topt < mvaEval.size(); topt++)  {
              if (topt < mvaEval.size() - 1) cout << mvaEval[topt] << ", ";
              else cout << mvaEval[topt];
            }
            cout << endl;
            topTagger->clear(); 

            int leadMuIndex = 0;
            float leadMuPt = -1;
            int leadEleIndex = 0;
            float leadElePt = -1;

            for( unsigned int muonIndex = 0; muonIndex < Muons.size(); muonIndex++ )
            {
                Ptr<flashgg::Muon> muon = Muons[muonIndex];

                if(muon->pt()>leadMuPt)
                {
                    leadMuPt = muon->pt();
                    leadMuIndex = muonIndex;
                }
            }

            for( unsigned int eleIndex = 0; eleIndex < Electrons.size(); eleIndex++ )
            {
                Ptr<flashgg::Electron> ele = Electrons[eleIndex];

                if(ele->pt()>leadElePt)
                {
                    leadElePt = ele->pt();
                    leadEleIndex = eleIndex;
                }
            }

            if(leadMuPt>=leadElePt)
            {
                lepton_leadPt_ = Muons[leadMuIndex]->pt();
                lepton_leadEta_ = Muons[leadMuIndex]->eta();
                lepton_leadPhi_ = Muons[leadMuIndex]->phi();
            }
            else
            {
                lepton_leadPt_ = Electrons[leadEleIndex]->pt();
                lepton_leadEta_ = Electrons[leadEleIndex]->eta();
                lepton_leadPhi_ = Electrons[leadEleIndex]->phi();
            }

            float mvaValue = DiphotonMva_-> EvaluateMVA( "BDT" );
            int catNumber = -1;
            catNumber = chooseCategory( mvaValue , debug_);  

          if(debug_)
                cout << "I'm going to check selections, mva value: " << mvaValue << endl;

            if(SplitDiLeptEv_ && lepPt.size()>1 && njet_ >= DiLeptonJetThreshold_ && njets_btagmedium_ >= DiLeptonbJetThreshold_ && mvaValue > DiLeptonMVAThreshold_ ) 
            // Check DiLepton selection and assigne to purest cat if splitting is True
            {    catNumber = 0;
                 if(debug_)
                    cout << "DiLepton event with: " << njet_ << "jets, (threshold " << DiLeptonJetThreshold_ << ") " << njets_btagmedium_ << " bjets, (threshold " << DiLeptonbJetThreshold_  << ")" << mvaValue << " mva (threshold " << DiLeptonMVAThreshold_ << endl;

            }
            
            else if(lepPt.size()==1 && njet_ >= jetsNumberThreshold_ && njets_btagmedium_ >= bjetsNumberThreshold_) 
            // Check single lepton selections
            {    catNumber = chooseCategory( mvaValue, debug_ );  
               if(debug_)
                    cout << "Single lepton event with: "<< njet_ << " jets, (threshold " << DiLeptonJetThreshold_ << ") " << njets_btagmedium_ << " bjets, (threshold " << DiLeptonbJetThreshold_  << ")" << mvaValue << " mva (thresholds "  << endl;
            }

            if(debug_)
            { 
                cout << "TTHLeptonicTag -- MVA iput variables: " << endl;
                cout << "--------------------------------------" << endl;
                cout << "Lead and sublead photon eta " << leadeta_ << " " << subleadeta_ << endl;
                cout << "Lead and sublead photon pt/m " << leadptom_ << " " << subleadptom_ << endl;
                cout << "Lead and sublead photon IdMVA " << leadIDMVA_ << " " << subleadIDMVA_ << endl;
                cout << "Lead and sublead photon PSV " << leadPSV_ << " " << subleadPSV_ << endl;
                cout << "Photon delta phi " << deltaphi_ << endl;
                cout << "Number of jets " << nJets_ << endl;
                cout << "Number of b-jets " << nJets_bTagMedium_  << endl;
                cout << "Pt of the three leading jets " << jet_pt1_ << " " << jet_pt2_ << " " << jet_pt3_ << endl;
                cout << "Eta of the three leading jets " << jet_eta1_ << " " << jet_eta2_ << " " << jet_eta3_ << endl;
                cout << "Two highest bTag scores " << bTag1_ << " " << bTag2_ << endl;
                cout << "MetPt " << MetPt_ << endl;
                cout << "Lepton pT and Eta " << lepton_leadPt_ << " " << lepton_leadEta_ << endl;
                cout << "--------------------------------------" << endl;
                cout << "TTHLeptonicTag -- output MVA value " << mvaValue << " " << DiphotonMva_-> EvaluateMVA( "BDT" ) << ", category " << catNumber << endl;
            }

            if(catNumber!=-1)
            {
                TTHLeptonicTag tthltags_obj( dipho, mvares );
                tthltags_obj.setCategoryNumber(catNumber);

                for( unsigned int i = 0; i < tagJets.size(); ++i )
                {
                    tthltags_obj.includeWeightsByLabel( *tagJets[i] , "JetBTagReshapeWeight");
                }


                for( unsigned int i = 0; i < Muons.size(); ++i )
                    tthltags_obj.includeWeights( *Muons.at(i));

               for( unsigned int i = 0; i < Electrons.size(); ++i )
                    tthltags_obj.includeWeights( *Electrons.at(i));

                tthltags_obj.includeWeights( *dipho );
                tthltags_obj.setJets( tagJets );
                tthltags_obj.setMuons( Muons );
                tthltags_obj.setElectrons( Electrons );
                tthltags_obj.setDiPhotonIndex( diphoIndex );
                tthltags_obj.setSystLabel( systLabel_ );
                tthltags_obj.setMvaRes(mvaValue);
                tthltags_obj.setLepPt( lepPt );
                tthltags_obj.setLepE( lepE );
                tthltags_obj.setLepEta( lepEta );
                tthltags_obj.setLepPhi( lepPhi );
                tthltags_obj.setLepType( lepType );

                tthltags_obj.setLeadPrompt(-999);
                tthltags_obj.setLeadMad(-999);
                tthltags_obj.setLeadPythia(-999);
                tthltags_obj.setLeadPassFrix(-999);
                tthltags_obj.setLeadSimpleMomID(-999);
                tthltags_obj.setLeadSimpleMomStatus(-999);
                tthltags_obj.setLeadMomID(-999);
                tthltags_obj.setLeadMomMomID(-999);
                tthltags_obj.setLeadSmallestDr(-999);

                tthltags_obj.setSubleadPrompt(-999);
                tthltags_obj.setSubleadMad(-999);
                tthltags_obj.setSubleadPythia(-999);
                tthltags_obj.setSubleadPassFrix(-999);
                tthltags_obj.setSubleadSimpleMomID(-999);
                tthltags_obj.setSubleadSimpleMomStatus(-999);
                tthltags_obj.setSubleadMomID(-999);
                tthltags_obj.setSubleadMomMomID(-999);
                tthltags_obj.setSubleadSmallestDr(-999);

                // Custom stuff
                if( ! evt.isRealData() ) {
                    evt.getByToken( genParticleToken_, genParticles );
                    int nGoodEls(0), nGoodMus(0), nGoodElsFromTau(0), nGoodMusFromTau(0), nGoodTaus(0);
                    for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                        int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                        double pt = genParticles->ptrAt( genLoop )->p4().pt();
                        int status = genParticles->ptrAt( genLoop )->status(); 
                        bool isPromptFinalState = genParticles->ptrAt( genLoop )->isPromptFinalState();
                        bool isPromptDecayed = genParticles->ptrAt( genLoop )->isPromptDecayed();
                        bool isDirectPromptTauDecayProductFinalState = genParticles->ptrAt( genLoop )->isDirectPromptTauDecayProductFinalState();              
                        if (pt < 20) continue;
                        if (abs(pdgid) == 11 || abs(pdgid) == 13 || abs(pdgid) == 15) {
                            //cout << "Found a gen lepton/tau with pT > 20" << endl;
                            //cout << "pdgid: " << pdgid << endl;
                            //cout << "pt: " << pt << endl;
                            //cout << "status: " << status << endl;
                            //cout << "isPromptFinalState: " << isPromptFinalState << endl;
                            //cout << "isPromptDecayed: " << isPromptDecayed << endl;
                            //cout << "isDirectPromptTauDecayProductFinalState: " << isDirectPromptTauDecayProductFinalState << endl;
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
                    }
                    bool lead_photon_is_electron(false), sublead_photon_is_electron(false);
                    double lead_photon_eta = dipho->leadingPhoton()->eta();
                    double lead_photon_phi = dipho->leadingPhoton()->phi();
                    double sublead_photon_eta = dipho->subLeadingPhoton()->eta();
                    double sublead_photon_phi = dipho->subLeadingPhoton()->phi();
                    for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                        int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                        if (abs(pdgid) != 11) continue;
                        if (genParticles->ptrAt( genLoop )->p4().pt() < 15) continue;
                        if (!genParticles->ptrAt( genLoop )->isPromptFinalState()) continue;
                        double electron_eta = genParticles->ptrAt( genLoop )->p4().eta();
                        double electron_phi = genParticles->ptrAt( genLoop )->p4().phi();

                        double deltaR_lead = deltaR(electron_eta, electron_phi, lead_photon_eta, lead_photon_phi);
                        double deltaR_sublead = deltaR(electron_eta, electron_phi, sublead_photon_eta, sublead_photon_phi);

                        //cout << "Delta R with leading photon: " << deltaR_lead << endl;
                        //cout << "Delta R with subleading photon: " << deltaR_sublead << endl;

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
                    double lead_photon_closest_match = dipho->leadingPhoton()->genMatchType() == 1 ? deltaR(lead_photon_eta, lead_photon_phi, dipho->leadingPhoton()->matchedGenPhoton()->eta(), dipho->leadingPhoton()->matchedGenPhoton()->phi()) : 999;
                    double sublead_photon_closest_match = dipho->subLeadingPhoton()->genMatchType() == 1 ? deltaR(sublead_photon_eta, sublead_photon_phi, dipho->subLeadingPhoton()->matchedGenPhoton()->eta(), dipho->subLeadingPhoton()->matchedGenPhoton()->phi()) : 999; 

                    //double lead_photon_closest_match = dipho->leadingPhoton()->genMatchType() == 1 ? sqrt( pow(lead_photon_eta - dipho->leadingPhoton()->matchedGenPhoton()->eta(), 2) + pow(lead_photon_phi - dipho->leadingPhoton()->matchedGenPhoton()->phi(), 2)) : 999;
                    //double sublead_photon_closest_match = dipho->subLeadingPhoton()->genMatchType() == 1 ? sqrt( pow(sublead_photon_eta - dipho->subLeadingPhoton()->matchedGenPhoton()->eta(), 2) + pow(sublead_photon_phi - dipho->subLeadingPhoton()->matchedGenPhoton()->phi(), 2)) : 999;
                    double lead_photon_closest_match_pt = dipho->leadingPhoton()->genMatchType() == 1 ? dipho->leadingPhoton()->pt() : -1;
                    double sublead_photon_closest_match_pt = dipho->subLeadingPhoton()->genMatchType() == 1 ? dipho->subLeadingPhoton()->pt() : -1;

                    for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                        int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                        if (abs(pdgid) != 22) continue;
                        if (genParticles->ptrAt( genLoop )->p4().pt() < 10) continue;
                        if (!genParticles->ptrAt( genLoop )->isPromptFinalState()) continue;
                        double gen_photon_candidate_eta = genParticles->ptrAt( genLoop )->p4().eta();
                        double gen_photon_candidate_phi = genParticles->ptrAt( genLoop )->p4().phi();


                        double deltaR_lead = deltaR(gen_photon_candidate_eta, gen_photon_candidate_phi, lead_photon_eta, lead_photon_phi);
                        //double deltaR_lead = sqrt( pow(gen_photon_candidate_eta - lead_photon_eta, 2) + pow(gen_photon_candidate_phi - lead_photon_phi, 2));
                        double deltaR_sublead = deltaR(gen_photon_candidate_eta, gen_photon_candidate_phi, sublead_photon_eta, sublead_photon_phi); 
                        //double deltaR_sublead = sqrt( pow(gen_photon_candidate_eta - sublead_photon_eta, 2) + pow(gen_photon_candidate_phi - sublead_photon_phi, 2));

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


                    tthltags_obj.setnGoodEls(nGoodEls);
                    tthltags_obj.setnGoodElsFromTau(nGoodElsFromTau);
                    tthltags_obj.setnGoodMus(nGoodMus);
                    tthltags_obj.setnGoodMusFromTau(nGoodMusFromTau);
                    tthltags_obj.setnGoodTaus(nGoodTaus);
                    tthltags_obj.setLeadPhotonType(lead_photon_type);
                    tthltags_obj.setSubleadPhotonType(sublead_photon_type);
                    tthltags_obj.setLeadPhotonClosestDeltaR(lead_photon_closest_match);
                    tthltags_obj.setSubleadPhotonClosestDeltaR(sublead_photon_closest_match);
                    tthltags_obj.setLeadPhotonClosestPt(lead_photon_closest_match_pt);
                    tthltags_obj.setSubleadPhotonClosestPt(sublead_photon_closest_match_pt);
                }
                else {
                    tthltags_obj.setnGoodEls(-1);
                    tthltags_obj.setnGoodElsFromTau(-1);
                    tthltags_obj.setnGoodMus(-1);
                    tthltags_obj.setnGoodMusFromTau(-1);
                    tthltags_obj.setnGoodTaus(-1);
                    tthltags_obj.setLeadPhotonType(-1);
                    tthltags_obj.setSubleadPhotonType(-1);
                    tthltags_obj.setLeadPhotonClosestDeltaR(-1);
                    tthltags_obj.setSubleadPhotonClosestDeltaR(-1);
                    tthltags_obj.setLeadPhotonClosestPt(-1);
                    tthltags_obj.setSubleadPhotonClosestPt(-1);
                } 


                tthltags_obj.setRand(myRandLeptonic->Rndm());
                tthltags_obj.setMT( sqrt(2 * lepton_leadPt_ * MetPt_ * ( 1 - cos( lepton_leadPhi_ - MetPhi_ ) ) ) );
                tthltags_obj.setMetPt(MetPt_);
                tthltags_obj.setMetPhi(MetPhi_); 


                if (mvaEval.size() > 0) {
                  tthltags_obj.setTopTagScore(mvaEval[0] != -99 ? mvaEval[0] : -1);
                  tthltags_obj.setTopTagTopMass(mvaEval[4]);
                  tthltags_obj.setTopTagWMass(mvaEval[8]);
                } else
                  {
                  tthltags_obj.setTopTagScore(-999);
                  tthltags_obj.setTopTagTopMass(-999);
                  tthltags_obj.setTopTagWMass(-999);
                  }

                tthltags_obj.setNBLoose( njets_btagloose_ );
                tthltags_obj.setNBMedium( njets_btagmedium_ );
                tthltags_obj.setNBTight( njets_btagtight_ );

                tthltags_obj.setNMuonLoose( -1 );
                tthltags_obj.setNMuonMedium( Muons.size() );
                tthltags_obj.setNMuonTight( MuonsTight.size() );

                tthltags_obj.setNElecLoose( -1 );
                tthltags_obj.setNElecMedium( Electrons.size() );
                tthltags_obj.setNElecTight( ElectronsTight.size() );

                tthltags_obj.setMuonLeadIso( Muons.size() > 0 ? (Muons[0]->pfIsolationR04().sumChargedHadronPt + max(0., Muons[0]->pfIsolationR04().sumNeutralHadronEt + Muons[0]->pfIsolationR04().sumPhotonEt - 0.5*Muons[0]->pfIsolationR04().sumPUPt)) : -999);
                tthltags_obj.setMuonSubleadIso( Muons.size() > 1 ? (Muons[1]->pfIsolationR04().sumChargedHadronPt + max(0., Muons[1]->pfIsolationR04().sumNeutralHadronEt + Muons[1]->pfIsolationR04().sumPhotonEt - 0.5*Muons[1]->pfIsolationR04().sumPUPt)) : -999);

                tthltags->push_back( tthltags_obj );
 
                if( ! evt.isRealData() )
                {
                    TagTruthBase truth_obj;
                    truth_obj.setGenPV( higgsVtx );
                    if ( stage0cat.isValid() ) 
                    {   truth_obj.setHTXSInfo( *( stage0cat.product() ),
                                               *( stage1cat.product() ),
                                               *( njets.product() ),
                                               *( pTH.product() ),
                                               *( pTV.product() ) );
                    } else if ( htxsClassification.isValid() ) {
                        truth_obj.setHTXSInfo( htxsClassification->stage0_cat,
                                               htxsClassification->stage1_cat_pTjet30GeV,
                                               htxsClassification->jets30.size(),
                                               htxsClassification->p4decay_higgs.pt(),
                                               htxsClassification->p4decay_V.pt() );

                    } else {
                        truth_obj.setHTXSInfo( 0, 0, 0, 0., 0. );
                    }
                    truths->push_back( truth_obj );
                    tthltags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, idx++ ) ) );

                    int gp_lead_index = GenPhoIndex(genParticles, dipho->leadingPhoton(), -1);
                    int gp_sublead_index = GenPhoIndex(genParticles, dipho->subLeadingPhoton(), gp_lead_index);
                    vector<int> leadFlags; leadFlags.clear();
                    vector<int> subleadFlags; subleadFlags.clear();

                    if (gp_lead_index != -1) {
                       const edm::Ptr<reco::GenParticle> gp_lead = genParticles->ptrAt(gp_lead_index);                
                       leadFlags = IsPromptAfterOverlapRemove(genParticles, gp_lead);
                       tthltags->back().setLeadPrompt(leadFlags[0]);
                       tthltags->back().setLeadMad(leadFlags[1]);
                       tthltags->back().setLeadPythia(leadFlags[2]);
                       tthltags->back().setLeadPassFrix(leadFlags[3]);
                       tthltags->back().setLeadSimpleMomID(leadFlags[4]);
                       tthltags->back().setLeadSimpleMomStatus(leadFlags[5]);
                       tthltags->back().setLeadMomID(leadFlags[6]);
                       tthltags->back().setLeadMomMomID(leadFlags[7]);
                       tthltags->back().setLeadSmallestDr(NearestDr(genParticles, &(*gp_lead)));

                       cout << "leadPrompt: " << leadFlags[0] << endl;
                       } 

                   if (gp_sublead_index != -1) {
                       const edm::Ptr<reco::GenParticle> gp_sublead = genParticles->ptrAt(gp_sublead_index);
                       subleadFlags = IsPromptAfterOverlapRemove(genParticles, gp_sublead);
                       tthltags->back().setSubleadPrompt(subleadFlags[0]);
                       tthltags->back().setSubleadMad(subleadFlags[1]);
                       tthltags->back().setSubleadPythia(subleadFlags[2]);
                       tthltags->back().setSubleadPassFrix(subleadFlags[3]);
                       tthltags->back().setSubleadSimpleMomID(subleadFlags[4]);
                       tthltags->back().setSubleadSimpleMomStatus(subleadFlags[5]);
                       tthltags->back().setSubleadMomID(subleadFlags[6]);
                       tthltags->back().setSubleadMomMomID(subleadFlags[7]);
                       tthltags->back().setSubleadSmallestDr(NearestDr(genParticles, &(*gp_sublead)));
                       }

                }
            }
        }//diPho loop end !
        evt.put( std::move( tthltags ) );
        evt.put( std::move( truths ) );
    }

}
typedef flashgg::TTHLeptonicTagProducer FlashggTTHLeptonicTagProducer;
DEFINE_FWK_MODULE( FlashggTTHLeptonicTagProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

