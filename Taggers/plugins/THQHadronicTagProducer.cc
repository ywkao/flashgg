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
#include "flashgg/DataFormats/interface/THQHadronicTag.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/PatCandidates/interface/MET.h"

#include "flashgg/DataFormats/interface/Met.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"

#include "DataFormats/Math/interface/deltaR.h"

//include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "flashgg/Taggers/interface/SemiLepTopQuark.h"
//#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "flashgg/Taggers/interface/FoxWolfram.hpp"

#include "flashgg/DataFormats/interface/PDFWeightObject.h"
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
#include <vector>
#include <algorithm>
#include <string>
#include <utility>
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMVA/Reader.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "flashgg/Taggers/interface/TTH_DNN_Helper.h"
#include "flashgg/Taggers/interface/THQ_BDT_Helper.h"
#include "TCanvas.h"
#include <map>
#include <typeinfo>
using namespace std;
using namespace edm;


namespace flashgg {

class THQHadronicTagProducer : public EDProducer
{

public:
    typedef math::XYZPoint Point;
    THQHadronicTagProducer( const ParameterSet & );
    ~THQHadronicTagProducer();

private:
    std::string processId_;
//    edm::EDGetTokenT< LHEEventProduct > token_lhe;
//    int  chooseCategory( float, float);
    int  chooseCategory( float );
    void produce( Event &, const EventSetup & ) override;
    static bool sortByValue(const std::pair<int, double> &, const std::pair<int, double> &);
    const reco::GenParticle* getMother( const reco::GenParticle &part );
    std::vector<std::pair<int, double>> sortVector(const std::vector<double>);

    std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > tokenJets_;
    EDGetTokenT<View<DiPhotonCandidate> > diPhotonToken_;
    std::vector<edm::InputTag> inputTagJets_;
    EDGetTokenT<View<Electron> > electronToken_;
    EDGetTokenT<View<flashgg::Muon> > muonToken_;
    EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;
    EDGetTokenT<View<Photon> > photonToken_;
    EDGetTokenT<View<reco::Vertex> > vertexToken_;
    EDGetTokenT<View<flashgg::Met> > METToken_;
    EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
    EDGetTokenT<View<reco::GenJet> > genJetToken_;
    edm::EDGetTokenT<vector<flashgg::PDFWeightObject> > weightToken_;
    EDGetTokenT<double> rhoTag_;
    string systLabel_;
    FileInPath BDT_nrb_xml_file_;
    FileInPath BDT_smh_xml_file_;
    ///afs/cern.ch/work/y/ykao/tPrimeExcessHgg/CMSSW_10_6_8/src/tprimetH/mva/
    //mva_smh_m600->BookMVA("BDT", "./mva/Hadronic_tprime_SMH_varSet2_sigM600_bdt.xml");


    typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

    //Thresholds
    double muonPtThreshold_;
    double muonEtaThreshold_;
    vector<double> electronEtaThresholds_;
    double electronPtThreshold_;
    double leadPhoOverMassThreshold_;
    double subleadPhoOverMassThreshold_;
    double MVAThreshold_;
    double deltaRLepPhoThreshold_;
    double deltaRJetLepThreshold_;

    double deltaRJetLeadPhoThreshold_;
    double deltaRJetSubLeadPhoThreshold_;

    double jetsNumberThreshold_;
    double bjetsNumberThreshold_;
    double jetPtThreshold_;
    double jetEtaThreshold_;

    vector<double> bDiscriminator_;
    string bTag_;
    double muPFIsoSumRelThreshold_;
    double PhoMVAThreshold_;
    double DeltaRTrkElec_;

    double deltaRPhoElectronThreshold_;
    double deltaMassElectronZThreshold_;
    double DeltaRbjetfwdjet_;
    double DeltaRtHchainfwdjet_;

    bool hasGoodElec = false;
    bool hasVetoElec = false;
    bool hasGoodMuons = false;

    vector<double> boundaries;
    bool use_MVAs_;
    bool use_tthVstHDNN_;
    bool use_tthVstHBDT_;

    TLorentzVector G1, G2;  //temp solution: make met, bjet & jprime global TLorentzVectors

    //Variables for TMVA and Likelihood
    edm::Ptr<flashgg::Jet> fwdJet1;
    edm::Ptr<flashgg::Jet> bJet1;
    float lepton_ch_;
    TLorentzVector b1;
    float dipho_pt_ ;
    float dipho_leadPtOvermass_;
    float dipho_subleadPtOvermass_;
    float dipho_leadEta_;
    float dipho_subleadEta_;
    float dipho_leadIDMVA_;
    float dipho_subleadIDMVA_;
    float dipho_lead_haspixelseed_;
    float dipho_sublead_haspixelseed_;
    float n_jets_;
    float n_bjets_;
    float n_centraljets_;

    float dRbjetfwdjet_ ;
    float dRleadphobjet_ ;
    float dRsubleadphobjet_ ;
    float dRleadphofwdjet_ ;
    float dRsubleadphofwdjet_ ;
    float bjet1_discr_;
    float bjet2_discr_;
    float bjet3_discr_;
    float bjet1_pt_;
    float bjet2_pt_;
    float bjet3_pt_;
    float bjet1_eta_;
    float bjet2_eta_;
    float bjet3_eta_;
    float jet1_pt_;
    float jet2_pt_;
    float jet3_pt_;
    float jet4_pt_;
    float jet1_eta_;
    float jet2_eta_;
    float jet3_eta_;
    float jet4_eta_;
    float jet1_discr_;
    float jet2_discr_;
    float jet3_discr_;
    float jet4_discr_;
    float fwdJet1_pt_;
    float fwdJet1_eta_;
    float fwdJet1_discr_;

    float minPhoID_;
    float maxPhoID_;
    float maxBTagVal_;
    float secondMaxBTagVal_;
    float MVAscore_tHqVsttHDNN;

    THQ_BDT_Helper *tprimeTagger_nrb;
    THQ_BDT_Helper *tprimeTagger_smh;

    int counter = 0; //mytool
    bool debug  = false; //mytool
    bool debug_ = false;
    struct GreaterByPt
    {
    public:
        bool operator()( edm::Ptr<flashgg::Jet> lh, edm::Ptr<flashgg::Jet> rh ) const
        {
            return lh->pt() > rh->pt();
        };
    };

    struct GreaterByEta
    {
    public:
        bool operator()( edm::Ptr<flashgg::Jet> lh, edm::Ptr<flashgg::Jet> rh ) const
        {
            return fabs(lh->eta()) > fabs(rh->eta());
        };
    };

    struct GreaterByBTagging
    {
    public:
        GreaterByBTagging(std::string urName, std::string urName1):
            urName( urName ), urName1(urName1)
        {
        }

        bool operator()( edm::Ptr<flashgg::Jet> lh, edm::Ptr<flashgg::Jet> rh ) const
        {
            return (lh->bDiscriminator(urName.data()) + lh->bDiscriminator(urName1.data())) > (rh->bDiscriminator(urName.data()) + rh->bDiscriminator(urName1.data())) ;
        };
    private:
        const std::string urName, urName1;
    };


    //int LeptonType;
};

THQHadronicTagProducer::THQHadronicTagProducer( const ParameterSet &iConfig ) :
    processId_( iConfig.getParameter<string>("processId") ),
    diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
    electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
    muonToken_( consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
    mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag> ( "MVAResultTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
    METToken_( consumes<View<flashgg::Met> >( iConfig.getParameter<InputTag> ( "METTag" ) ) ),
    genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
    genJetToken_ ( consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag> ( "GenJetTag" ) ) ),
    weightToken_( consumes<vector<flashgg::PDFWeightObject> >( iConfig.getUntrackedParameter<InputTag>( "WeightTag", InputTag( "flashggPDFWeightObject" ) ) ) ),
    rhoTag_( consumes<double>( iConfig.getParameter<InputTag>( "rhoTag" ) ) ),
    systLabel_( iConfig.getParameter<string> ( "SystLabel" ) )
{

/*ps    if(processId_.find("thq") != std::string::npos or processId_.find("thw") != std::string::npos) {
        token_lhe = consumes<LHEEventProduct>( InputTag( "externalLHEProducer" )  );
    }
*/
    double default_deltaMassElectronZThreshold_ = 5.;
    vector<double> default_electronEtaCuts_;
    muonEtaThreshold_ = iConfig.getParameter<double>( "muonEtaThreshold" );
    muonPtThreshold_ = iConfig.getParameter<double>( "muonPtThreshold" );
    electronEtaThresholds_ = iConfig.getParameter<vector<double > >( "electronEtaThresholds" );
    electronPtThreshold_ = iConfig.getParameter<double>( "electronPtThreshold" );
    leadPhoOverMassThreshold_ = iConfig.getParameter<double>( "leadPhoOverMassThreshold" );
    subleadPhoOverMassThreshold_ = iConfig.getParameter<double>( "subleadPhoOverMassThreshold" );
    MVAThreshold_ = iConfig.getParameter<double>( "MVAThreshold" );
    deltaRLepPhoThreshold_ = iConfig.getParameter<double>( "deltaRLepPhoThreshold" );
    deltaRJetLepThreshold_ = iConfig.getParameter<double>( "deltaRJetLepThreshold" );
    jetsNumberThreshold_ = iConfig.getParameter<double>( "jetsNumberThreshold" );
    bjetsNumberThreshold_ = iConfig.getParameter<double>( "bjetsNumberThreshold" );
    jetPtThreshold_ = iConfig.getParameter<double>( "jetPtThreshold" );
    jetEtaThreshold_ = iConfig.getParameter<double>( "jetEtaThreshold" );

    deltaRJetLeadPhoThreshold_ = iConfig.getParameter<double>( "deltaRJetLeadPhoThreshold" );
    deltaRJetSubLeadPhoThreshold_ = iConfig.getParameter<double>( "deltaRJetSubLeadPhoThreshold" );

    bDiscriminator_ = iConfig.getParameter<vector<double > >( "bDiscriminator" );
    bTag_ = iConfig.getParameter<string>( "bTag" );

    muPFIsoSumRelThreshold_ = iConfig.getParameter<double>( "muPFIsoSumRelThreshold" );
    PhoMVAThreshold_ = iConfig.getParameter<double>( "PhoMVAThreshold" );
    DeltaRTrkElec_ = iConfig.getParameter<double>( "DeltaRTrkElec" );

    deltaRPhoElectronThreshold_ = iConfig.getParameter<double>( "deltaRPhoElectronThreshold" );
    deltaMassElectronZThreshold_ = iConfig.getUntrackedParameter<double>( "deltaMassElectronZThreshold_", default_deltaMassElectronZThreshold_ );
    DeltaRbjetfwdjet_ = iConfig.getParameter<double>( "DeltaRbjetfwdjet" );
    DeltaRtHchainfwdjet_ = iConfig.getParameter<double>( "DeltaRtHchainfwdjet" );
            
    boundaries = iConfig.getParameter<vector<double > >( "Boundaries" );
    use_MVAs_ = iConfig.getParameter<bool> ( "use_MVAs" );
    use_tthVstHDNN_ = iConfig.getParameter<bool> ( "use_tthVstHDNN" );
    use_tthVstHBDT_ = iConfig.getParameter<bool> ( "use_tthVstHBDT" );
    BDT_nrb_xml_file_ = iConfig.getParameter<edm::FileInPath> ( "tprime_bdt_nrb_xmlfile" );
    BDT_smh_xml_file_ = iConfig.getParameter<edm::FileInPath> ( "tprime_bdt_smh_xmlfile" );

    if(use_MVAs_) {
        tprimeTagger_nrb = new THQ_BDT_Helper(BDT_nrb_xml_file_.fullPath());
        tprimeTagger_smh = new THQ_BDT_Helper(BDT_smh_xml_file_.fullPath());
    }

    for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
        auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
        tokenJets_.push_back(token);
    }
    produces<vector<THQHadronicTag> >();
}

THQHadronicTagProducer::~THQHadronicTagProducer() {
}

int THQHadronicTagProducer::chooseCategory( float mvavalue )
{
    for(int n = 0 ; n < ( int )boundaries.size() ; n++ ) {
        if( ( double ) mvavalue > boundaries[boundaries.size() - n - 1] ) { return n; }
    }
    return -1;
}
/*int THQHadronicTagProducer::chooseCategory( float mvavalue1, float mvavalue2 )
    {
    int n;
        for( n = 0 ; n < ( int )boundaries.size() ; n++ ) {
            if( ( double )mvavalue1 > boundaries[boundaries.size() - n - 1] && ( double )mvavalue2 > boundaries[boundaries.size() - n - 1] ) { return n; }
        }
        return -1;
     }
*/

bool THQHadronicTagProducer::sortByValue(const std::pair<int, double> &pair1, const std::pair<int, double> &pair2)
{
    return fabs(pair1.second) < fabs(pair2.second);
}

std::vector<std::pair<int, double>> THQHadronicTagProducer::sortVector(const std::vector<double> v)
{
    std::vector<std::pair<int, double>> v2;
    for(unsigned int i=0; i<v.size(); ++i)
        v2.push_back( std::pair<int, double>(i, v[i]) );
    std::sort(v2.begin(), v2.end(), sortByValue);
    return v2;
}

const reco::GenParticle* THQHadronicTagProducer::getMother( const reco::GenParticle &part )
{
    const reco::GenParticle *mom = &part;

    while( mom->numberOfMothers() > 0 ) {
        for( unsigned int j=0;j<mom->numberOfMothers();++j ) {
            mom = dynamic_cast<const reco::GenParticle*>(mom->mother(j));
            if( mom->pdgId() != part.pdgId() ) return mom;
        }
    }
    
    return mom;
}

void THQHadronicTagProducer::produce( Event &evt, const EventSetup & )
{
    if(debug) {printf("[check] inside THQHadronicTagProducer::produce, processId = %s, counter = %d\n", processId_.c_str(), counter); counter += 1;}

    JetCollectionVector Jets( inputTagJets_.size() );
    for( size_t j = 0; j < inputTagJets_.size(); ++j ) evt.getByToken( tokenJets_[j], Jets[j] );

    std::unique_ptr<vector<THQHadronicTag> >  thqhtags( new vector<THQHadronicTag> );
    edm::Handle<double>                       rho          ; evt.getByToken( rhoTag_         , rho          );
    Handle<View<flashgg::DiPhotonCandidate> > diPhotons    ; evt.getByToken( diPhotonToken_  , diPhotons    );
    Handle<View<flashgg::Muon> >              theMuons     ; evt.getByToken( muonToken_      , theMuons     );
    Handle<View<flashgg::Electron> >          theElectrons ; evt.getByToken( electronToken_  , theElectrons );
    Handle<View<flashgg::DiPhotonMVAResult> > mvaResults   ; evt.getByToken( mvaResultToken_ , mvaResults   );
    Handle<View<reco::Vertex> >               vertices     ; evt.getByToken( vertexToken_    , vertices     );
    Handle<View<flashgg::Met> >               METs         ; evt.getByToken( METToken_       , METs         );
    Handle<View<reco::GenParticle> >          genParticles;
    Handle<View<reco::GenJet> >               genJets;

    assert( diPhotons->size() == mvaResults->size() );

    float rho_ = *rho;
    bool photonSelection = false;
    double idmva1 = 0.;
    double idmva2 = 0.;
    for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ ) {
        if(debug) printf("[check] diphoIndex = %d/%d\n", diphoIndex, diPhotons->size());

        hasGoodElec  = false;
        hasVetoElec  = false;
        hasGoodMuons = false;

        unsigned int jetCollectionIndex = diPhotons->ptrAt( diphoIndex )->jetCollectionIndex();

        edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );
        edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );

        flashgg::THQHadronicTag thqhtags_obj( dipho, mvares );

        if( dipho->leadingPhoton()->pt() < ( dipho->mass() )*leadPhoOverMassThreshold_ ) continue;
        if( dipho->subLeadingPhoton()->pt() < ( dipho->mass() )*subleadPhoOverMassThreshold_ ) continue;
        idmva1 = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
        idmva2 = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
        if(debug) printf("[check] idmva1 = %.3f, idmva2 = %.3f\n", idmva1, idmva2);
        if( idmva1 < PhoMVAThreshold_ || idmva2 < PhoMVAThreshold_ ) continue;
        //if( (! evt.isRealData()) && (idmva1 < PhoMVAThreshold_ || idmva2 < PhoMVAThreshold_) ) continue; // for ntuple production, data skips the condition
        if(debug) printf("[check] survive idmva cuts!\n");

        //if( mvares->result < MVAThreshold_ ) continue;            //DiPho_MVA

        photonSelection = true;

        G1.SetPtEtaPhiM( diPhotons->ptrAt( diphoIndex )->leadingPhoton()->pt(),
                         diPhotons->ptrAt( diphoIndex )->leadingPhoton()->eta(),
                         diPhotons->ptrAt( diphoIndex )->leadingPhoton()->phi() ,
                         0 );

        G2.SetPtEtaPhiM( diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->pt(),
                         diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->eta(),
                         diPhotons->ptrAt( diphoIndex )->subLeadingPhoton()->phi(),
                         0 );

        if( METs->size() != 1 ) std::cout << "WARNING - #MET is not 1" << std::endl;
        Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );
        thqhtags_obj.setRECOMET(theMET);

        std::vector<edm::Ptr<flashgg::Muon> > LooseMu15  = selectMuons( theMuons->ptrs(), dipho, vertices->ptrs(), muonEtaThreshold_, muonPtThreshold_, 0.15, deltaRLepPhoThreshold_, deltaRLepPhoThreshold_ );
        std::vector<edm::Ptr<flashgg::Muon> > LooseMu25  = selectMuons( theMuons->ptrs(), dipho, vertices->ptrs(), muonEtaThreshold_, muonPtThreshold_, 0.25, deltaRLepPhoThreshold_, deltaRLepPhoThreshold_ );
        std::vector<edm::Ptr<flashgg::Muon> > LooseMu200 = selectMuons( theMuons->ptrs(), dipho, vertices->ptrs(), muonEtaThreshold_, muonPtThreshold_, 2.  , deltaRLepPhoThreshold_, deltaRLepPhoThreshold_ );
        std::vector<edm::Ptr<flashgg::Muon> > MediumMu15 = selectMuons( theMuons->ptrs(), dipho, vertices->ptrs(), muonEtaThreshold_, muonPtThreshold_, 0.15, deltaRLepPhoThreshold_, deltaRLepPhoThreshold_ );
        std::vector<edm::Ptr<flashgg::Muon> > MediumMu25 = selectMuons( theMuons->ptrs(), dipho, vertices->ptrs(), muonEtaThreshold_, muonPtThreshold_, 0.25, deltaRLepPhoThreshold_, deltaRLepPhoThreshold_ );
        std::vector<edm::Ptr<flashgg::Muon> > TightMuo15 = selectMuons( theMuons->ptrs(), dipho, vertices->ptrs(), muonEtaThreshold_, muonPtThreshold_, 0.15, deltaRLepPhoThreshold_, deltaRLepPhoThreshold_ );
        std::vector<edm::Ptr<flashgg::Muon> > TightMuo25 = selectMuons( theMuons->ptrs(), dipho, vertices->ptrs(), muonEtaThreshold_, muonPtThreshold_, 0.25, deltaRLepPhoThreshold_, deltaRLepPhoThreshold_ );
        std::vector<edm::Ptr<flashgg::Muon> > goodMuons  = muPFIsoSumRelThreshold_== 0.15 ? TightMuo15 : TightMuo25 ;

        std::vector<int> looseMus_PassTight;
        for(auto mu: LooseMu200) looseMus_PassTight.push_back( std::find( goodMuons.begin() , goodMuons.end() , mu ) != goodMuons.end() );

        std::vector<edm::Ptr<Electron> > vetoNonIsoElectrons = selectStdElectrons(theElectrons->ptrs(), dipho, vertices->ptrs(), electronPtThreshold_, electronEtaThresholds_, 0, 4, deltaRPhoElectronThreshold_, DeltaRTrkElec_, deltaMassElectronZThreshold_, rho_, true );
        std::vector<edm::Ptr<Electron> > looseElectrons      = selectStdElectrons(theElectrons->ptrs(), dipho, vertices->ptrs(), electronPtThreshold_, electronEtaThresholds_, 0, 3, deltaRPhoElectronThreshold_, DeltaRTrkElec_, deltaMassElectronZThreshold_, rho_, true );
        std::vector<edm::Ptr<Electron> > vetoElectrons       = selectStdElectrons(theElectrons->ptrs(), dipho, vertices->ptrs(), electronPtThreshold_, electronEtaThresholds_, 0, 0, deltaRPhoElectronThreshold_, DeltaRTrkElec_, deltaMassElectronZThreshold_, rho_, true );
        std::vector<edm::Ptr<Electron> > mediumElectrons     = selectStdElectrons(theElectrons->ptrs(), dipho, vertices->ptrs(), electronPtThreshold_, electronEtaThresholds_, 0, 2, deltaRPhoElectronThreshold_, DeltaRTrkElec_, deltaMassElectronZThreshold_, rho_, true );
        std::vector<edm::Ptr<Electron> > goodElectrons       = selectStdElectrons(theElectrons->ptrs(), dipho, vertices->ptrs(), electronPtThreshold_, electronEtaThresholds_, 0, 1, deltaRPhoElectronThreshold_, DeltaRTrkElec_, deltaMassElectronZThreshold_, rho_, true );

        std::vector<int> vetoNonIsoElectrons_PassTight;
        std::vector<int> vetoNonIsoElectrons_PassVeto;
        for(auto ele : vetoNonIsoElectrons ) {
            vetoNonIsoElectrons_PassTight.push_back( std::find( goodElectrons.begin() , goodElectrons.end() , ele ) != goodElectrons.end() );
            vetoNonIsoElectrons_PassVeto.push_back( std::find( vetoElectrons.begin() , vetoElectrons.end() , ele ) != vetoElectrons.end() );
        }

        if((goodElectrons.size() + goodMuons.size()) > 0) continue; // selection of hadronic channel
        if(debug) printf("[check] survive lepton veto!\n");

        //hasGoodElec = ( goodElectrons.size() == 1 );
        //hasVetoElec = ( vetoElectrons.size() > 0 );
        //hasGoodMuons = ( goodMuons.size() == 1 );

        //LeptonType = 0; //1 : electron, 2:muon
        //if( hasGoodMuons && !hasVetoElec) LeptonType = 2;

        float ht = 0;
        float dRPhoLeadJet    = 0;
        float dRPhoSubLeadJet = 0;
        double minDrLepton    = 999.;
        int njets_btagloose_  = 0;
        int njets_btagmedium_ = 0;
        int njets_btagtight_  = 0;
        std::vector<edm::Ptr<flashgg::Jet> > SelJetVect;
        std::vector<edm::Ptr<flashgg::Jet> > SelJetVect_EtaSorted;
        std::vector<edm::Ptr<flashgg::Jet> > SelJetVect_PtSorted;
        std::vector<edm::Ptr<flashgg::Jet> > SelJetVect_BSorted;
        std::vector<edm::Ptr<flashgg::Jet> > MediumBJetVect, MediumBJetVect_PtSorted, BJetVect;
        std::vector<edm::Ptr<flashgg::Jet> > LooseBJetVect, LooseBJetVect_PtSorted;
        std::vector<edm::Ptr<flashgg::Jet> > TightBJetVect, TightBJetVect_PtSorted;
        std::vector<edm::Ptr<flashgg::Jet> > centraljet;
        std::vector<edm::Ptr<flashgg::Jet> > forwardjet;
	    std::vector<float> bDiscr_bjets;
        std::vector<float> bDiscr_jets;
        std::vector<float> bDiscr_fwdjets;
        for( unsigned int candIndex_outer = 0; candIndex_outer < Jets[jetCollectionIndex]->size() ; candIndex_outer++ ) {
            edm::Ptr<flashgg::Jet> thejet = Jets[jetCollectionIndex]->ptrAt( candIndex_outer );
            //std::vector <float> minDrLepton_ele;
            //std::vector <float> minDrLepton_muon;

            if( !thejet->passesPuJetId( dipho ) ) continue;
            if( fabs( thejet->eta() ) > jetEtaThreshold_ ) continue;
            if( thejet->pt() < jetPtThreshold_ ) continue;

            dRPhoLeadJet    = deltaR( thejet->eta(), thejet->phi(), dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ;
            dRPhoSubLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->subLeadingPhoton()->superCluster()->eta(),dipho->subLeadingPhoton()->superCluster()->phi() );

            if( dRPhoLeadJet < deltaRJetLeadPhoThreshold_ || dRPhoSubLeadJet < deltaRJetSubLeadPhoThreshold_ ) continue;

            float bDiscriminatorValue; //= -2.;
            if(bTag_ == "pfDeepCSV") bDiscriminatorValue = thejet->bDiscriminator("pfDeepCSVJetTags:probb") + thejet->bDiscriminator("pfDeepCSVJetTags:probbb");
            else bDiscriminatorValue = thejet->bDiscriminator( bTag_ );

            if( bDiscriminatorValue > bDiscriminator_[0] ) {
                LooseBJetVect_PtSorted.push_back( thejet );
                LooseBJetVect.push_back( thejet );
                njets_btagloose_++;
            }
            if( bDiscriminatorValue > bDiscriminator_[1] ) {
                MediumBJetVect.push_back( thejet );
                MediumBJetVect_PtSorted.push_back( thejet );
                njets_btagmedium_++;
            }
            if( bDiscriminatorValue > bDiscriminator_[2] ) {
                TightBJetVect_PtSorted.push_back( thejet );
                TightBJetVect.push_back( thejet );
                njets_btagtight_++;
            }

            ht+=thejet->pt();

            if( fabs( thejet->eta() ) < 1) centraljet.push_back( thejet );
            if( fabs( thejet->eta() ) > 1) forwardjet.push_back( thejet );

            SelJetVect.push_back( thejet );
            SelJetVect_EtaSorted.push_back( thejet );
            SelJetVect_PtSorted.push_back( thejet );
            bDiscr_jets.push_back( thejet->bDiscriminator("pfDeepCSVJetTags:probb") + thejet->bDiscriminator("pfDeepCSVJetTags:probbb") );

            if(use_MVAs_)
            {
                tprimeTagger_nrb->addJet(thejet->pt(), thejet->eta(), thejet->phi(), thejet->mass(), bDiscriminatorValue);
                tprimeTagger_smh->addJet(thejet->pt(), thejet->eta(), thejet->phi(), thejet->mass(), bDiscriminatorValue);
            }
        }//end of jets loop


        if(debug) printf("[debug] num leps = %lu, num_jets = %lu \n", goodElectrons.size() + goodMuons.size(), SelJetVect.size());

        //Calculate scalar sum of jets
        std::sort( LooseBJetVect_PtSorted.begin()  , LooseBJetVect_PtSorted.end()  , GreaterByPt()  );
        std::sort( MediumBJetVect_PtSorted.begin() , MediumBJetVect_PtSorted.end() , GreaterByPt()  );
        std::sort( TightBJetVect_PtSorted.begin()  , TightBJetVect_PtSorted.end()  , GreaterByPt()  );
        std::sort( SelJetVect_EtaSorted.begin()    , SelJetVect_EtaSorted.end()    , GreaterByEta() );
        std::sort( SelJetVect_PtSorted.begin()     , SelJetVect_PtSorted.end()     , GreaterByPt()  );

	    for(unsigned int bjetsindex = 0 ; bjetsindex < LooseBJetVect.size(); bjetsindex++){
	        if(LooseBJetVect[bjetsindex] !=  SelJetVect_EtaSorted[0] ){
	            BJetVect.push_back( LooseBJetVect[bjetsindex] );
	            bDiscr_bjets.push_back( LooseBJetVect[bjetsindex]->bDiscriminator("pfDeepCSVJetTags:probb") + LooseBJetVect[bjetsindex]->bDiscriminator("pfDeepCSVJetTags:probbb") );
	        }
	    }
	    LooseBJetVect.clear();
	    std::sort(BJetVect.begin(),BJetVect.end(), GreaterByBTagging("pfDeepCSVJetTags:probb", "pfDeepCSVJetTags:probbb"));
	    std::sort(bDiscr_bjets.begin(), bDiscr_bjets.end(), std::greater<float>());

        if(SelJetVect.size() < jetsNumberThreshold_ || BJetVect.size() < bjetsNumberThreshold_) continue;	
        if(debug) printf("[check] survive jets criteria!\n");

        bjet1_discr_   = -999;
        bjet2_discr_   = -999;
        bjet3_discr_   = -999;
        bjet1_pt_      = -999;
        bjet2_pt_      = -999;
        bjet3_pt_      = -999;
        bjet1_eta_     = -999;
        bjet2_eta_     = -999;
        bjet3_eta_     = -999;
        jet1_pt_       = -999;
        jet2_pt_       = -999;
        jet3_pt_       = -999;
        jet4_pt_       = -999;
        jet1_eta_      = -999;
        jet2_eta_      = -999;
        jet3_eta_      = -999;
        jet4_eta_      = -999;
        jet1_discr_    = -999;
        jet2_discr_    = -999;
        jet3_discr_    = -999;
        jet4_discr_    = -999;
        fwdJet1_pt_    = -999;
        fwdJet1_eta_   = -999;
        fwdJet1_discr_ = -999;

//------------------------------------ Likelihood and MVA -----------------------------------------//

        fwdJet1 = SelJetVect_EtaSorted[0];
        bJet1 = BJetVect[0];
        b1.SetPtEtaPhiE(0., 0., 0., 0.);
        b1.SetPtEtaPhiE(bJet1->pt(), bJet1->eta(), bJet1->phi(), bJet1->energy());

        dipho_pt_ = dipho->pt(); 
        dipho_leadPtOvermass_ = dipho->leadingPhoton()->pt()/dipho->mass();
        dipho_subleadPtOvermass_ = dipho->subLeadingPhoton()->pt()/dipho->mass();
        dipho_leadEta_ = dipho->leadingPhoton()->eta();
        dipho_subleadEta_ = dipho->subLeadingPhoton()->eta();
        dipho_leadIDMVA_ = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
        dipho_subleadIDMVA_ = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
        dipho_lead_haspixelseed_ = dipho->leadingPhoton()->hasPixelSeed();
        dipho_sublead_haspixelseed_ = dipho->subLeadingPhoton()->hasPixelSeed(); 
        n_jets_ = SelJetVect.size();
        n_bjets_ = BJetVect.size();
        n_centraljets_ = centraljet.size();	
        bjet1_pt_ = bJet1->pt();
        bjet1_eta_ = bJet1->eta();
        bjet1_discr_ = bDiscr_bjets.at(0);
        
        if(SelJetVect.size()>0){
            jet1_pt_ = SelJetVect[0]->pt();
            jet1_eta_ = SelJetVect[0]->eta();
            jet1_discr_= SelJetVect[0]->bDiscriminator("pfDeepCSVJetTags:probb") + SelJetVect[0]->bDiscriminator("pfDeepCSVJetTags:probbb");
        }
        if(SelJetVect.size()>1){
            jet2_pt_ = SelJetVect[1]->pt();
            jet2_eta_ = SelJetVect[1]->eta();
            jet2_discr_= SelJetVect[1]->bDiscriminator("pfDeepCSVJetTags:probb") + SelJetVect[1]->bDiscriminator("pfDeepCSVJetTags:probbb");
        }
        if(SelJetVect.size()>2){
            jet3_pt_ = SelJetVect[2]->pt();
            jet3_eta_ = SelJetVect[2]->eta();
            jet3_discr_= SelJetVect[2]->bDiscriminator("pfDeepCSVJetTags:probb") + SelJetVect[2]->bDiscriminator("pfDeepCSVJetTags:probbb");
        }
        if(SelJetVect.size()>3){
            jet4_pt_ = SelJetVect[3]->pt();
            jet4_eta_ = SelJetVect[3]->eta();
            jet4_discr_= SelJetVect[3]->bDiscriminator("pfDeepCSVJetTags:probb") + SelJetVect[3]->bDiscriminator("pfDeepCSVJetTags:probbb");
        }
        if(BJetVect.size()>0){
            bjet1_pt_ = BJetVect[0]->pt();
            bjet1_eta_ = BJetVect[0]->eta();
            bjet1_discr_= BJetVect[0]->bDiscriminator("pfDeepCSVJetTags:probb") + BJetVect[0]->bDiscriminator("pfDeepCSVJetTags:probbb");
        }
        if(BJetVect.size()>1){
            bjet2_pt_ = BJetVect[1]->pt();
            bjet2_eta_ = BJetVect[1]->eta();
            bjet2_discr_= BJetVect[1]->bDiscriminator("pfDeepCSVJetTags:probb") + BJetVect[1]->bDiscriminator("pfDeepCSVJetTags:probbb");
        }
        if(BJetVect.size()>2){
            bjet3_pt_ = BJetVect[2]->pt();
            bjet3_eta_ = BJetVect[2]->eta();
            bjet3_discr_= BJetVect[2]->bDiscriminator("pfDeepCSVJetTags:probb") + BJetVect[2]->bDiscriminator("pfDeepCSVJetTags:probbb");
        }
        if(SelJetVect_EtaSorted.size()>0){
            fwdJet1_pt_ = SelJetVect_EtaSorted[0]->pt();
            fwdJet1_eta_ = SelJetVect_EtaSorted[0]->eta();
            fwdJet1_discr_= SelJetVect_EtaSorted[0]->bDiscriminator("pfDeepCSVJetTags:probb") + SelJetVect_EtaSorted[0]->bDiscriminator("pfDeepCSVJetTags:probbb");
        }
        bDiscr_fwdjets.push_back(fwdJet1_discr_ );

        dRbjetfwdjet_       = deltaR( b1.Eta() , b1.Phi() , fwdJet1->eta() , fwdJet1->phi() );
        dRleadphobjet_      = deltaR( G1.Eta() , G1.Phi(), b1.Eta() , b1.Phi());
        dRsubleadphobjet_   = deltaR( G2.Eta() , G2.Phi(), b1.Eta() , b1.Phi());
        dRleadphofwdjet_    = deltaR( G1.Eta() , G1.Phi(), fwdJet1->eta() , fwdJet1->phi());
        dRsubleadphofwdjet_ = deltaR( G2.Eta() , G2.Phi(), fwdJet1->eta() , fwdJet1->phi());

        minPhoID_           = TMath::Min( dipho_leadIDMVA_, dipho_subleadIDMVA_);
        maxPhoID_           = TMath::Max( dipho_leadIDMVA_, dipho_subleadIDMVA_);

        maxBTagVal_         = bDiscr_bjets.size() > 0 ? bDiscr_bjets[0] : -1.;
        secondMaxBTagVal_   = bDiscr_bjets.size() > 1 ? bDiscr_bjets[1]: -1.;

        double mva_value_nrb = -999;
        double mva_value_smh = -999;
        if(use_MVAs_) {
            tprimeTagger_nrb->addPhoton(dipho->leadingPhoton()->pt(), dipho->leadingPhoton()->eta(), dipho->leadingPhoton()->phi(), 0., dipho_leadIDMVA_);
            tprimeTagger_nrb->addPhoton(dipho->subLeadingPhoton()->pt(), dipho->subLeadingPhoton()->eta(), dipho->subLeadingPhoton()->phi(), 0., dipho_subleadIDMVA_);
            tprimeTagger_nrb->addPhotonPixelSeed(dipho_lead_haspixelseed_, dipho_sublead_haspixelseed_);
            tprimeTagger_nrb->addNbjets((float) MediumBJetVect.size());
            tprimeTagger_nrb->addHt(ht);
            tprimeTagger_nrb->addMet(theMET->getCorPt());
            mva_value_nrb = tprimeTagger_nrb->EvalMVA();
            tprimeTagger_nrb->clear();

            tprimeTagger_smh->addPhoton(dipho->leadingPhoton()->pt(), dipho->leadingPhoton()->eta(), dipho->leadingPhoton()->phi(), 0., dipho_leadIDMVA_);
            tprimeTagger_smh->addPhoton(dipho->subLeadingPhoton()->pt(), dipho->subLeadingPhoton()->eta(), dipho->subLeadingPhoton()->phi(), 0., dipho_subleadIDMVA_);
            tprimeTagger_smh->addPhotonPixelSeed(dipho_lead_haspixelseed_, dipho_sublead_haspixelseed_);
            tprimeTagger_smh->addNbjets((float) MediumBJetVect.size());
            tprimeTagger_smh->addHt(ht);
            tprimeTagger_smh->addMet(theMET->getCorPt());
            mva_value_smh = tprimeTagger_smh->EvalMVA();
            tprimeTagger_smh->clear();
        }

        // select events that pass specified bdt scores
        if(mva_value_nrb < 0.56) continue;
        if(mva_value_smh < 0.52) continue;

        /*---------------------------------------------------------------------------------------
        # Evaluate MVA
        # DNN input variables
        # Evaluate DNN
        int catnum = -1;
        catnum = chooseCategory( MVAscore_tHqVsttHBDT );
        catnum = chooseCategory( mvares->result );
        catnum = chooseCategory( idmva1, idmva2 );
        ---------------------------------------------------------------------------------------*/
        if(debug_){
            cout<<"jet1_eta_"<<jet1_eta_<<endl;
            cout<<"jet2_eta_"<<jet2_eta_<<endl;
            cout<<"jet3_eta_"<<jet3_eta_<<endl;
            cout<<"jet4_eta_"<<jet4_eta_<<endl;
            cout<<"jet1_pt_"<<jet1_pt_<<endl;
            cout<<"jet2_pt_"<<jet2_pt_<<endl;
            cout<<"jet3_pt_"<<jet3_pt_<<endl;
            cout<<"jet4_pt_"<<jet4_pt_<<endl;
            cout<<"jet1_discr_"<<jet1_discr_<<endl;
            cout<<"jet2_discr_"<<jet2_discr_<<endl;
            cout<<"jet3_discr_"<<jet3_discr_<<endl;
            cout<<"jet4_discr_"<<jet4_discr_<<endl;
            cout<<"bjet1_pt_"<<bjet1_pt_<<endl;
            cout<<"bjet2_pt_"<<bjet2_pt_<<endl;
            cout<<"bjet3_pt_"<<bjet3_pt_<<endl;
            cout<<"bjet1_eta_"<<bjet1_eta_<<endl;
            cout<<"bjet2_eta_"<<bjet2_eta_<<endl;
            cout<<"bjet3_eta_"<<bjet3_eta_<<endl;
            cout<<"bjet1_discr_"<<bjet1_discr_<<endl;
            cout<<"bjet2_discr_"<<bjet2_discr_<<endl;
            cout<<"bjet3_discr_"<<bjet3_discr_<<endl;
            cout<<"fwdJet1_pt_= "<<fwdJet1_pt_<<endl;
            cout<<"fwdJet1_pt_= "<<fwdJet1_eta_<<endl;
            cout<<"fwdJet1_discr_= "<<fwdJet1_discr_<<endl;
        }

//-------------------------------------------------------------------------------------------------//

        thqhtags_obj.setHT(ht);

        if( photonSelection ) {
            thqhtags_obj.setMVAscore_nrb(mva_value_nrb);
            thqhtags_obj.setMVAscore_smh(mva_value_smh);
            thqhtags_obj.setrho(rho_);

            //thqhtags_obj.setLeptonType(LeptonType);
            thqhtags_obj.includeWeights( *dipho );
            thqhtags_obj.photonWeights = dipho->leadingPhoton()->centralWeight()*dipho->subLeadingPhoton()->centralWeight() ;
            thqhtags_obj.setJets( SelJetVect_PtSorted , SelJetVect_EtaSorted);
            thqhtags_obj.setBJets( BJetVect );
            thqhtags_obj.nCentralJets = centraljet.size();
            thqhtags_obj.nForwardJets = forwardjet.size();
            thqhtags_obj.setcentraljet( centraljet );
            thqhtags_obj.setforwardjet( forwardjet );
            thqhtags_obj.setdRbjetfwdjet( dRbjetfwdjet_ ) ;
            thqhtags_obj.setdRleadphobjet( dRleadphobjet_ );
            thqhtags_obj.setdRsubleadphobjet( dRsubleadphobjet_ );
            thqhtags_obj.setdRleadphofwdjet( dRleadphofwdjet_ );
            thqhtags_obj.setdRsubleadphofwdjet( dRsubleadphofwdjet_ );
            thqhtags_obj.setmvaresult ( mvares->result ) ;     //diphoton mva
		    thqhtags_obj.setbDiscriminatorValue( bDiscr_bjets, bDiscr_jets, bDiscr_fwdjets );
            //thqhtags_obj.setCategoryNumber( catnum );
            thqhtags_obj.bTagWeight     = 1.0;
            thqhtags_obj.bTagWeightDown = 1.0;
            thqhtags_obj.bTagWeightUp   = 1.0;
            for( auto j : SelJetVect_PtSorted ) {
                thqhtags_obj.includeWeights( *j );

                //if( j->hasWeight("JetBTagCutWeightCentral") ){
                //    thqhtags_obj.bTagWeight *= j->weight( "JetBTagCutWeightCentral" );
                //    thqhtags_obj.bTagWeightDown *= j->weight( "JetBTagCutWeightDown01sigma" );
                //    thqhtags_obj.bTagWeightUp *= j->weight( "JetBTagCutWeightUp01sigma" );
                //}
                //else
                //  cout << "BTag weight is not set in jet" << endl;
            }

            thqhtags_obj.setVertices( vertices->ptrs() );
            thqhtags_obj.setMuons( goodMuons , looseMus_PassTight , LooseMu25.size() , LooseMu15.size() , MediumMu25.size() , MediumMu15.size() , TightMuo25.size() , TightMuo15.size() );
            thqhtags_obj.setElectrons( goodElectrons , vetoNonIsoElectrons_PassTight , vetoNonIsoElectrons_PassVeto , looseElectrons.size() , vetoElectrons.size() , mediumElectrons.size() , goodElectrons.size() );
            thqhtags_obj.setDiPhotonIndex( diphoIndex );
            thqhtags_obj.setSystLabel( systLabel_ );
            thqhtags_obj.nMedium_bJets = MediumBJetVect.size();
            thqhtags_obj.nLoose_bJets  = LooseBJetVect_PtSorted.size();
            thqhtags_obj.nTight_bJets  = TightBJetVect_PtSorted.size();

//------------------------------------ MC-Truth Info -----------------------------------------//
//            if ( ! evt.isRealData() ) {
//                evt.getByToken( genParticleToken_, genParticles );
//                evt.getByToken( genJetToken_, genJets );
//
//                //--- gen met ---//
//                TLorentzVector nu_lorentzVector, allnus_LorentzVector, promptnus_LorentzVector;
//
//                for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
//                    edm::Ptr<reco::GenParticle> part = genParticles->ptrAt( genLoop );
//                    bool fid_cut = (abs(part->eta())<5.0 && part->status()==1) ? 1 : 0;
//                    bool isNu = (abs(part->pdgId())==12 || abs(part->pdgId())==14 || abs(part->pdgId())==16) ? 1 : 0;
//                    if (!fid_cut || !isNu) continue;
//                    if( part->isPromptFinalState() || part->isDirectPromptTauDecayProductFinalState()) {
//                        nu_lorentzVector.SetPtEtaPhiE(  part->pt() , part->eta() , part->phi() , part->energy() );
//                        promptnus_LorentzVector+=nu_lorentzVector;
//                    }
//                    else {
//                        nu_lorentzVector.SetPtEtaPhiE(  part->pt() , part->eta() , part->phi() , part->energy() );
//                        allnus_LorentzVector+=nu_lorentzVector;
//                    }
//                }
//                thqhtags_obj.setMETPtEtaPhiE( "allPromptNus", promptnus_LorentzVector.Pt(), promptnus_LorentzVector.Eta(), promptnus_LorentzVector.Phi(), promptnus_LorentzVector.Energy() );
//                thqhtags_obj.setMETPtEtaPhiE( "allNus", allnus_LorentzVector.Pt(), allnus_LorentzVector.Eta(), allnus_LorentzVector.Phi(), allnus_LorentzVector.Energy() );
//                thqhtags_obj.setMETPtEtaPhiE( "genMetTrue", theMET->genMET()->pt(), theMET->genMET()->eta(), theMET->genMET()->phi(), theMET->genMET()->energy() );
//                thqhtags_obj.setStage1recoTag( DiPhotonTagBase::stage1recoTag::RECO_THQ_LEP );
//
//                // MC truth matching
//                std::vector<int> vec_pdgId_register;
//                std::vector<int> vec_index_register;
//                std::vector<double> vec_deltaR_register;
//                std::vector<int> gens_pdgId;
//                std::vector<int> gens_status;
//                std::vector<double> gens_pt;
//                std::vector<double> gens_eta;
//                std::vector<double> gens_phi;
//                std::vector<double> gens_mass;
//                std::vector<int> moms_pdgId;
//                std::vector<int> moms_status;
//                std::vector<double> moms_pt;
//                std::vector<double> moms_eta;
//                std::vector<double> moms_phi;
//                std::vector<double> moms_mass;
//                int counter_num_gen_particles = 0;
//                for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
//                    edm::Ptr<reco::GenParticle> part = genParticles->ptrAt( genLoop );
//                    //printf("[check] ");
//                    //printf("(%2d) ", genLoop);
//                    //printf("status = %5d, ", part->status());
//                    //printf("pdgId = %5d, ", part->pdgId());
//                    //printf("pt = %7.3f, ", part->pt());
//                    //printf("eta = %7.3f, ", part->eta());
//                    //printf("\n");
//
//                    const reco::GenParticle *mom = getMother( *part );
//
//                    bool fid_cut = ( abs(part->eta())<jetEtaThreshold_ ) ? 1 : 0;
//                    bool status_cut = ( part->status()==22 || part->status()==23 ) ? 1 : 0; //22: intermediate state, 23: outgoing particles
//                    if (!fid_cut || !status_cut) continue;
//                    counter_num_gen_particles += 1;
//
//                    //printf("[check] ");
//                    //printf("(%2d) ", genLoop);
//                    //printf("status = %5d, ", part->status());
//                    //printf("pdgId = %9d, ", part->pdgId());
//                    //printf("pt = %7.3f, ", part->pt());
//                    //printf("eta = %7.3f, ", part->eta());
//                    //printf("mom status = %9d, ", mom->status());
//                    //printf("pdgId = %5d, ", mom->pdgId());
//                    //printf("pt = %7.3f, ", mom->pt());
//                    //printf("eta = %7.3f, ", mom->eta());
//                    //printf("\n");
//                    
//                    // store gen parton & mom kinematic info
//                    gens_pdgId.push_back(part->pdgId());
//                    gens_status.push_back(part->status());
//                    gens_pt.push_back(part->pt());
//                    gens_eta.push_back(part->eta());
//                    gens_phi.push_back(part->phi());
//                    gens_mass.push_back(part->mass());
//
//                    moms_pdgId.push_back(mom->pdgId());
//                    moms_status.push_back(mom->status());
//                    moms_pt.push_back(mom->pt());
//                    moms_eta.push_back(mom->eta());
//                    moms_phi.push_back(mom->phi());
//                    moms_mass.push_back(mom->mass());
//
//                    // gen matching
//                    std::vector<double> vec_deltaR;
//	                for(unsigned int index = 0 ; index < SelJetVect_PtSorted.size(); index++){
//                        edm::Ptr<reco::Jet> jet = SelJetVect_PtSorted[index];
//                        double deltaR = sqrt( (jet->eta()-part->eta())*(jet->eta()-part->eta()) + (jet->phi()-part->phi())*(jet->phi()-part->phi()) );
//                        vec_deltaR.push_back(deltaR);
//                    }
//                    std::vector<std::pair<int, double>> vec_deltaR_sorted = sortVector(vec_deltaR);
//                    if(vec_deltaR_sorted[0].second > 0.4) continue;
//                    vec_pdgId_register.push_back(part->pdgId());
//                    vec_index_register.push_back(vec_deltaR_sorted[0].first);
//                    vec_deltaR_register.push_back(vec_deltaR_sorted[0].second);
//                    //printf("part->pdgId() = %d\n", part->pdgId());
//                    //printf("vec_deltaR_sorted[0].first = %d\n", vec_deltaR_sorted[0].first);
//                    //printf("vec_deltaR_sorted[0].second = %.2f\n", vec_deltaR_sorted[0].second);
//                }
//
//                printf("[check] counter_num_gen_particles = %d\n", counter_num_gen_particles);
//
//                unsigned int register_size = vec_deltaR_register.size();
//                int    pdgId_jet0  = register_size > 0 ? vec_pdgId_register[0] : -1;
//                int    pdgId_jet1  = register_size > 1 ? vec_pdgId_register[1] : -1;
//                int    pdgId_jet2  = register_size > 2 ? vec_pdgId_register[2] : -1;
//                int    index_jet0  = register_size > 0 ? vec_index_register[0] : -1;
//                int    index_jet1  = register_size > 1 ? vec_index_register[1] : -1;
//                int    index_jet2  = register_size > 2 ? vec_index_register[2] : -1;
//                double deltaR_jet0 = register_size > 0 ? vec_deltaR_register[0] : -1.;
//                double deltaR_jet1 = register_size > 1 ? vec_deltaR_register[1] : -1.;
//                double deltaR_jet2 = register_size > 2 ? vec_deltaR_register[2] : -1.;
//
//                //printf("[check] vec_deltaR_register size = %lu\n", vec_deltaR_register.size());
//                //printf("(jet0, jet1, jet2) = (%d, %d, %d)\n", index_jet0, index_jet1, index_jet2);
//                //printf("(jet0, jet1, jet2) = (%d, %d, %d)\n", pdgId_jet0, pdgId_jet1, pdgId_jet2);
//                //printf("(jet0, jet1, jet2) = (%.2f, %.2f, %.2f)\n", deltaR_jet0, deltaR_jet1, deltaR_jet2);
//                //printf("\n");
//                std::vector<int > jets_matched_pdgId = { pdgId_jet0, pdgId_jet1, pdgId_jet2 };
//                std::vector<int > jets_matched_index = { index_jet0, index_jet1, index_jet2 };
//                std::vector<double > jets_matched_deltaR = { deltaR_jet0, deltaR_jet1, deltaR_jet2 };
//                thqhtags_obj.setMatchedPdgId  ( jets_matched_pdgId  );
//                thqhtags_obj.setMatchedIndex  ( jets_matched_index  );
//                thqhtags_obj.setMatchedDeltaR ( jets_matched_deltaR );
//                thqhtags_obj.setGenPdgId      ( gens_pdgId          );
//                thqhtags_obj.setGenStatus     ( gens_status         );
//                thqhtags_obj.setGenPt         ( gens_pt             );
//                thqhtags_obj.setGenEta        ( gens_eta            );
//                thqhtags_obj.setGenPhi        ( gens_phi            );
//                thqhtags_obj.setGenMass       ( gens_mass           );
//                thqhtags_obj.setMomPdgId      ( moms_pdgId          );
//                thqhtags_obj.setMomStatus     ( moms_status         );
//                thqhtags_obj.setMomPt         ( moms_pt             );
//                thqhtags_obj.setMomEta        ( moms_eta            );
//                thqhtags_obj.setMomPhi        ( moms_phi            );
//                thqhtags_obj.setMomMass       ( moms_mass           );
//            }
//--------------------------------------------------------------------------------------------//

            thqhtags->push_back( thqhtags_obj );

            if (evt.isRealData()) {
                thqhtags_obj.setStage1recoTag( DiPhotonTagBase::stage1recoTag::RECO_THQ_LEP );
                thqhtags->push_back( thqhtags_obj ); //FIXME at next iteration!!
            }

        } else {
            if(false) std::cout << " THQHadronicTagProducer NO TAG " << std::endl;
        }// end of photon if-else statement
    }//diPho loop end !
    evt.put( std::move( thqhtags ) );
}// end of THQHadronicTagProducer::produce
}
typedef flashgg::THQHadronicTagProducer FlashggTHQHadronicTagProducer;
DEFINE_FWK_MODULE( FlashggTHQHadronicTagProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
