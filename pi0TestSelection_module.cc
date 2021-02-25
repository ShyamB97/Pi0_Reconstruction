////////////////////////////////////////////////////////////////////////
// Class:       pi0TestSelection
// Plugin Type: analyzer (art v2_07_03)
// File:        pi0TestSelection_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/OpFlash.h"

//protoDUNE analysis headers
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"

//dunetpc headers
#include "dunetpc/dune/DuneObj/ProtoDUNEBeamEvent.h"

//ROOT includes
#include <TTree.h>
#include <TH1D.h>

//Use this to get CNN number
#include "lardata/ArtDataHelper/MVAReader.h"

namespace protoana {
  class pi0TestSelection;
}


class protoana::pi0TestSelection : public art::EDAnalyzer {
public:

  explicit pi0TestSelection(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  pi0TestSelection(pi0TestSelection const &) = delete;
  pi0TestSelection(pi0TestSelection &&) = delete;
  pi0TestSelection & operator = (pi0TestSelection const &) = delete;
  pi0TestSelection & operator = (pi0TestSelection &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Custom functions.
  double CNNScoreCalculator(anab::MVAReader<recob::Hit,4> &hitResults, const std::vector< art::Ptr< recob::Hit > > &hits, unsigned int &n);

private:

  // fcl parameters
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  std::string fOpFlashTag;


  // local variables
  bool fVerbose;

  int nomialMomentum = 1;
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils; // get BeamLineUtils class
  TTree *fOutTree = new TTree;
  std::vector<int> pdgs; // particle ids
  double CNNScore; // global CNN track/em like score.
  std::vector<double> showerCNNScores; // CNN score per shower
  unsigned int nHits; // number of collection plane hits
  double E_MC; // Energy of MC particle (beam) in GeV
  std::vector<double> E_Showers; // Total energy deposited per shower

  unsigned int eventID;
};


protoana::pi0TestSelection::pi0TestSelection(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fVerbose(p.get<bool>("Verbose")),
  fBeamlineUtils(p.get<fhicl::ParameterSet>("BeamlineUtils"))
{ }

// Calculates the CNN score of the event. Used to determine if a an event is track or shower like.
double protoana::pi0TestSelection::CNNScoreCalculator(anab::MVAReader<recob::Hit,4> &hitResults, const std::vector< art::Ptr< recob::Hit > > &hits, unsigned int &n)
{
  double score = 0;
  double cnn_track = 0;
  double cnn_em = 0;

  // Calculate the score per hit than take the average
  for(unsigned int h = 0; h < n; h++)
  {
    std::array<float,4> cnn_out = hitResults.getOutput( hits[h] );
    cnn_track = cnn_out[ hitResults.getIndex("track") ];
    cnn_em = cnn_out[ hitResults.getIndex("em") ];

    score += cnn_em / (cnn_em + cnn_track);
  }
  score = (nHits > 0) ? ( score / n ) : -999.; // posts -999. if there were no hits
  return score;
}

/*
void Prod4BeamParticle()
{
  std::cout << "we have MC data" << std::endl;
  auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("generator");
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent> > beamVec;
  art::fill_ptr_vector(beamVec, beamHandle);
  const beam::ProtoDUNEBeamEvent& beamEvent = *(beamVec.at(0));

  std::cout << "Timing Trigger: " << beamEvent.GetTimingTrigger() << std::endl;  // print out timing tigger
  std::cout << "matched with database: " << beamEvent.CheckIsMatched() << std::endl; // check if beam event was matched with existing database
  std::cout << "got beamEvent!" << std::endl;
  
  pdgs = fBeamlineUtils.GetPID(beamEvent, nomialMomentum); // get ID of particles
  std::cout << "got pdg from beam event" << std::endl;

  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt, fPFParticleTag);
  std::cout << "got beam PFP" << std::endl;

  std::cout << "checking PFP vector size" << std::endl;
  unsigned int numOfBeamParticles = beamParticles.size();
  if(numOfBeamParticles != 1)
  {
    std::cout << "number of beam particles: " << numOfBeamParticles << std::endl;
    if(numOfBeamParticles == 0)
    {
      std::cout << "no beam particles, moving on..." << std::endl;
      return;
    }
  }

  const recob::PFParticle* beamParticle = beamParticles[0];
  std::cout << "assigned beam particle" << std::endl;

  //Verifying pdg code from beamEvent to PFP, should be equal for MC data
  int pdgFromPFP = beamParticle->PdgCode();
  std::cout << "pdg code from PFP: " << pdgFromPFP << std::endl;
}
*/

void protoana::pi0TestSelection::beginJob()
{
  // intiialize output root file
  art::ServiceHandle<art::TFileService> tfs;
  fOutTree = tfs->make<TTree>("tree", "");
  //Once we are done, stream results into the ROOT file
  fOutTree->Branch("pdgs", &pdgs);
  fOutTree->Branch("EventID", &eventID);
  fOutTree->Branch("CNNScore", &CNNScore);
  fOutTree->Branch("showerCNNScore", &showerCNNScores);
  fOutTree->Branch("CollectionPlaneHits", &nHits);
}

void protoana::pi0TestSelection::analyze(art::Event const & evt)
{
  std::cout << "module running..." << std::endl;
  // clear any outputs that are lists
  pdgs.clear();
  showerCNNScores.clear();

  //Initialise util classes
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNETrackUtils trkUtil;
  protoana::ProtoDUNETruthUtils truthUtil;
  protoana::ProtoDUNEShowerUtils showerUtil;
  protoana::ProtoDUNETrackUtils trackUtil;

  eventID = evt.id().event();

  // Get only Beam particle by checking the Beam slices
  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt, fPFParticleTag);
  if(beamParticles.size() != 1)
  {
    std::cout << "there shouldn't be more than one beam particle" << std::endl;
    return;
  }
  auto beamParticle = beamParticles[0];
  pdgs.push_back(beamParticle->PdgCode()); // get ID of particles
  
  //=============================================================================================//
  // Get the beam Daughter tracks, showers.
  std::cout << "Beam daughters: " << beamParticle->NumDaughters() << std::endl;
  const std::vector<const recob::Shower*> daughterShowers = pfpUtil.GetPFParticleDaughterShowers(*beamParticle, evt, fPFParticleTag, fShowerTag);

  std::cout << "Number of daughter showers: " << daughterShowers.size() << std::endl;

  const std::vector<const recob::Track*> daughterTracks = pfpUtil.GetPFParticleDaughterTracks(*beamParticle, evt, fPFParticleTag, fTrackerTag);

  std::cout << "Number of daughter tracks: " << daughterTracks.size() << std::endl;

  //=============================================================================================//
  // Compute CNN scores

  // Helper to get hits and the 4 associated CNN outputs
  // CNN Outputs: EM, Track, Michel, Empty
  // outputNames: track, em, none, michel
  anab::MVAReader<recob::Hit,4> hitResults(evt, "emtrkmichelid:emtrkmichel" );

  // First get the CNN score per event (shower + tracks)
  // we want the hits from the collection plane only. i.e. plane 2 (I think)
  const std::vector< art::Ptr< recob::Hit > > hits = pfpUtil.GetPFParticleHitsFromPlane_Ptrs( *beamParticle, evt, fPFParticleTag, 2 );
  nHits = hits.size();

  CNNScore = CNNScoreCalculator(hitResults, hits, nHits); // compute the CNN_score
  std::cout << "global CNN track-em like score: " << CNNScore << std::endl;
  std::cout << "number of collection plane hits: " << nHits << std::endl;

  //Now get the CNN score per shower
  for(auto shower : daughterShowers)
  {
    const std::vector< art::Ptr<recob::Hit> > showerHits = showerUtil.GetRecoShowerArtHits(*shower, evt, fShowerTag);
    unsigned int num = showerHits.size();
    double score = CNNScoreCalculator(hitResults, showerHits, num);
    showerCNNScores.push_back(score);

  }

  if(!evt.isRealData())
  {
    //=============================================================================================//
    // Get MC particle and the energy
    const simb::MCParticle* trueParticle = 0x0;
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    trueParticle = truthUtil.GetMCParticleFromPFParticle(clockData, *beamParticle, evt, fPFParticleTag);
    if(trueParticle != 0x0)
    {
      std::cout << "MC particle matched by energy" << std::endl;
    }
    std::cout << "MC particle id: " << trueParticle->PdgCode() << std::endl;

    E_MC = trueParticle->E();
    std::cout << "MC particle energy (GeV): " << E_MC << std::endl;

    //=============================================================================================//
    // Get the energy per shower

    auto const detProp =  art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);

    for(auto shower : daughterShowers)
    {
      //const std::vector< const recob::Hit * > showerHits = showerUtil.GetRecoShowerHits(*shower, evt, fShowerTag);
      //std::vector<double> energyPerPlane = showerUtil.EstimateEnergyFromHitCharge(clockData, detProp, showerHits, calo::CalorimetryAlg caloAlg);
      //double totalEnergy = energyPerPlane[0] + energyPerPlane[1] + energyPerPlane[2];
      //E_showers.push_back(totalEnergy);
      const std::vector<double> E = shower->Energy();
      std::cout << "energy size?: " << E.size() << std::endl;
      //E_shower.push_back(shower->Energy()); // note sure how accurate this is...
    }

  }
  else
  {
    std::cout << "we have real data?" << std::endl;
  }

  fOutTree->Fill();
}

void protoana::pi0TestSelection::endJob()
{

}

DEFINE_ART_MODULE(protoana::pi0TestSelection)
