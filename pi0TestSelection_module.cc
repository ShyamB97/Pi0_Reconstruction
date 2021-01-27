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
{

}

void protoana::pi0TestSelection::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fOutTree = tfs->make<TTree>("tree", "");
  fOutTree->Branch("pdgs", &pdgs);
  fOutTree->Branch("EventID", &eventID);
}

void protoana::pi0TestSelection::analyze(art::Event const & evt)
{
  std::cout << "module running..." << std::endl;

  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNETrackUtils trkUtil;

  eventID = evt.id().event();

  if(!evt.isRealData())
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
		std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt, fPFParticleTag);

		if(beamParticles.size() > 1)
		{
			std::cout << "We have more tha one Beam particle!" << std::endl;
			return;
		}
		auto beamParticle = beamParticles[0];
		//Verifying pdg code from beamEvent to PFP, should be equal for MC data
		int pdgFromPFP = beamParticle->PdgCode();
		if(pdgFromPFP == pdgs[0])
		{
			std::cout << "PFP code matches beamline code" << std::endl;
		}
  }

  fOutTree->Fill();
}

void protoana::pi0TestSelection::endJob()
{

}

DEFINE_ART_MODULE(protoana::pi0TestSelection)
