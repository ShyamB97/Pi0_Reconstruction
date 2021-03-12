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
#include "larcore/Geometry/Geometry.h"

//protoDUNE analysis headers
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNECalibration.h"

//dunetpc headers
#include "dunetpc/dune/DuneObj/ProtoDUNEBeamEvent.h"

//ROOT includes
#include <TTree.h>
#include <TH1D.h>
#include "Math/Vector3D.h"

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
  std::vector<double> CNNScoreCalculator(anab::MVAReader<recob::Hit,4> &hitResults, const std::vector< art::Ptr< recob::Hit > > &hits, unsigned int &n);
  std::vector<double> StartHitQuantityCalculator(TVector3 &hitStart, TVector3 &hit, TVector3 &direction);
  void reset();

  private:

  // fcl parameters, order matters!
  protoana::ProtoDUNECalibration calibration_SCE;
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fHitTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  std::string fOpFlashTag;
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;

  // local variables
  bool fVerbose;

  int nomialMomentum = 1;
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils; // get BeamLineUtils class
  TTree *fOutTree = new TTree;
  std::vector<int> pdgs; // particle ids

  std::vector<int> pandoraTags; // track/shower like tag from pandora

  std::vector<double> emScore;
  std::vector<double> trackScore;
  std::vector<double> CNNScore; // CNN score per shower

  // shower start position
  std::vector<double> startPosX;
  std::vector<double> startPosY;
  std::vector<double> startPosZ;

  // shower direction
  std::vector<double> dirX;
  std::vector<double> dirY;
  std::vector<double> dirZ;

  std::vector<double> energy; // reco shower energy in ???
  std::vector<double> energyMC; // mc shower energy in GeV
  std::vector<int> nHits; // number of collection plane hits

  // quantity used to calculate the number of start hits
  std::vector<std::vector<double>> hitRadial;
  std::vector<std::vector<double>> hitLongitudinal;

  // beam start position
  double beamStartPosX;
  double beamStartPosY;
  double beamStartPosZ;

  // beam end position
  double beamEndPosX;
  double beamEndPosY;
  double beamEndPosZ;

  // true start position
  std::vector<double> trueStartPosX;
  std::vector<double> trueStartPosY;
  std::vector<double> trueStartPosZ;

  // true end position
  std::vector<double> trueEndPosX;
  std::vector<double> trueEndPosY;
  std::vector<double> trueEndPosZ;

  unsigned int eventID;
};


protoana::pi0TestSelection::pi0TestSelection(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  calibration_SCE(p.get<fhicl::ParameterSet>("CalibrationParsSCE")),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fHitTag(p.get<std::string>("HitTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fVerbose(p.get<bool>("Verbose")),
  fBeamlineUtils(p.get<fhicl::ParameterSet>("BeamlineUtils"))
{ }


// Calculates the CNN score of the event. Used to determine if a an event is track or shower like.
std::vector<double> protoana::pi0TestSelection::CNNScoreCalculator(anab::MVAReader<recob::Hit,4> &hitResults, const std::vector< art::Ptr< recob::Hit > > &hits, unsigned int &n)
{
  std::vector<double> output{};
  double score = 0;
  double mean_em = 0;
  double mean_track = 0;
  double cnn_track = 0;
  double cnn_em = 0;

  // Calculate the score per hit than take the average
  for(unsigned int h = 0; h < n; h++)
  {
    std::array<float,4> cnn_out = hitResults.getOutput( hits[h] );
    cnn_track = cnn_out[ hitResults.getIndex("track") ];
    cnn_em = cnn_out[ hitResults.getIndex("em") ];

    score += cnn_em / (cnn_em + cnn_track);
    mean_em += cnn_em;
    mean_track += cnn_track;
  }
  score = (n > 0) ? ( score / n ) : -999; // posts -999. if there were no hits
  mean_em = (n > 0) ? ( mean_em / n ) : -999;
  mean_track = (n > 0) ? ( mean_track / n ) : -999;

  output.push_back(score);
  output.push_back(mean_em);
  output.push_back(mean_track);

  return output;
}


// calculates the quantities for determining hits close to the shower start
std::vector<double> protoana::pi0TestSelection::StartHitQuantityCalculator(TVector3 &hitStart, TVector3 &hit, TVector3 &direction)
{
  std::vector<double> output{};  // make the return type a pair!
  TVector3 cross = (hit - hitStart).Cross(hitStart); // rsin(theta) compare to cylinder radius
  double dot = (hit - hitStart).Dot(hitStart); // rcos(theta) compare to cylinder length
  output.push_back(cross.Mag());
  output.push_back(dot);
  return output;
}

void protoana::pi0TestSelection::beginJob()
{
  // intiialize output root file
  art::ServiceHandle<art::TFileService> tfs;
  fOutTree = tfs->make<TTree>("tree", "");
  //Once we are done, stream results into the ROOT file
  fOutTree->Branch("pdgs", &pdgs);
  fOutTree->Branch("EventID", &eventID);
  fOutTree->Branch("pandoraTag", &pandoraTags);
  fOutTree->Branch("reco_daughter_PFP_emScore", &emScore);
  fOutTree->Branch("reco_daughter_PFP_trackScore", &trackScore);
  fOutTree->Branch("CNNScore", &CNNScore);
  fOutTree->Branch("reco_daughter_PFP_nHits_collection", &nHits);
  fOutTree->Branch("reco_daughter_allShower_startX", &startPosX);
  fOutTree->Branch("reco_daughter_allShower_startY", &startPosY);
  fOutTree->Branch("reco_daughter_allShower_startZ", &startPosZ);
  fOutTree->Branch("reco_daughter_allShower_dirX", &dirX);
  fOutTree->Branch("reco_daughter_allShower_dirY", &dirY);
  fOutTree->Branch("reco_daughter_allShower_dirZ", &dirZ);
  fOutTree->Branch("reco_daughter_allShower_energy", &energy);
  fOutTree->Branch("hitRadial", &hitRadial);
  fOutTree->Branch("hitLongitudinal", &hitLongitudinal);
}


void protoana::pi0TestSelection::reset()
{
  pdgs.clear();
  pandoraTags.clear();

  emScore.clear();
  trackScore.clear();
  CNNScore.clear();

  startPosX.clear();
  startPosY.clear();
  startPosZ.clear();

  dirX.clear();
  dirY.clear();
  dirZ.clear();

  energy.clear();
  //energyMC.clear();
  nHits.clear();
  hitRadial.clear();
  hitLongitudinal.clear();

  //trueStartPosX.clear();
  //trueStartPosY.clear();
  //trueStartPosZ.clear();
  
  //trueEndPosX.clear();
  //trueEndPosY.clear();
  //trueEndPosZ.clear();

}


void protoana::pi0TestSelection::analyze(art::Event const & evt)
{
  std::cout << "module running..." << std::endl;
  // clear any outputs that are lists
  reset();

  //Initialise util classes
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNETrackUtils trkUtil;
  protoana::ProtoDUNETruthUtils truthUtil;
  protoana::ProtoDUNEShowerUtils showerUtil;
  protoana::ProtoDUNETrackUtils trackUtil;

  std::cout << "getting handle" << std::endl;
  auto pfpVec = evt.getValidHandle< std::vector< recob::PFParticle > >( fPFParticleTag );
  std::cout << "got handle" << std::endl;
  eventID = evt.id().event();

  std::cout << "getting clockData" << std::endl;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  std::cout << "got clockData" << std::endl;
  auto const detProp =  art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);

  // Get only Beam particle by checking the Beam slices
  std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt, fPFParticleTag);
  if(beamParticles.size() != 1)
  {
    std::cout << "there shouldn't be more than one beam particle" << std::endl;
    return;
  }
  auto beamParticle = beamParticles[0];
  pdgs.push_back(beamParticle->PdgCode()); // get ID of particles

  const recob::Track* beamTrack = 0x0;
  beamTrack = pfpUtil.GetPFParticleTrack(*beamParticle, evt, fPFParticleTag, fTrackerTag);

  if(!beamTrack)
  {
    std::cout<< "no beam track found, moving on" << std::endl;
    beamStartPosX = -999;
    beamStartPosY = -999;
    beamStartPosZ = -999;
    beamEndPosX = -999;
    beamEndPosY = -999;
    beamEndPosZ = -999;
  }
  else
  {
    beamStartPosX = beamTrack->Trajectory().Start().X();
    beamStartPosY = beamTrack->Trajectory().Start().Y();
    beamStartPosZ = beamTrack->Trajectory().Start().Z();
    beamEndPosX = beamTrack->Trajectory().End().X();
    beamEndPosY = beamTrack->Trajectory().End().Y();
    beamEndPosZ = beamTrack->Trajectory().End().Z();

    // if the beam enters from the opposite direction
    if(beamStartPosZ > beamEndPosZ)
    {
      beamStartPosX = beamTrack->Trajectory().End().X();
      beamStartPosY = beamTrack->Trajectory().End().Y();
      beamStartPosZ = beamTrack->Trajectory().End().Z();
      beamEndPosX = beamTrack->Trajectory().Start().X();
      beamEndPosY = beamTrack->Trajectory().Start().Y();
      beamEndPosZ = beamTrack->Trajectory().Start().Z();
    }
  }

  anab::MVAReader<recob::Hit,4> hitResults(evt, "emtrkmichelid:emtrkmichel" );

  // get daughterPFP objects using for loops
  for( size_t daughterID : beamParticle->Daughters() )
  {
    const recob::PFParticle * daughterPFP = &(pfpVec->at( daughterID ));
    // determine if they are track like or shower like using pandora
    // then fill a vector containing this data: 11 = shower 13 = track
    if(pfpUtil.IsPFParticleShowerlike(*daughterPFP, evt, fPFParticleTag, fShowerTag))
    {
      pandoraTags.push_back(11);
    }
    if(pfpUtil.IsPFParticleTracklike(*daughterPFP, evt, fPFParticleTag, fTrackerTag))
    {
      pandoraTags.push_back(13);
    }

    // number of collection plane hits
    auto hits = pfpUtil.GetPFParticleHitsFromPlane_Ptrs( *beamParticle, evt, fPFParticleTag, 2 );
    unsigned int num = hits.size();
    nHits.push_back(num);
    std::cout << "PFP hits: " << num << std::endl;
    
    // calculate cnn score
    std::vector<double> cnnOutput = CNNScoreCalculator(hitResults, hits, num);
    CNNScore.push_back(cnnOutput[0]);
    emScore.push_back(cnnOutput[1]);
    trackScore.push_back(cnnOutput[2]);

    // force daughters to be shower
    const recob::Shower* shower = 0x0;
    std::cout << "Getting shower" << std::endl;
    try
    {
      shower =	pfpUtil.GetPFParticleShower(*daughterPFP, evt, fPFParticleTag, "pandora2Shower");
      art::FindManyP<recob::SpacePoint> spFromHits(hits, evt, fHitTag);

      if(shower)
      {
        std::cout << "got shower" << std::endl;

        std::cout << "getting start and direction" << std::endl;
        TVector3 showerStart = shower->ShowerStart();
        TVector3 showerDir = shower->Direction();
        std::cout << "got start and direction" << std::endl;

        startPosX.push_back(showerStart.X());
        startPosY.push_back(showerStart.Y());
        startPosZ.push_back(showerStart.Z());

        dirX.push_back(showerDir.X());
        dirY.push_back(showerDir.Y());
        dirZ.push_back(showerDir.Z());

        std::vector<double> hitRad;
        std::vector<double> hitLong;
        for(unsigned int n = 0; n < hits.size(); n++)
        {
          std::vector<art::Ptr<recob::SpacePoint>> sps = spFromHits.at(n);

          if(!sps.empty())
          {
            TVector3 hitPoint(sps[0]->XYZ()[0], sps[0]->XYZ()[1], sps[0]->XYZ()[2]);

            std::vector<double> startHitQuantities = StartHitQuantityCalculator(showerStart, hitPoint, showerDir);

            hitRad.push_back(startHitQuantities[0]);
            hitLong.push_back(startHitQuantities[1]);
          }
          else
          {
            hitRad.push_back(-999);
            hitLong.push_back(-999);
          }
        }

        hitRadial.push_back(hitRad);
        hitLongitudinal.push_back(hitLong);

        std::vector<double> x_vec, y_vec, z_vec;
        double total_y = 0;
        int n_good_y = 0;
        std::vector<art::Ptr<recob::Hit>> good_hits;
        for(unsigned int i = 0; i < hits.size(); i++)
        {
          auto hit = hits[i];
          // skip any hits not on the ocllection plane (shouldn't be anyways)
          if(hit->View() != 2)
          {
            continue;
          }

          good_hits.push_back(hit);

          double shower_hit_x = detProp.ConvertTicksToX(hit->PeakTime(), hit->WireID().Plane, hit->WireID().TPC, 0);
          double shower_hit_z = geom->Wire(hit->WireID()).GetCenter().Z();

          x_vec.push_back(shower_hit_x);
          z_vec.push_back(shower_hit_z);

          std::vector<art::Ptr<recob::SpacePoint>> sps = spFromHits.at(i);
          if (!sps.empty())
          {
            y_vec.push_back(sps[0]->XYZ()[1]);
            total_y += y_vec.back();
            n_good_y++;
          }
          else
          {
            y_vec.push_back(-999.);
          }
        }

        if(n_good_y < 1)
        {
          std::cout << "could not reconstruct energy" << std::endl;
          energy.push_back(-999);
        }
        else
        {
          double total_energy = 0;
          for(unsigned int j = 0; j < good_hits.size(); j++)
          {
            auto good_hit = good_hits[j];
            
            if(good_hit->View() != 2)
            {
              continue;
            }

            if (y_vec[j] < -100.)
            {
              y_vec[j] = total_y / n_good_y;
            }
            total_energy += calibration_SCE.HitToEnergy(good_hit, x_vec[j], y_vec[j], z_vec[j]);
          }
          energy.push_back(total_energy);
        }
      }
      else
      {
        startPosX.push_back(-999);
        startPosY.push_back(-999);
        startPosZ.push_back(-999);
        dirX.push_back(-999);
        dirY.push_back(-999);
        dirZ.push_back(-999);
        energy.push_back(-999);
      }
    }
    catch( const cet::exception &e )
    {
      std::cout << "couldn't get shower object! Moving on" << std::endl;
      continue;
    }
  }

  // store any MC reated information here i.e. MC truth info
  if(!evt.isRealData())
  {
    for( size_t daughterID : beamParticle->Daughters() )
    {
      const recob::PFParticle * daughterPFP = &(pfpVec->at( daughterID ));
      //const sim::ParticleList & plist = pi_serv->ParticleList();
      std::cout << "getting shared hits" << std::endl;
      protoana::MCParticleSharedHits match = truthUtil.GetMCParticleByHits( clockData, *daughterPFP, evt, fPFParticleTag, fHitTag );
      std::cout << "got shared hits" << std::endl;
      if(match.particle)
      {
        std::cout << "we have matched the MC particles!" << std::endl;
      }
      else
      {
        std::cout << "MC particle not matched" << std::endl;
      }

    }
  }

  fOutTree->Fill();
}

void protoana::pi0TestSelection::endJob()
{

}

DEFINE_ART_MODULE(protoana::pi0TestSelection)
