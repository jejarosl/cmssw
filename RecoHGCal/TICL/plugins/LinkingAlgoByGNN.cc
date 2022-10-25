
#include <cmath>
#include <string>
#include "RecoHGCal/TICL/plugins/LinkingAlgoByGNN.h"

#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"

using namespace ticl;
using namespace cms::Ort;

LinkingAlgoByGNN::LinkingAlgoByGNN(const edm::ParameterSet &conf)
    : LinkingAlgoBase(conf),
      del_tk_ts_layer1_(conf.getParameter<double>("delta_tk_ts_layer1")),
      del_tk_ts_int_(conf.getParameter<double>("delta_tk_ts_interface")),
      del_ts_em_had_(conf.getParameter<double>("delta_ts_em_had")),
      del_ts_had_had_(conf.getParameter<double>("delta_ts_had_had")),
      timing_quality_threshold_(conf.getParameter<double>("track_time_quality_threshold")),
      cutTk_(conf.getParameter<std::string>("cutTk")) {}

LinkingAlgoByGNN::~LinkingAlgoByGNN() {}

void LinkingAlgoByGNN::initialize(const HGCalDDDConstants *hgcons,
                                                 const hgcal::RecHitTools rhtools,
                                                 const edm::ESHandle<MagneticField> bfieldH,
                                                 const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;

  bfield_ = bfieldH;
  propagator_ = propH;
}
void LinkingAlgoByGNN::linkTracksters(const edm::Handle<std::vector<reco::Track>> tkH,
                                                     const edm::ValueMap<float> &tkTime,
                                                     const edm::ValueMap<float> &tkTimeErr,
                                                     const edm::ValueMap<float> &tkTimeQual,
                                                     const std::vector<reco::Muon> &muons,
                                                     const edm::Handle<std::vector<Trackster>> tsH,
                                                     std::vector<TICLCandidate> &resultLinked,
                                                     std::vector<TICLCandidate> &chargedHadronsFromTk, const ONNXRuntime* cache) {
	std::cout << "Linking Algo by GNN " << std::endl;
  const auto &tracks = *tkH;
  const auto &tracksters = *tsH;

  auto bFieldProd = bfield_.product();
  const Propagator &prop = (*propagator_);
  const std::vector<std::string> input_names = {"features", "edge_index"};

  const auto N = 5;
  const auto shapeFeatures = 16;
  
	FloatArrays data;
	std::vector<std::vector<int64_t>> input_shapes;

	input_shapes.push_back({N, shapeFeatures});
	for(int i = 0; i != N; i++){
	data.emplace_back(shapeFeatures,0.1);
  //const  std::vector<float> vec1(N, 0.1);
	}
	
	data.emplace_back(3*N, 0);
	data.emplace_back(3*N, 1);
  
	input_shapes.push_back({2, 3*N});
	
  std::vector<float> outputs(3*N, 0.);  // init as all zeros

	outputs = cache->run(input_names, data, input_shapes)[0];
  std::cout << "OUPUTS SHAPE 	 " << outputs.size() << std::endl;

	for(auto & values : outputs){
		
		std::cout << "Val " << values << std::endl;
	
	}

}  // linkTracksters

void LinkingAlgoByGNN::fillPSetDescription(edm::ParameterSetDescription &desc) {
  desc.add<std::string>("cutTk",
                        "1.48 < abs(eta) < 3.0 && pt > 1. && quality(\"highPurity\") && "
                        "hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 5");
  desc.add<double>("delta_tk_ts_layer1", 0.02);
  desc.add<double>("delta_tk_ts_interface", 0.03);
  desc.add<double>("delta_ts_em_had", 0.03);
  desc.add<double>("delta_ts_had_had", 0.03);
  desc.add<double>("track_time_quality_threshold", 0.5);
  LinkingAlgoBase::fillPSetDescription(desc);
}
