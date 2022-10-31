
#include <cmath>
#include <string>
#include "RecoHGCal/TICL/plugins/LinkingAlgoByGNN.h"

#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"

#include <bits/stdc++.h>

using namespace std;
using namespace ticl;
using namespace cms::Ort;

class Graph {
  // A function used by DFS
  void DFSUtil(int v);

public:
  map<int, bool> visited;
  map<int, list<int>> adj;
  std::vector<std::vector<int>> connected_components;
  // function to add an edge to graph
  void addEdge(int v, int w);

  // prints DFS traversal of the complete graph
  void DFS();
};

void Graph::addEdge(int v, int w) {
  adj[v].push_back(w);  // Add w to v�s list.
}

void Graph::DFSUtil(int v) {
  // Mark the current node as visited and print it
  visited[v] = true;
  std::cout << v << std::endl;

  // Recur for all the vertices adjacent to this vertex
  list<int>::iterator i;
  for (auto i = adj[v].begin(); i != adj[v].end(); ++i)
    if (!visited[*i]) {
      connected_components.back().push_back(*i);
      std::cout << "Pushed back " << *i << std::endl;
      DFSUtil(*i);
    }
}

// The function to do DFS traversal. It uses recursive
// DFSUtil()
void Graph::DFS() {
  // Call the recursive helper function to print DFS
  // traversal starting from all vertices one by one
  for (auto i : adj)
    if (visited[i.first] == false) {
      std::cout << "Emplaced back: " << i.first << std::endl;
      connected_components.emplace_back(1, i.first);  // {i.first}
      std::cout << "Starting DFS from node: " << i.first << std::endl;
      DFSUtil(i.first);
    }
}

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
                                      //const edm::Handle<std::vector<TICLGraph>> graph,
                                      std::vector<TICLCandidate> &resultLinked,
                                      std::vector<TICLCandidate> &chargedHadronsFromTk,
                                      const ONNXRuntime *cache) {
  std::cout << "Linking Algo by GNN " << std::endl;
  const auto &tracks = *tkH;
  const auto &tracksters = *tsH;

  auto bFieldProd = bfield_.product();
  const Propagator &prop = (*propagator_);
  std::cout << "LINKING ALGO GNN" << std::endl;

  const std::vector<std::string> input_names = {"features", "edge_index"};

  long int N = tracksters.size();
  const auto shapeFeatures = 16;

  FloatArrays data;
  std::vector<std::vector<int64_t>> input_shapes;

  std::vector<float> features;

  // Print out info about tracksters
  std::cout << "Number of tracksters in event: " << N << std::endl;
  for (unsigned i = 0; i < tracksters.size(); ++i) {
    const auto &ts = tracksters[i];
    std::cout << "Trackster " << i << "--------------------" << std::endl;
    std::cout << "Barycenter X: " << ts.barycenter().x() << std::endl;
    std::cout << "Barycenter Y: " << ts.barycenter().y() << std::endl;
    std::cout << "Barycenter Z: " << ts.barycenter().z() << std::endl;

    Vector eigenvector0 = ts.eigenvectors(0);
    std::cout << "eVector0 X: " << eigenvector0.x() << std::endl;
    std::cout << "eVector0 Y: " << eigenvector0.y() << std::endl;
    std::cout << "eVector0 Z: " << eigenvector0.z() << std::endl;

    std::array<float, 3> eigenvalues = ts.eigenvalues();
    std::cout << "EV1: " << eigenvalues[0] << std::endl;
    std::cout << "EV2: " << eigenvalues[1] << std::endl;
    std::cout << "EV3: " << eigenvalues[2] << std::endl;

    std::array<float, 3> sigmasPCA = ts.sigmasPCA();
    std::cout << "sigmaPCA1: " << sigmasPCA[0] << std::endl;
    std::cout << "sigmaPCA2: " << sigmasPCA[1] << std::endl;
    std::cout << "sigmaPCA3: " << sigmasPCA[2] << std::endl;

    // size
    std::cout << "Size: " << ts.vertices().size() << std::endl;

    std::cout << "Raw Energy: " << ts.raw_energy() << std::endl;
    std::cout << "Raw EM Energy: " << ts.raw_em_energy() << std::endl;

    // candidate labels
    /* Python:
    
    in_candidate = [-1 for i in range(trk_data[ev].NTracksters)]
            for indx, cand in enumerate(cand_data[ev].tracksters_in_candidate):
                for ts in cand:
                    in_candidate[ts] = indx
    
    */
    features.push_back(ts.barycenter().x());
    features.push_back(ts.barycenter().y());
    features.push_back(ts.barycenter().z());
    features.push_back(eigenvector0.x());
    features.push_back(eigenvector0.y());
    features.push_back(eigenvector0.z());
    features.push_back(eigenvalues[0]);
    features.push_back(eigenvalues[1]);
    features.push_back(eigenvalues[2]);
    features.push_back(sigmasPCA[0]);
    features.push_back(sigmasPCA[1]);
    features.push_back(sigmasPCA[2]);
    features.push_back(ts.vertices().size());
    features.push_back(ts.raw_energy());
    features.push_back(ts.raw_em_energy());
    features.push_back(-1.0);  // set candidates index to -1 as of now; remove later (after changing the model file)

    std::cout << "--------------------" << std::endl;
  }

  input_shapes.push_back({1, N, shapeFeatures});
  data.emplace_back(features);

  // Creating Edges: uncomment when have a Graph as an input

  //std::vector<float_t> edges_src;
  //std::vector<float_t> edges_dst;
  //for (int i = 0; i < N; i++){
  //  for (auto & i_neighbour : graph.node_linked_inners[i]){
  //    // Create an edge between the tracksters
  //    edges_src.push_back(i_neighbour);
  //    edges_dst.push_back(i);
  //  }
  //}

  // Create fully connected graph for testing
  std::vector<float> edges_src;
  std::vector<float> edges_dst;

  for (int i = 0; i < N; i++) {
    std::cout << "i: " << i << std::endl;
    for (int j = i; j < N; j++) {
      std::cout << "j: " << j << std::endl;
      edges_src.push_back(static_cast<float>(i));
      edges_dst.push_back(static_cast<float>(j));
    }
  }

  long unsigned int numEdges = edges_src.size();
  input_shapes.push_back({1, 2, static_cast<int>(numEdges)});
  std::cout << "Num edges: " << numEdges << std::endl;

  data.emplace_back(edges_src);
  for (auto &dst : edges_dst) {
    data.back().push_back(dst);
  }

  std::vector<float> edge_predictions = cache->run(input_names, data, input_shapes)[0];

  std::cout << "Network output shape is " << edge_predictions.size() << std::endl;

  for (long unsigned int i = 0; i < edge_predictions.size(); i++) {
    std::cout << "Network output for edge " << data[1][i] << "-" << data[1][numEdges + i]
              << " is: " << edge_predictions[i] << std::endl;
  }

  // Create a graph
  Graph g;
  const auto classification_threshold = 0.7;

  // Building a predicted graph
  for (long unsigned int i = 0; i < numEdges; i++) {
    if (edge_predictions[i] >= classification_threshold) {
      auto src = data[1][i];
      auto dst = data[1][numEdges + i];
      // Make undirectional
      g.addEdge(src, dst);
      g.addEdge(dst, src);
    }
  }

  std::cout << "HERE 8" << std::endl;

  std::cout << "Following is Depth First Traversal" << std::endl;
  std::cout << "Connected components are: " << std::endl;
  g.DFS();

  int i = 0;
  std::vector<TICLCandidate> connectedCandidates;

  for (auto &component : g.connected_components) {
    TICLCandidate tracksterCandidate;
    for (auto &trackster_id : component) {
      std::cout << "Component " << i << ": trackster id " << trackster_id << std::endl;
      tracksterCandidate.addTrackster(edm::Ptr<Trackster>(tsH, trackster_id));
    }
    i++;
    connectedCandidates.push_back(tracksterCandidate);
  }

  // The final candidates are passed to `resultLinked`
  resultLinked.insert(std::end(resultLinked), std::begin(connectedCandidates), std::end(connectedCandidates));

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
