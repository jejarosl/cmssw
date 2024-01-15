// Author : Jekaterina Jaroslavceva - jejarosl@cern.ch
// Date : 01/2024
/*
TICL plugin for graph neural network-based superclustering in HGCAL, primarily focusing on hadronic interactions. 
DNN designed and trained by Jekaterina Jaroslavceva.

Inputs are CLUE3D tracksters. Outputs are superclusters (as vectors of IDs of trackster).

Algorithm description :
#TODO
*/


#include <map>
#include <list>

#include "RecoHGCal/TICL/interface/TracksterLinkingAlgoBase.h"
#include "RecoHGCal/TICL/interface/TICLGraph.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoHGCal/TICL/plugins/TracksterLinkingbyGNN.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"

using namespace ticl;


std::unique_ptr<TICLGraph> produce_ticl_graph(const auto tracksters) {

  TICLLayerTile tracksterTilePos;
  TICLLayerTile tracksterTileNeg;

  for (size_t id_t = 0; id_t < tracksters.size(); ++id_t) {
    auto t = tracksters[id_t];
    if (t.barycenter().eta() > 0.) {
      tracksterTilePos.fill(t.barycenter().eta(), t.barycenter().phi(), id_t);
    } else if (t.barycenter().eta() < 0.) {
      tracksterTileNeg.fill(t.barycenter().eta(), t.barycenter().phi(), id_t);
    }
  }

  std::vector<Node> allNodes;

  for (size_t id_t = 0; id_t < tracksters.size(); ++id_t) {
    auto t = tracksters[id_t];

    Node tNode(id_t);

    auto bary = t.barycenter();
    double del = 0.1;

    double eta_min = std::max(abs(bary.eta()) - del, (double)TileConstants::minEta);
    double eta_max = std::min(abs(bary.eta()) + del, (double)TileConstants::maxEta);

    if (bary.eta() > 0.) {
      std::array<int, 4> search_box =
          tracksterTilePos.searchBoxEtaPhi(eta_min, eta_max, bary.phi() - del, bary.phi() + del);
      if (search_box[2] > search_box[3]) {
        search_box[3] += TileConstants::nPhiBins;
      }

      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          auto &neighbours = tracksterTilePos[tracksterTilePos.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
          for (auto n : neighbours) {
            if (tracksters[n].barycenter().z() < bary.z()) {
              tNode.addInner(n);
            } else if (tracksters[n].barycenter().z() > bary.z()) {
              tNode.addOuter(n);
            }
          }
        }
      }
    }

    else if (bary.eta() < 0.) {
      std::array<int, 4> search_box =
          tracksterTileNeg.searchBoxEtaPhi(eta_min, eta_max, bary.phi() - del, bary.phi() + del);
      if (search_box[2] > search_box[3]) {
        search_box[3] += TileConstants::nPhiBins;
      }

      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          auto &neighbours = tracksterTileNeg[tracksterTileNeg.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
          for (auto n : neighbours) {
            if (abs(tracksters[n].barycenter().z()) < abs(bary.z())) {
              tNode.addInner(n);
            } else if (abs(tracksters[n].barycenter().z()) > abs(bary.z())) {
              tNode.addOuter(n);
            }
          }
        }
      }
    }
    allNodes.push_back(tNode);
  }
  return std::make_unique<TICLGraph>(allNodes);
}


// DFS Graph
class Graph {
  // A function used by DFS
  void DFSUtil(int v);

public:
  std::map<int, bool> visited;
  std::map<int, std::list<int>> adj;
  std::vector<std::vector<int>> connected_components;
  // function to add an edge to graph
  void addEdge(int v, int w);

  // prints DFS traversal of the complete graph
  void DFS();
};

void Graph::addEdge(int v, int w) {
  adj[v].push_back(w);  // Add w to v’s list.
}

void Graph::DFSUtil(int v) {
  // Mark the current node as visited and print it
  visited[v] = true;
  //std::cout << v << std::endl;

  // Recur for all the vertices adjacent to this vertex
  std::list<int>::iterator i;
  for (auto i = adj[v].begin(); i != adj[v].end(); ++i)
    if (!visited[*i]) {
      connected_components.back().push_back(*i);
      std::cout << "Pushed back " << *i << std::endl;
      DFSUtil(*i);
    }
}

// The function to do DFS traversal. It uses recursive DFSUtil()
void Graph::DFS() {
  // Call the recursive helper function to print DFS
  // traversal starting from all vertices one by one
  for (auto i : adj)
    if (visited[i.first] == false) {
      std::cout << "Emplaced back: " << i.first << std::endl;
      connected_components.emplace_back(1, i.first);
      std::cout << "Starting DFS from node: " << i.first << std::endl;
      DFSUtil(i.first);
    }
}

TracksterLinkingbyGNN::TracksterLinkingbyGNN(const edm::ParameterSet& ps, edm::ConsumesCollector iC, cms::Ort::ONNXRuntime const* onnxRuntime)
      : TracksterLinkingAlgoBase(ps, iC, onnxRuntime),
      nnVersion_(ps.getParameter<std::string>("nnVersion")),
      nnWorkingPoint_(ps.getParameter<double>("nnWorkingPoint")),
      deltaEtaWindow_(ps.getParameter<double>("deltaEtaWindow")),
      deltaPhiWindow_(ps.getParameter<double>("deltaPhiWindow"))
{
  assert(onnxRuntime_ && "TracksterLinkingbySuperClustering : ONNXRuntime was not provided, the model should have been set in onnxModelPath in the plugin config");
} 

/**
 * resultTracksters : should be all tracksters (including those not in Superclusters (SC)).
 * outputSuperclusters : indices into resultsTracksters. Should include tracksters not in SC as one-element vectors.
 * linkedTracksterIdToInputTracksterId : map indices from output to input.
*/
void TracksterLinkingbyGNN::linkTracksters(const Inputs& input, std::vector<Trackster>& resultTracksters,
                    std::vector<std::vector<unsigned int>>& outputSuperclusters,
                    std::vector<std::vector<unsigned int>>& linkedTracksterIdToInputTracksterId) {  
    
    std::cout << "Linking Algo by GNN" << std::endl;
    
    // Create trackster graph
    auto TICL_graph_result = produce_ticl_graph(input.tracksters);

    // Access the TICLGraph object through the unique pointer
    TICLGraph* ticlGraph = TICL_graph_result.get();
  
    const std::vector<std::string> input_names = {"features", "edge_index", "adj", "trackster_index"};
  
    long int N = input.tracksters.size();
    const auto shapeFeatures = 15;
  
  cms::Ort::FloatArrays data;
    std::vector<std::vector<int64_t>> input_shapes;
  
    std::vector<float> features;
    // For nearest neighbour finding
    std::vector<float> barycenters_x;
    std::vector<float> barycenters_y;
    std::vector<float> barycenters_z;
    
  
    // Print out info about tracksters
    std::cout << "Number of tracksters in event: " << N << std::endl;
    for (unsigned i = 0; i < N; ++i) {
      const auto ts = input.tracksters[i];
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
      
      barycenters_x.push_back(ts.barycenter().x());
      barycenters_y.push_back(ts.barycenter().y());
      barycenters_z.push_back(ts.barycenter().z());
  
      std::cout << "--------------------" << std::endl;
    }
  
    input_shapes.push_back({1, N, shapeFeatures});
    data.emplace_back(features);
  
    // Creating Edges of the input graphs
    std::vector<float> edges_src;
    std::vector<float> edges_dst;
    for (int i = 0; i < N; i++){
      for (auto & i_neighbour : ticlGraph->getNode(i).getInner()){
        // Create an edge between the tracksters
        edges_src.push_back(static_cast<float>(i_neighbour));
        edges_dst.push_back(static_cast<float>(i));
      }
    }
  
    auto numEdges = static_cast<int>(edges_src.size());
    std::cout << "Num edges: " << numEdges << std::endl;
      
    // Do not run the network if no edges in the event or less than two tracksters present
    if (numEdges == 0 || N < 2){
      linkedTracksterIdToInputTracksterId.resize(N);    
      std::vector<unsigned int> linkedTracksters;
      
      for (int trackster_id = 0; trackster_id < N; trackster_id++) {
        linkedTracksterIdToInputTracksterId[trackster_id].push_back(trackster_id);
        linkedTracksters.push_back(resultTracksters.size());
        resultTracksters.push_back(input.tracksters[trackster_id]);
        outputSuperclusters.push_back(linkedTracksters);
      }
      return;
    }
    
    input_shapes.push_back({1, 2, numEdges});
  
    data.emplace_back(edges_src);
    for (auto &dst : edges_dst) {
      data.back().push_back(dst);
    }
  
    
    // Creating Adjacency matrix and trackster index
    std::vector<float> A;
    std::vector<float> trackster_index;
    
    // Add self-loops
    /*
    for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
        if (i == j){
          A.push_back(1.);
          trackster_index.push_back(1.);
        }
        else {
          A.push_back(0.);
          trackster_index.push_back(0.);
        }
      }
    }
    */
    
    // Add undirectional edges
    for (int i=0; i<numEdges; i++){
      // for source
      for (int j=0; j<N; j++){
        if (j == edges_src[i]){A.push_back(1.);} //static_cast<float>(1.0)
        else {A.push_back(0.);}
        if (j == edges_dst[i]){trackster_index.push_back(1.);}
        else {trackster_index.push_back(0.);}
      }
      
      // for dst
      for (int j=0; j<N; j++){
        if (j == edges_dst[i]){A.push_back(1.);}
        else {A.push_back(0.);}
        if (j == edges_src[i]){trackster_index.push_back(1.);}
        else {trackster_index.push_back(0.);}
      }
    }
    
    //input_shapes.push_back({1, 2*numEdges+N, N});
    input_shapes.push_back({1, 2*numEdges, N});
    data.emplace_back(A);
    //input_shapes.push_back({1, 2*numEdges+N, N});
    input_shapes.push_back({1, 2*numEdges, N});
    data.emplace_back(trackster_index);
    
    std::cout << "Adj size: " << A.size() << std::endl;
  
    // Get network output
    std::vector<float> edge_predictions =  onnxRuntime_->run(input_names, data, input_shapes)[0];
  
    for (int i = 0; i < static_cast<int>(edge_predictions.size()); i++) {
      std::cout << "Network output for edge " << data[1][i] << "-" << data[1][numEdges + i]
                << " is: " << edge_predictions[i] << std::endl;
    }
  
    // Create a graph
    Graph g;
    const auto classification_threshold = 0.7;
  
    // Self-loop for not connected nodes
    for (int i = 0; i < N; i++){
      g.addEdge(i, i);
    }
    // Building a predicted graph
    for (int i = 0; i < numEdges; i++) {
      if (edge_predictions[i] >= classification_threshold) {
        auto src = data[1][i];
        auto dst = data[1][numEdges + i];
        // Make undirectional
        g.addEdge(src, dst);
        g.addEdge(dst, src);
      }
    }
  
    std::cout << "Following is Depth First Traversal" << std::endl;
    std::cout << "Connected components are: " << std::endl;
    g.DFS();
  
    int i = 0;
    
    linkedTracksterIdToInputTracksterId.resize(g.connected_components.size());
  
    for (auto &component : g.connected_components) {
    
      std::vector<unsigned int> linkedTracksters;
      Trackster outTrackster;
      
      for (auto &trackster_id : component) {
        std::cout << "Component " << i << ": trackster id " << trackster_id << std::endl;
        linkedTracksterIdToInputTracksterId[i].push_back(trackster_id);
        outTrackster.mergeTracksters(input.tracksters[trackster_id]);
      }
      i++;
      linkedTracksters.push_back(resultTracksters.size());
      resultTracksters.push_back(outTrackster);
      // Store the linked tracksters
      outputSuperclusters.push_back(linkedTracksters);
    }    
    std::cout << "Done" << std::endl;
}

void TracksterLinkingbyGNN::fillPSetDescription(edm::ParameterSetDescription &desc) {
  TracksterLinkingAlgoBase::fillPSetDescription(desc); // adds algo_verbosity
  desc.add<edm::FileInPath>("onnxModelPath", edm::FileInPath("RecoHGCal/TICL/data/tf_models/gnn.onnx"))
    ->setComment("Path to GNN (as ONNX model)");
  desc.add<std::string>("nnVersion", "gnn-v1")
    ->setComment("GNN version tag.");
  desc.add<double>("nnWorkingPoint", 0.91)
     ->setComment("Working point of GNN (in [0, 1]). DNN score above WP will attempt to supercluster.");
  desc.add<double>("deltaEtaWindow", 0.1)
     ->setComment("Size of delta eta window to consider for superclustering. Seed-candidate pairs outside this window are not considered for DNN inference.");
  desc.add<double>("deltaPhiWindow", 0.5)
     ->setComment("Size of delta phi window to consider for superclustering. Seed-candidate pairs outside this window are not considered for DNN inference.");
}