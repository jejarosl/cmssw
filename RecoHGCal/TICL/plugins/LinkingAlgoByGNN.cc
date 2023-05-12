#include <cmath>
#include <string>
#include "RecoHGCal/TICL/plugins/LinkingAlgoByGNN.h"

#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"

#include <bits/stdc++.h>

#define HGCAL_LEN 300.0

using namespace std;
using namespace ticl;
using namespace cms::Ort;

std::vector<float> calculate_dist(const std::vector<std::vector<float>>& lcs_i, const std::vector<std::vector<float>>& lcs_j) {
    std::vector<float> dists;
    
    for (const auto& point_i : lcs_i) {
        float x_i = point_i[0];
        float y_i = point_i[1];
        float z_i = point_i[2];
        
        for (const auto& point_j : lcs_j) {
            float x_j = point_j[0];
            float y_j = point_j[1];
            float z_j = point_j[2];
            
            float distance = std::sqrt(std::pow((x_i - x_j), 2) + std::pow((y_i - y_j), 2) + std::pow((z_i - z_j), 2));
            dists.push_back(distance);
        }
    }
    
    return dists;
}

// Graph search of connected components, as a post-processing of the network output
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
  adj[v].push_back(w);  // Add w to v’s list.
}

void Graph::DFSUtil(int v) {
  // Mark the current node as visited and print it
  visited[v] = true;

  // Recursion for all the vertices adjacent to this vertex
  list<int>::iterator i;
  for (auto i = adj[v].begin(); i != adj[v].end(); ++i)
    if (!visited[*i]) {
      connected_components.back().push_back(*i);
      //std::cout << "Pushed back " << *i << std::endl;
      DFSUtil(*i);
    }
}

// The function to do DFS traversal. It uses recursive DFSUtil()
void Graph::DFS() {
  // Call the recursive helper function to print DFS
  // traversal starting from all vertices one by one
  for (auto i : adj)
    if (visited[i.first] == false) {
      //std::cout << "Emplaced back: " << i.first << std::endl;
      connected_components.emplace_back(1, i.first);  // {i.first}
      //std::cout << "Starting DFS from node: " << i.first << std::endl;
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
                                                     std::vector<TICLCandidate> &resultLinked,
                                                     std::vector<TICLCandidate> &chargedHadronsFromTk,
                                                     std::vector<double>& prop_tracks_x,
                                                     std::vector<double>& prop_tracks_y,
                                                     std::vector<double>& prop_tracks_z,
                                                     std::vector<double>& prop_tracks_eta,
                                                     std::vector<double>& prop_tracks_phi,
                                                     std::vector<double>& prop_tracks_px,
                                                     std::vector<double>& prop_tracks_py,
                                                     std::vector<double>& prop_tracks_pz,
                                                     std::vector<bool>& masked_tracks, 
                                                     const TICLGraph &ticlGraph,
                                                     const std::vector<reco::CaloCluster>& layerClusters,
                                                     const ONNXRuntime *cache)  {
  std::cout << "Linking Algo by GNN" << std::endl;
  // GNN input names
  const std::vector<std::string> input_names = {"features", "edge_index", "SRC", "DST", "edge_features"};
  // Array of data to be filled as a network input. Should be a float array of flattened values.
  FloatArrays data;
  // Network input shapes.
  std::vector<std::vector<int64_t>> input_shapes;
  // Network feature shape.
  const auto shapeFeatures = 28;
  const auto num_edge_features = 7;
  
  // TEST mode if 1, otherwise RUN mode. Test mode provides the network with a test input.
  const auto TEST = 0;
  // DEBUG mode if 1, otherwise STANDARD mode. Debugging mode prints out the event details.
  const auto DEBUG = 0;

  
  if (TEST == 1){
  

    std::cout << "Providing GNN with a test input." << std::endl;
  
    // Test input
    const auto N = 10;
    
    input_shapes.push_back({1, N, shapeFeatures});
    input_shapes.push_back({1, 2, 3*N});
    input_shapes.push_back({1, N, N});
    input_shapes.push_back({1, N, N});
    
    data.emplace_back(N*shapeFeatures, 0.1);
    data.emplace_back(3*N, 0);
    
    for (int i=0; i < 3*N; i++){
      data.back().push_back(1);
    }
    
    // Creating Adjacency matrix and trackster index: unity matrices of size N x N
    std::vector<float> A;
    std::vector<float> trackster_index;
    
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
    
    data.emplace_back(A);
    data.emplace_back(trackster_index);
  
  
    std::vector<float> edge_predictions = cache->run(input_names, data, input_shapes)[0];
  
    std::cout << "Network output shape is " << edge_predictions.size() << std::endl;
  
    for (int i = 0; i < static_cast<int>(edge_predictions.size()); i++) {
      std::cout << "Network output for edge " << data[1][i] << "-" << data[1][3*N + i]
                << " is: " << edge_predictions[i] << std::endl;
    }
  }
  
  else{
    const auto &tracks = *tkH;
    const auto &tracksters = *tsH;
  
    auto bFieldProd = bfield_.product();
    const Propagator &prop = (*propagator_);
  
    long int N = tracksters.size();
    
    if (DEBUG == 1){
      // Print out info about tracksters
      std::cout << "Number of tracksters in event: " << N << std::endl;
    }
    
    if (N < 2){
      // do not run the network - return the original tracksters
      // TODO: the same if zero edges
      if (DEBUG == 1){
        std::cout << "Number of tracksters less than 2 - no linking is done." << std::endl;
      }
      std::vector<TICLCandidate> connectedCandidates;
      TICLCandidate tracksterCandidate;
    
      for (int trackster_id = 0; trackster_id < N; trackster_id++) {
        tracksterCandidate.addTrackster(edm::Ptr<Trackster>(tsH, trackster_id));     
      }
      connectedCandidates.push_back(tracksterCandidate);
    
      // The final candidates are passed to `resultLinked`
      resultLinked.insert(std::end(resultLinked), std::begin(connectedCandidates), std::end(connectedCandidates));
      return;
    }

    std::vector<float> features;
    std::vector<float> edge_features;
    // For nearest neighbour finding
    //std::vector<float> barycenters_x;
    //std::vector<float> barycenters_y;
    //std::vector<float> barycenters_z;

    for (unsigned i = 0; i < N; ++i) {
      const auto &ts = tracksters[i];
      
      const Vector &barycenter = ts.barycenter();
      if (DEBUG == 1){
        std::cout << "Trackster " << i << "--------------------" << std::endl;
        std::cout << "Barycenter X: " << barycenter.x() << std::endl;
        std::cout << "Barycenter Y: " << barycenter.y() << std::endl;
        std::cout << "Barycenter Z: " << barycenter.z() << std::endl;
      }
  
      const Vector &eigenvector0 = ts.eigenvectors(0);
      if (DEBUG == 1){
        std::cout << "eVector0 X: " << eigenvector0.x() << std::endl;
        std::cout << "eVector0 Y: " << eigenvector0.y() << std::endl;
        std::cout << "eVector0 Z: " << eigenvector0.z() << std::endl;
      }
  
      const std::array<float, 3> &eigenvalues = ts.eigenvalues();
      if (DEBUG == 1){
        std::cout << "EV1: " << eigenvalues[0] << std::endl;
        std::cout << "EV2: " << eigenvalues[1] << std::endl;
        std::cout << "EV3: " << eigenvalues[2] << std::endl;
      }
  
      const std::array<float, 3> &sigmasPCA = ts.sigmasPCA();
      if (DEBUG == 1){
        std::cout << "sigmaPCA1: " << sigmasPCA[0] << std::endl;
        std::cout << "sigmaPCA2: " << sigmasPCA[1] << std::endl;
        std::cout << "sigmaPCA3: " << sigmasPCA[2] << std::endl;
      }
      
      // Get representative points of the trackster

      const std::vector<unsigned int> &vertices_indices = ts.vertices();
      
      const auto max_z_lc_it = std::max_element(
      vertices_indices.begin(),
      vertices_indices.end(),
      [&layerClusters](const int &a, const int &b) {
        return layerClusters[a].z() > layerClusters[b].z();
        }
      );
  
      const auto min_z_lc_it = std::min_element(
        vertices_indices.begin(),
        vertices_indices.end(),
        [&layerClusters](const int &a, const int &b) {
          return layerClusters[a].z() > layerClusters[b].z();
        }
      );
      
      
      int num_hits = 0;
      std::set<float> lc_z_unique;
      for (auto & lc_idx : ts.vertices()){
        auto lc = layerClusters[lc_idx];
        num_hits += lc.size();
        lc_z_unique.insert(lc.z());
      }
      auto length = lc_z_unique.size();
      
      
      const reco::CaloCluster &min_z_lc = layerClusters[*min_z_lc_it];
      const reco::CaloCluster &max_z_lc = layerClusters[*max_z_lc_it];
  
      // size and energy
      if (DEBUG == 1){
        std::cout << "Size: " << ts.vertices().size() << std::endl;
        std::cout << "Raw Energy: " << ts.raw_energy() << std::endl;
        std::cout << "Raw EM Energy: " << ts.raw_em_energy() << std::endl;
        std::cout << "Length: " << length << std::endl;
        std::cout << "Min Z: " << min_z_lc.z() << std::endl;
        std::cout << "Max Z: " << max_z_lc.z() << std::endl;
        std::cout << "Time: " << ts.time() << std::endl;
      }
      

      // Features
      features.push_back(barycenter.x());         // 0:  barycenter x
      features.push_back(barycenter.y());         // 1:  barycenter y
      features.push_back(barycenter.z());         // 2:  barycenter z
      
      features.push_back(barycenter.eta());       // 3: trackster_barycenter_eta
      features.push_back(barycenter.phi());       // 4: trackster_barycenter_phi
      
      features.push_back(eigenvector0.x());       // 5: eVector0_x
      features.push_back(eigenvector0.y());       // 6: eVector0_y
      features.push_back(eigenvector0.z());       // 7: eVector0_z
      
      features.push_back(eigenvalues[0]);         // 8: EV1
      features.push_back(eigenvalues[1]);         // 9: EV2
      features.push_back(eigenvalues[2]);         // 10: EV3
      
      features.push_back(sigmasPCA[0]);           // 11: sigmaPCA1
      features.push_back(sigmasPCA[1]);           // 12: sigmaPCA2
      features.push_back(sigmasPCA[2]);           // 13: sigmaPCA3
      
      features.push_back(ts.vertices().size());   // 14: number of LCs per trackster
      features.push_back(num_hits);               // 15: number of hits per trackster
      
      features.push_back(ts.raw_energy());        // 16: raw_energy
      features.push_back(ts.raw_em_energy());     // 17: raw_em_energy
      
      features.push_back(ts.id_probabilities(0)); // 18: photon probability
      features.push_back(ts.id_probabilities(1)); // 19: electron probability
      features.push_back(ts.id_probabilities(2)); // 20: muon probability
      features.push_back(ts.id_probabilities(3)); // 21: neutral_pion probability
      features.push_back(ts.id_probabilities(4)); // 22: charged_hadron probability
      features.push_back(ts.id_probabilities(5)); // 23: neutral_hadron probability
      
      features.push_back(min_z_lc.z());           // 24: z_min (minimum z of constituent LCs)
      features.push_back(max_z_lc.z());           // 25: z_max (maximum z of constituent LCs)
      features.push_back(length);                 // 26: trackster length from minimum to maximum element in terms of layers
      
      features.push_back(ts.time());              // 27: time
            
      //barycenters_x.push_back(ts.barycenter().x());
      //barycenters_y.push_back(ts.barycenter().y());
      //barycenters_z.push_back(ts.barycenter().z());
      
      if (DEBUG == 1){
        std::cout << "--------------------" << std::endl;
      }
    }
  
    input_shapes.push_back({1, N, shapeFeatures});
    data.emplace_back(features);
  
  
    std::vector<float> edges_src;
    std::vector<float> edges_dst;
    for (int i = 0; i < N; i++){
      const auto &ts_i = tracksters[i];
      
      std::vector<float> vertex_energy_i;
      std::vector<float> vertex_i_x;
      std::vector<float> vertex_i_y;
      std::vector<float> vertex_i_z;
      
      const Vector &eigenvector0_i = ts_i.eigenvectors(0);
      
      for (auto & lc_idx : ts_i.vertices()){
        auto lc = layerClusters[lc_idx];
        vertex_energy_i.push_back(lc.energy());
        vertex_i_x.push_back(lc.x());
        vertex_i_y.push_back(lc.y());
        vertex_i_z.push_back(lc.z()); 
      }
      
      for (auto & i_neighbour : ticlGraph.getNode(i).getInner()){
      
        const auto &ts_j = tracksters[i_neighbour];
        
        std::vector<float> vertex_energy_j;
        std::vector<float> vertex_j_x;
        std::vector<float> vertex_j_y;
        std::vector<float> vertex_j_z;
        
        const Vector &eigenvector0_j = ts_j.eigenvectors(0);
        
        for (auto & lc_idx : ts_j.vertices()){
          auto lc = layerClusters[lc_idx];
          vertex_energy_j.push_back(lc.energy());
          vertex_j_x.push_back(lc.x());
          vertex_j_y.push_back(lc.y());
          vertex_j_z.push_back(lc.z()); 
        }
        
        // Create an edge between the tracksters
        edges_src.push_back(static_cast<float>(i_neighbour));
        edges_dst.push_back(static_cast<float>(i));
        
        
        // TODO: add edge features
        auto delta_en = std::abs(ts_i.raw_energy() - ts_j.raw_energy());
        auto delta_z = std::abs(ts_i.barycenter().z() - ts_j.barycenter().z());
        
        // Find energy thresholds
        float e_thr_i = 0.01 * *std::max_element(vertex_energy_i.begin(), vertex_energy_i.end());
        float e_thr_j = 0.01 * *std::max_element(vertex_energy_j.begin(), vertex_energy_j.end());
        
        // Filter vertices based on energy threshold
        std::vector<std::vector<float>> lcs_i;
        for (int k = 0; k < (int) vertex_energy_i.size(); k++) {
            float x = vertex_i_x[k];
            float y = vertex_i_y[k];
            float z = vertex_i_z[k];
            float e = vertex_energy_i[k];
            if (e >= e_thr_i)
                lcs_i.push_back({x, y, z});
        }
        
        std::vector<std::vector<float>> lcs_j;
        for (int k = 0; k < (int) vertex_energy_j.size(); k++) {
            float x = vertex_j_x[k];
            float y = vertex_j_y[k];
            float z = vertex_j_z[k];
            float e = vertex_energy_j[k];
            if (e >= e_thr_j)
                lcs_j.push_back({x, y, z});
        }
        
        // Calculate distances
        auto dists = calculate_dist(lcs_i, lcs_j);
        float min_dist = *std::min_element(dists.begin(), dists.end());
        float max_dist = *std::max_element(dists.begin(), dists.end());
        
        // Calculate delta_r 
        float delta_r = std::sqrt(std::pow((ts_i.barycenter().x() - ts_j.barycenter().x()), 2) + 
                                   std::pow((ts_i.barycenter().y() - ts_j.barycenter().y()), 2));
        
        // Calculate spatial compatibility
        float prod = eigenvector0_i.x() * eigenvector0_j.x() +
                      eigenvector0_i.y() * eigenvector0_j.y() +
                      eigenvector0_i.z() * eigenvector0_j.z();
        prod = (prod > 1) ? 1 : prod;
        prod = (prod < 0) ? 0 : prod;
        float spatial_comp = std::acos(prod);
        
        // Calculate time compatibility
        float time_comp = 0.0;
        if (ts_i.time() > -99 && ts_j.time() > -99)
            time_comp = std::abs(ts_i.time() - ts_j.time());
        else
            time_comp = -10.0;
          
        
        edge_features.push_back(delta_en);   // 0 - energy_difference
        edge_features.push_back(delta_z);    // 1 - delta_z
        edge_features.push_back(min_dist);   // 2 - min_dist
        edge_features.push_back(max_dist);   // 3 - max_dist
        edge_features.push_back(delta_r);    // 4 - delta_r
        edge_features.push_back(prod);       // 5 - spatial_comp
        edge_features.push_back(time_comp);  // 6 - time_comp
        
        if (DEBUG==1){
          std::cout << "en_diff: " << delta_en << ", delta_z: " << delta_z
                  << ", min_dist: " << min_dist << ", max_dist: " << max_dist
                  << ", delta_r: " << delta_r
                  << ", spatial_comp: " << spatial_comp << ", time_comp: " << time_comp << std::endl;
        }
      }
      
      /* 
      // Nearest Neighbour connection
      if (len(ticlGraph.getNode(i).getInner()) == 0 && len(ticlGraph.getNode(i).getOuter()) == 0){
       
        const auto pos_i = [barycenters_x[i], barycenters_y[i], barycenters_z[i]];
        const auto d_least = 1000;
        for (auto k=0; k< len(barycenters_x; k++){
            if (k == i){
                continue;
            }
            pos_k = [barycenters_x[k], barycenters_y[k], barycenters_z[k]];
            del_pos = pos_k - pos_i;
            d_squared = del_pos[0]**2 + del_pos[1]**2 + del_pos[2]**2;
            if (d_squared < d_least){
                d_least = d_squared;
                i_least = k;
            }
        }
                
        const auto nearest_id = findNearestNeighbour(i, b_x, b_y, b_z)
        edges_src.push_back(static_cast<float>(nearest_id));
        edges_dst.push_back(static_cast<float>(i));
        }
      */
    }
  
    /*
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
    */
  
    auto numEdges = static_cast<int>(edges_src.size());
    
    if (numEdges < 1){
      std::cout << "No edges for the event - no linking is done." << std::endl;
      // do not run the network - return the original tracksters
      std::vector<TICLCandidate> connectedCandidates;
      TICLCandidate tracksterCandidate;
    
      for (int trackster_id = 0; trackster_id < N; trackster_id++) {
        tracksterCandidate.addTrackster(edm::Ptr<Trackster>(tsH, trackster_id));     
      }
      connectedCandidates.push_back(tracksterCandidate);
    
      // The final candidates are passed to `resultLinked`
      resultLinked.insert(std::end(resultLinked), std::begin(connectedCandidates), std::end(connectedCandidates));
      return;
    }
    
    input_shapes.push_back({1, 2, numEdges});
    //std::cout << "Num edges: " << numEdges << std::endl;
  
    data.emplace_back(edges_src);
    for (auto &dst : edges_dst) {
      data.back().push_back(dst);
    }
  
    
    // Creating Adjacency matrix and trackster index
    std::vector<float> A;
    std::vector<float> trackster_index;

    // self-loops
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
    
    for (int i=0; i<numEdges; i++){
      // for source
      for (int j=0; j<N; j++){
        if (j == edges_src[i]){A.push_back(1.);}
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
    
    input_shapes.push_back({1, 2*numEdges + N, N});
    data.emplace_back(A);
    input_shapes.push_back({1, 2*numEdges + N, N});
    data.emplace_back(trackster_index);
    
    input_shapes.push_back({1, numEdges, num_edge_features});
    data.emplace_back(edge_features);
    
    if (DEBUG == 1){
      std::cout << "Adj size: " << A.size() << std::endl;
    }
   
    std::vector<float> edge_predictions = cache->run(input_names, data, input_shapes)[0];
  
    if (DEBUG == 1){
      std::cout << "Network output shape is " << edge_predictions.size() << std::endl;
      for (int i = 0; i < static_cast<int>(edge_predictions.size()); i++) {
        std::cout << "Network output for edge " << data[1][i] << "-" << data[1][numEdges + i]
                  << " is: " << edge_predictions[i] << std::endl;
      }
    }
  
    // Create a graph
    Graph g;
    const auto classification_threshold = 0.9;
  
    // Self-loop for not connected nodes.
    for (int i = 0; i < N; i++){
      g.addEdge(i, i);
    }
    // Building a predicted graph.
    for (int i = 0; i < numEdges; i++) {
      if (edge_predictions[i] >= classification_threshold) {
        auto src = data[1][i];
        auto dst = data[1][numEdges + i];
        // Make undirectional
        g.addEdge(src, dst);
        g.addEdge(dst, src);
      }
    }
  
    //std::cout << "Following Depth First Traversal" << std::endl;
    //std::cout << "Connected components are: " << std::endl;
    g.DFS();
  
    int i = 0;
    std::vector<TICLCandidate> connectedCandidates;
  
    for (auto &component : g.connected_components) {
      TICLCandidate tracksterCandidate;
      for (auto &trackster_id : component) {
        //std::cout << "Component " << i << ": trackster id " << trackster_id << std::endl;
        tracksterCandidate.addTrackster(edm::Ptr<Trackster>(tsH, trackster_id));
      }
      i++;
      connectedCandidates.push_back(tracksterCandidate);
    }
  
    // The final candidates are passed to `resultLinked`
    resultLinked.insert(std::end(resultLinked), std::begin(connectedCandidates), std::end(connectedCandidates));
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

