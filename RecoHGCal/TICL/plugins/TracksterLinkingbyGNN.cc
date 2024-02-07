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

std::vector<float> calculate_dist(const std::vector<std::vector<float>> &lcs_i, const std::vector<std::vector<float>> &lcs_j)
{
    std::vector<float> dists;

    for (const auto &point_i : lcs_i)
    {
        float x_i = point_i[0];
        float y_i = point_i[1];
        float z_i = point_i[2];

        for (const auto &point_j : lcs_j)
        {
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
class Graph
{
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

void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w); // Add w to v�s list.
}

void Graph::DFSUtil(int v)
{
    // Mark the current node as visited and print it
    visited[v] = true;

    // Recursion for all the vertices adjacent to this vertex
    std::list<int>::iterator i;
    for (auto i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[*i])
        {
            connected_components.back().push_back(*i);
            // std::cout << "Pushed back " << *i << std::endl;
            DFSUtil(*i);
        }
}

// The function to do DFS traversal. It uses recursive DFSUtil()
void Graph::DFS()
{
    // Call the recursive helper function to print DFS
    // traversal starting from all vertices one by one
    for (auto i : adj)
        if (visited[i.first] == false)
        {
            // std::cout << "Emplaced back: " << i.first << std::endl;
            connected_components.emplace_back(1, i.first); // {i.first}
            // std::cout << "Starting DFS from node: " << i.first << std::endl;
            DFSUtil(i.first);
        }
}

/*
std::unique_ptr<TICLGraph> produce_ticl_graph(const MultiVectorManager<Trackster> &tracksters)
{

    TICLLayerTile tracksterTilePos;
    TICLLayerTile tracksterTileNeg;

    for (size_t id_t = 0; id_t < tracksters.size(); ++id_t)
    {
        auto t = tracksters[id_t];
        if (t.barycenter().eta() > 0.)
        {
            tracksterTilePos.fill(t.barycenter().eta(), t.barycenter().phi(), id_t);
        }
        else if (t.barycenter().eta() < 0.)
        {
            tracksterTileNeg.fill(t.barycenter().eta(), t.barycenter().phi(), id_t);
        }
    }

    std::vector<Node> allNodes;

    for (size_t id_t = 0; id_t < tracksters.size(); ++id_t)
    {
        auto t = tracksters[id_t];

        Node tNode(id_t);

        auto bary = t.barycenter();
        double del = 0.1;

        double eta_min = std::max(abs(bary.eta()) - del, (double)TileConstants::minEta);
        double eta_max = std::min(abs(bary.eta()) + del, (double)TileConstants::maxEta);

        if (bary.eta() > 0.)
        {
            std::array<int, 4> search_box =
                tracksterTilePos.searchBoxEtaPhi(eta_min, eta_max, bary.phi() - del, bary.phi() + del);
            if (search_box[2] > search_box[3])
            {
                search_box[3] += TileConstants::nPhiBins;
            }

            for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i)
            {
                for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i)
                {
                    auto &neighbours = tracksterTilePos[tracksterTilePos.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
                    for (auto n : neighbours)
                    {
                        if (tracksters[n].barycenter().z() < bary.z())
                        {
                            tNode.addInner(n);
                        }
                        else if (tracksters[n].barycenter().z() > bary.z())
                        {
                            tNode.addOuter(n);
                        }
                    }
                }
            }
        }

        else if (bary.eta() < 0.)
        {
            std::array<int, 4> search_box =
                tracksterTileNeg.searchBoxEtaPhi(eta_min, eta_max, bary.phi() - del, bary.phi() + del);
            if (search_box[2] > search_box[3])
            {
                search_box[3] += TileConstants::nPhiBins;
            }

            for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i)
            {
                for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i)
                {
                    auto &neighbours = tracksterTileNeg[tracksterTileNeg.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
                    for (auto n : neighbours)
                    {
                        if (abs(tracksters[n].barycenter().z()) < abs(bary.z()))
                        {
                            tNode.addInner(n);
                        }
                        else if (abs(tracksters[n].barycenter().z()) > abs(bary.z()))
                        {
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
*/

std::unique_ptr<TICLGraph> produce_ticl_graph(const MultiVectorManager<Trackster> &trackstersclue3d) {

  //std::vector<Trackster> trackstersclue3d_sorted(trackstersclue3d);
  //std::sort(trackstersclue3d_sorted.begin(), trackstersclue3d_sorted.end(), [](Trackster& t1, Trackster& t2){return t1.barycenter().z() < t2.barycenter().z();});

  TICLLayerTile tracksterTilePos;
  TICLLayerTile tracksterTileNeg;

  for (size_t id_t = 0; id_t < trackstersclue3d.size(); ++id_t) {
    auto t = trackstersclue3d[id_t];
    if (t.barycenter().eta() > 0.) {
      tracksterTilePos.fill(t.barycenter().eta(), t.barycenter().phi(), id_t);
    } else if (t.barycenter().eta() < 0.) {
      tracksterTileNeg.fill(t.barycenter().eta(), t.barycenter().phi(), id_t);
    }
  }

  std::vector<Node> allNodes;

  for (size_t id_t = 0; id_t < trackstersclue3d.size(); ++id_t) {
    auto t = trackstersclue3d[id_t];

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
            if (trackstersclue3d[n].barycenter().z() < bary.z()) {
              tNode.addInner(n);
            } else if (trackstersclue3d[n].barycenter().z() > bary.z()) {
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
            if (abs(trackstersclue3d[n].barycenter().z()) < abs(bary.z())) {
              tNode.addInner(n);
            } else if (abs(trackstersclue3d[n].barycenter().z()) > abs(bary.z())) {
              tNode.addOuter(n);
            }
          }
        }
      }
    }
    allNodes.push_back(tNode);
  }
  auto resultGraph = std::make_unique<TICLGraph>(allNodes);

  return resultGraph;
}

TracksterLinkingbyGNN::TracksterLinkingbyGNN(const edm::ParameterSet &ps, edm::ConsumesCollector iC, cms::Ort::ONNXRuntime const *onnxRuntime)
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
void TracksterLinkingbyGNN::linkTracksters(const Inputs &input, std::vector<Trackster> &resultTracksters,
                                           std::vector<std::vector<unsigned int>> &outputSuperclusters,
                                           std::vector<std::vector<unsigned int>> &linkedTracksterIdToInputTracksterId)
{

    std::cout << "Linking Algo by GNN" << std::endl;

    // Create trackster graph
    auto TICL_graph_result = produce_ticl_graph(input.tracksters);

    // Access the TICLGraph object through the unique pointer
    TICLGraph *ticlGraph = TICL_graph_result.get();

    // GNN input names
    // const std::vector<std::string> input_names = {"features", "edge_index", "SRC", "DST", "edge_features"};
    const std::vector<std::string> input_names = {"features", "edge_index", "edge_features"};
    // Array of data to be filled as a network input. Should be a float array of flattened values.
    cms::Ort::FloatArrays data;
    // Network input shapes.
    std::vector<std::vector<int64_t>> input_shapes;
    // Network feature shape.
    const auto shapeFeatures = 28;
    const auto num_edge_features = 7;

    // DEBUG mode if 1, otherwise STANDARD mode. Debugging mode prints out the event details.
    const auto DEBUG = 1;

    auto const &tracksters = input.tracksters;
    auto const &layerClusters = input.layerClusters;
    long int N = input.tracksters.size();

    if (DEBUG == 1)
    {
        // Print out info about tracksters
        std::cout << "Number of tracksters in event: " << N << std::endl;
    }

    if (N < 2)
    {
        // do not run the network - return the original tracksters
        if (DEBUG == 1)
        {
            std::cout << "Number of tracksters less than 2 - no linking is done." << std::endl;
        }
        linkedTracksterIdToInputTracksterId.resize(N);
        std::vector<unsigned int> linkedTracksters;

        for (int trackster_id = 0; trackster_id < N; trackster_id++)
        {
            linkedTracksterIdToInputTracksterId[trackster_id].push_back(trackster_id);
            linkedTracksters.push_back(resultTracksters.size());
            resultTracksters.push_back(input.tracksters[trackster_id]);
            outputSuperclusters.push_back(linkedTracksters);
        }
        return;
    }

    std::vector<float> features;
    std::vector<float> edge_features;

    for (unsigned i = 0; i < N; ++i)
    {

        const auto ts = tracksters[i];
        const Vector &barycenter = ts.barycenter();
        const Vector eigenvector0 = ts.eigenvectors(0);
        const std::array<float, 3> eigenvalues = ts.eigenvalues();
        const std::array<float, 3> sigmasPCA = ts.sigmasPCA();

        if (DEBUG == 1)
        {
            std::cout << "Trackster " << i << "--------------------" << std::endl;
            std::cout << "Barycenter X: " << barycenter.x() << std::endl;
            std::cout << "Barycenter Y: " << barycenter.y() << std::endl;
            std::cout << "Barycenter Z: " << barycenter.z() << std::endl;
            
            std::cout << "Barycenter Y: " << barycenter.eta() << std::endl;
            std::cout << "Barycenter Z: " << barycenter.phi() << std::endl;

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
        }

        // Get representative points of the trackster

        const std::vector<unsigned int> &vertices_indices = ts.vertices();

        const auto max_z_lc_it = std::max_element(
            vertices_indices.begin(),
            vertices_indices.end(),
            [&layerClusters](const int &a, const int &b)
            {
                return layerClusters[a].z() > layerClusters[b].z();
            });

        const auto min_z_lc_it = std::min_element(
            vertices_indices.begin(),
            vertices_indices.end(),
            [&layerClusters](const int &a, const int &b)
            {
                return layerClusters[a].z() > layerClusters[b].z();
            });

        int num_hits = 0;
        std::set<float> lc_z_unique;
        for (auto &lc_idx : ts.vertices())
        {
            auto lc = layerClusters[lc_idx];
            num_hits += lc.size();
            lc_z_unique.insert(lc.z());
        }
        auto length = lc_z_unique.size();

        const reco::CaloCluster &min_z_lc = layerClusters[*min_z_lc_it];
        const reco::CaloCluster &max_z_lc = layerClusters[*max_z_lc_it];

        // size and energy
        if (DEBUG == 1)
        {
            std::cout << "Num LCs: " << ts.vertices().size() << std::endl;
            std::cout << "Num hits: " << ts.vertices().size() << std::endl;
            std::cout << "Raw Energy: " << ts.raw_energy() << std::endl;
            std::cout << "Raw EM Energy: " << ts.raw_em_energy() << std::endl;
            std::cout << "Length z layers: " << length << std::endl;
            std::cout << "Min Z: " << min_z_lc.z() << std::endl;
            std::cout << "Max Z: " << max_z_lc.z() << std::endl;
            std::cout << "Time: " << ts.time() << std::endl;

            std::cout << "prob photon: " << ts.id_probabilities(0) << std::endl;
            std::cout << "prob electrin: " << ts.id_probabilities(1) << std::endl;
            std::cout << "prob muon: " << ts.id_probabilities(2) << std::endl;
            std::cout << "prob neutral pion: " << ts.id_probabilities(3) << std::endl;
            std::cout << "prob charged hadr: " << ts.id_probabilities(4) << std::endl;
            std::cout << "prob neutral hadr: " << ts.id_probabilities(5) << std::endl;
        }

        // Features
        features.push_back(barycenter.x()); // 0:  barycenter x
        features.push_back(barycenter.y()); // 1:  barycenter y
        features.push_back(barycenter.z()); // 2:  barycenter z

        features.push_back(barycenter.eta()); // 3: trackster_barycenter_eta
        features.push_back(barycenter.phi()); // 4: trackster_barycenter_phi

        features.push_back(eigenvector0.x()); // 5: eVector0_x
        features.push_back(eigenvector0.y()); // 6: eVector0_y
        features.push_back(eigenvector0.z()); // 7: eVector0_z

        features.push_back(eigenvalues[0]); // 8: EV1
        features.push_back(eigenvalues[1]); // 9: EV2
        features.push_back(eigenvalues[2]); // 10: EV3

        features.push_back(sigmasPCA[0]); // 11: sigmaPCA1
        features.push_back(sigmasPCA[1]); // 12: sigmaPCA2
        features.push_back(sigmasPCA[2]); // 13: sigmaPCA3

        features.push_back(ts.vertices().size()); // 14: number of LCs per trackster
        features.push_back(num_hits);             // 15: number of hits per trackster

        features.push_back(ts.raw_energy());    // 16: raw_energy
        features.push_back(ts.raw_em_energy()); // 17: raw_em_energy

        features.push_back(ts.id_probabilities(0)); // 18: photon probability
        features.push_back(ts.id_probabilities(1)); // 19: electron probability
        features.push_back(ts.id_probabilities(2)); // 20: muon probability
        features.push_back(ts.id_probabilities(3)); // 21: neutral_pion probability
        features.push_back(ts.id_probabilities(4)); // 22: charged_hadron probability
        features.push_back(ts.id_probabilities(5)); // 23: neutral_hadron probability

        features.push_back(min_z_lc.z()); // 24: z_min (minimum z of constituent LCs)
        features.push_back(max_z_lc.z()); // 25: z_max (maximum z of constituent LCs)
        features.push_back(length);       // 26: trackster length from minimum to maximum element in terms of layers

        features.push_back(ts.time()); // 27: time

        if (DEBUG == 1)
        {
            std::cout << "--------------------" << std::endl;
        }
    }

    input_shapes.push_back({1, N, shapeFeatures});
    data.emplace_back(features);

    std::vector<float> edges_src;
    std::vector<float> edges_dst;
    for (int i = 0; i < N; i++)
    {
        const auto &ts_i = tracksters[i];

        std::vector<float> vertex_energy_i;
        std::vector<float> vertex_i_x;
        std::vector<float> vertex_i_y;
        std::vector<float> vertex_i_z;

        const Vector &eigenvector0_i = ts_i.eigenvectors(0);

        for (auto &lc_idx : ts_i.vertices())
        {
            auto lc = layerClusters[lc_idx];
            vertex_energy_i.push_back(lc.energy());
            vertex_i_x.push_back(lc.x());
            vertex_i_y.push_back(lc.y());
            vertex_i_z.push_back(lc.z());
        }

        for (auto &i_neighbour : ticlGraph->getNode(i).getInner())
        {

            const auto &ts_j = tracksters[i_neighbour];

            std::vector<float> vertex_energy_j;
            std::vector<float> vertex_j_x;
            std::vector<float> vertex_j_y;
            std::vector<float> vertex_j_z;

            const Vector &eigenvector0_j = ts_j.eigenvectors(0);

            for (auto &lc_idx : ts_j.vertices())
            {
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
            for (int k = 0; k < (int)vertex_energy_i.size(); k++)
            {
                float x = vertex_i_x[k];
                float y = vertex_i_y[k];
                float z = vertex_i_z[k];
                float e = vertex_energy_i[k];
                if (e >= e_thr_i)
                    lcs_i.push_back({x, y, z});
            }

            std::vector<std::vector<float>> lcs_j;
            for (int k = 0; k < (int)vertex_energy_j.size(); k++)
            {
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

            edge_features.push_back(delta_en);  // 0 - energy_difference
            edge_features.push_back(delta_z);   // 1 - delta_z
            edge_features.push_back(min_dist);  // 2 - min_dist
            edge_features.push_back(max_dist);  // 3 - max_dist
            edge_features.push_back(delta_r);   // 4 - delta_r
            edge_features.push_back(prod);      // 5 - spatial_comp
            edge_features.push_back(time_comp); // 6 - time_comp

            if (DEBUG == 1)
            {
                std::cout << "en_diff: " << delta_en << ", delta_z: " << delta_z
                          << ", min_dist: " << min_dist << ", max_dist: " << max_dist
                          << ", delta_r: " << delta_r
                          << ", spatial_comp: " << spatial_comp << ", time_comp: " << time_comp << std::endl;
            }
        }
    }

    auto numEdges = static_cast<int>(edges_src.size());

    if (numEdges < 1)
    {
        if (DEBUG == 1)
        {
            std::cout << "No edges for the event - no linking is done." << std::endl;
        }
        // do not run the network - return the original tracksters
        linkedTracksterIdToInputTracksterId.resize(N);
        std::vector<unsigned int> linkedTracksters;

        for (int trackster_id = 0; trackster_id < N; trackster_id++)
        {
            linkedTracksterIdToInputTracksterId[trackster_id].push_back(trackster_id);
            linkedTracksters.push_back(resultTracksters.size());
            resultTracksters.push_back(input.tracksters[trackster_id]);
            outputSuperclusters.push_back(linkedTracksters);
        }
        return;
    }

    input_shapes.push_back({1, 2, numEdges});
    // std::cout << "Num edges: " << numEdges << std::endl;

    data.emplace_back(edges_src);
    for (auto &dst : edges_dst)
    {
        data.back().push_back(dst);
    }

    input_shapes.push_back({1, numEdges, num_edge_features});
    data.emplace_back(edge_features);
    
    // DEBUGGING!

    std::cout << "Input Names: " << input_names.size() << std::endl;
    for (const auto& name : input_names) {
        std::cout << name << std::endl;
    }
    std::cout << "Input Shapes: " << input_shapes.size() <<  std::endl;
    for (const auto& shape : input_shapes) {
        std::cout << "Shape: [";
        for (size_t i = 0; i < shape.size(); ++i) {
            std::cout << shape[i];
            if (i < shape.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
    }
    std::cout << "Input Data: " << data.size() << std::endl;
    for (const auto& floatArray : data) {
        std::cout << "FloatArray: " << floatArray.size() << std::endl;
    }

    std::vector<float> edge_predictions = onnxRuntime_->run(input_names, data, input_shapes)[0];

    if (DEBUG == 1)
    {
        std::cout << "Network output shape is " << edge_predictions.size() << std::endl;
        for (int i = 0; i < static_cast<int>(edge_predictions.size()); i++)
        {
            std::cout << "Network output for edge " << data[1][i] << "-" << data[1][numEdges + i]
                        << " is: " << edge_predictions[i] << std::endl;
        }
    }

    // Create a graph
    Graph g;
    const auto classification_threshold = 0.94;

    // Self-loop for not connected nodes.
    for (int i = 0; i < N; i++)
    {
        g.addEdge(i, i);
    }
    // Building a predicted graph.
    for (int i = 0; i < numEdges; i++)
    {
        if (edge_predictions[i] >= classification_threshold)
        {
            auto src = data[1][i];
            auto dst = data[1][numEdges + i];
            // Make undirectional
            g.addEdge(src, dst);
            g.addEdge(dst, src);
        }
    }

    std::cout << "Following Depth First Traversal" << std::endl;
    std::cout << "Connected components are: " << std::endl;
    g.DFS();

    int id = 0;

    linkedTracksterIdToInputTracksterId.resize(g.connected_components.size());

    for (auto &component : g.connected_components)
    {
        std::vector<unsigned int> linkedTracksters;
        Trackster outTrackster;

        for (auto &trackster_id : component)
        {
            std::cout << "Component " << id << ": trackster id " << trackster_id << std::endl;
            linkedTracksterIdToInputTracksterId[id].push_back(trackster_id);
            outTrackster.mergeTracksters(input.tracksters[trackster_id]);
        }
        id++;
        linkedTracksters.push_back(resultTracksters.size());
        resultTracksters.push_back(outTrackster);
        // Store the linked tracksters
        outputSuperclusters.push_back(linkedTracksters);
    }
    std::cout << "Done." << std::endl;
}

void TracksterLinkingbyGNN::fillPSetDescription(edm::ParameterSetDescription &desc)
{
    TracksterLinkingAlgoBase::fillPSetDescription(desc); // adds algo_verbosity
    desc.add<edm::FileInPath>("onnxModelPath", edm::FileInPath("RecoHGCal/TICL/data/tf_models/gnn_linking_model.onnx"))
        ->setComment("Path to GNN (as ONNX model)"); // gnn_linking_model.onnx, mlp_linking_model.onnx or mlp_linking_model_with_edge_features.onnx
    desc.add<std::string>("nnVersion", "gnn-v1")
        ->setComment("GNN version tag.");
    desc.add<double>("nnWorkingPoint", 0.7)
        ->setComment("Working point of GNN (in [0, 1]). DNN score above WP will attempt to supercluster.");
    desc.add<double>("deltaEtaWindow", 0.1)
        ->setComment("Size of delta eta window to consider for superclustering. Seed-candidate pairs outside this window are not considered for DNN inference.");
    desc.add<double>("deltaPhiWindow", 0.5)
        ->setComment("Size of delta phi window to consider for superclustering. Seed-candidate pairs outside this window are not considered for DNN inference.");
}