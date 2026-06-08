#include "rng.h"
#include "postTS.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace tsPost {

void checkTsk(const int status, const char *context) {
  if (status < 0) {
    throw std::runtime_error(std::string(context) + ": " +
                             std::string(tsk_strerror(status)));
  }
}

void expandInbredSamplesInPlace(tsk_table_collection_t *tables,
                                const unsigned int ploidy) {
  if (tables == nullptr || ploidy <= 1) {
    return;
  }
  const double sequenceLength = tables->sequence_length;
  if (!(sequenceLength > 0.0)) {
    return;
  }
  
  std::vector<tsk_id_t> originalSamples;
  originalSamples.reserve(tables->nodes.num_rows);
  for (tsk_id_t nodeId = 0;
       nodeId < static_cast<tsk_id_t>(tables->nodes.num_rows);
       ++nodeId) {
    if ((tables->nodes.flags[nodeId] & TSK_NODE_IS_SAMPLE) != 0) {
      originalSamples.push_back(nodeId);
    }
  }
  if (originalSamples.empty()) {
    return;
  }
  
  const double smallestPositive = std::nextafter(0.0, 1.0);
  std::vector<double> minParentTimeByNode(
    tables->nodes.num_rows,
    std::numeric_limits<double>::infinity());
  for (tsk_size_t edgeId = 0; edgeId < tables->edges.num_rows; ++edgeId) {
    const tsk_id_t childId = tables->edges.child[edgeId];
    const tsk_id_t parentId = tables->edges.parent[edgeId];
    if (childId == TSK_NULL || parentId == TSK_NULL ||
        childId < 0 ||
        static_cast<tsk_size_t>(childId) >= tables->nodes.num_rows) {
      continue;
    }
    const double parentTime = tables->nodes.time[parentId];
    if (std::isfinite(parentTime) &&
        parentTime < minParentTimeByNode[childId]) {
      minParentTimeByNode[childId] = parentTime;
    }
  }
  
  for (const tsk_id_t sampleNodeId : originalSamples) {
    const double minParentTime = minParentTimeByNode[sampleNodeId];
    double internalTime = tables->nodes.time[sampleNodeId];
    if (!(internalTime > 0.0)) {
      internalTime = 1e-12;
      if (std::isfinite(minParentTime) && minParentTime > 0.0) {
        internalTime = std::min(internalTime, 0.5 * minParentTime);
      }
      if (!(internalTime > 0.0)) {
        internalTime = smallestPositive;
      }
      if (std::isfinite(minParentTime) && !(internalTime < minParentTime)) {
        internalTime = std::nextafter(minParentTime, 0.0);
      }
      if (!(internalTime > 0.0)) {
        internalTime = smallestPositive;
      }
    }
    
    tables->nodes.time[sampleNodeId] = internalTime;
    tables->nodes.flags[sampleNodeId]
      &= ~static_cast<tsk_flags_t>(TSK_NODE_IS_SAMPLE);
    const tsk_id_t population = tables->nodes.population[sampleNodeId];
    const tsk_id_t individual = tables->nodes.individual[sampleNodeId];
    tables->nodes.individual[sampleNodeId] = TSK_NULL;
    
    for (tsk_size_t mutationId = 0;
         mutationId < tables->mutations.num_rows;
         ++mutationId) {
      if (tables->mutations.node[mutationId] != sampleNodeId) {
        continue;
      }
      const double mutationTime = tables->mutations.time[mutationId];
      if (!tsk_is_unknown_time(mutationTime) && mutationTime < internalTime) {
        tables->mutations.time[mutationId] = internalTime;
      }
    }
    
    for (unsigned int copy = 0; copy < ploidy; ++copy) {
      const tsk_id_t childNodeId = tsk_node_table_add_row(&tables->nodes,
                                                          TSK_NODE_IS_SAMPLE,
                                                          0.0,
                                                          population,
                                                          individual,
                                                          nullptr,
                                                          0);
      checkTsk(static_cast<int>(childNodeId),
               "Failed to add duplicated inbred sample node");
      const tsk_id_t edgeId = tsk_edge_table_add_row(&tables->edges,
                                                     0.0,
                                                     sequenceLength,
                                                     sampleNodeId,
                                                     childNodeId,
                                                     nullptr,
                                                     0);
      checkTsk(static_cast<int>(edgeId),
               "Failed to add duplicated inbred sample edge");
    }
  }
}

void mutateTablesInPlace(tsk_table_collection_t *tables,
                         const double theta,
                         const uint64_t seed) {
  if (!(theta > 0.0) || !std::isfinite(theta)) {
    return;
  }
  if (tables == nullptr) {
    throw std::runtime_error("Table collection pointer is null");
  }
  if (!(tables->sequence_length > 0.0)) {
    throw std::runtime_error("Table collection has invalid sequence_length");
  }
  
  dqrng::rng64_t rng = alphasimrRng::createRng(seed);
  static const char ancestralState[] = "0";
  static const char derivedState[] = "1";
  const double sequenceLength = tables->sequence_length;
  
  for (tsk_size_t edgeId = 0; edgeId < tables->edges.num_rows; ++edgeId) {
    const tsk_id_t parent = tables->edges.parent[edgeId];
    const tsk_id_t child = tables->edges.child[edgeId];
    if (parent == TSK_NULL || child == TSK_NULL ||
        parent < 0 || child < 0) {
      continue;
    }
    
    const double left = tables->edges.left[edgeId];
    const double right = tables->edges.right[edgeId];
    const double span = right - left;
    if (!(span > 0.0)) {
      continue;
    }
    
    const double parentTime = tables->nodes.time[parent];
    const double childTime = tables->nodes.time[child];
    const double branch = parentTime - childTime;
    if (!(branch > 0.0) || !std::isfinite(branch)) {
      continue;
    }
    
    const double spanFraction = span / sequenceLength;
    const double lambda = theta * spanFraction * branch;
    if (!(lambda > 0.0) || !std::isfinite(lambda)) {
      continue;
    }
    
    const arma::uword nMut = alphasimrRng::samplePoisson(lambda, *rng);
    for (arma::uword i = 0; i < nMut; ++i) {
      double position = left + alphasimrRng::runif(*rng) * span;
      if (position >= sequenceLength) {
        position = std::nextafter(sequenceLength, 0.0);
      }
      const tsk_id_t siteId = tsk_site_table_add_row(&tables->sites,
                                                     position,
                                                     ancestralState,
                                                     1,
                                                     nullptr,
                                                     0);
      checkTsk(static_cast<int>(siteId), "Failed to add site row");
      
      double mutationTime = childTime + alphasimrRng::runif(*rng) * branch;
      if (!(mutationTime > childTime)) {
        mutationTime = std::nextafter(childTime, parentTime);
      }
      if (!(mutationTime < parentTime)) {
        mutationTime = std::nextafter(parentTime, childTime);
      }
      
      const tsk_id_t mutationId = tsk_mutation_table_add_row(&tables->mutations,
                                                             siteId,
                                                             child,
                                                             TSK_NULL,
                                                             mutationTime,
                                                             derivedState,
                                                             1,
                                                             nullptr,
                                                             0);
      checkTsk(static_cast<int>(mutationId), "Failed to add mutation row");
    }
  }
  
  (void)tsk_table_collection_drop_index(tables, 0);
  checkTsk(tsk_table_collection_sort(tables, nullptr, 0),
           "Failed to sort table collection after post-TS mutation");
  checkTsk(tsk_table_collection_build_index(tables, 0),
           "Failed to build index after post-TS mutation");
}

} // namespace tsPost
