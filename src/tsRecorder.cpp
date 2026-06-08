#include "simulator.h"
#include "postTS.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {

inline void checkTsk(int status, const char * context) {
  if (status < 0) {
    throw std::runtime_error(std::string(context) + ": " +
      std::string(tsk_strerror(status)));
  }
}

void rescaleTableCollectionTimes(tsk_table_collection_t * tables, double timeScale) {
  if (tables == nullptr || timeScale == 1.0) {
    return;
  }
  if (!(timeScale > 0.0) || !std::isfinite(timeScale)) {
    throw std::runtime_error("timeScale must be a finite positive value");
  }
  
  for (tsk_size_t i = 0; i < tables->nodes.num_rows; ++i) {
    tables->nodes.time[i] *= timeScale;
  }
  for (tsk_size_t i = 0; i < tables->mutations.num_rows; ++i) {
    const double t = tables->mutations.time[i];
    if (!tsk_is_unknown_time(t)) {
      tables->mutations.time[i] = t * timeScale;
    }
  }
  for (tsk_size_t i = 0; i < tables->migrations.num_rows; ++i) {
    tables->migrations.time[i] *= timeScale;
  }
  
  static const char generationUnits[] = "generations";
  checkTsk(tsk_table_collection_set_time_units(
             tables, generationUnits, sizeof(generationUnits) - 1),
           "Failed to set TS time_units");
}

} // namespace

TsRecorder::TsRecorder(double seqLengthBp, TsPositionMode positionMode,
                       bool inbred, unsigned int ploidy):
  pTables(nullptr),
  positionMode(positionMode),
  dSequenceLengthBp(seqLengthBp > 0.0 ? seqLengthBp : 1.0),
  bSamplesPreRegistered(false),
  bSimplified(false),
  bInbred(inbred),
  iPloidy(ploidy > 0 ? ploidy : 1) {
  
  pTables = new tsk_table_collection_t;
  checkTsk(tsk_table_collection_init(pTables, 0), "Failed to initialise tsk tables");
  pTables->sequence_length = positionMode == TsPositionMode::PHYSICAL_BP ?
    dSequenceLengthBp : 1.0;
}

TsRecorder::~TsRecorder() {
  if (pTables != nullptr) {
    tsk_table_collection_free(pTables);
    delete pTables;
    pTables = nullptr;
  }
}

void TsRecorder::preRegisterSamples(NodePtr * pSampleNodes, unsigned int nSamples) {
  if (bSamplesPreRegistered) {
    return;
  }
  if (pTables == nullptr) {
    throw std::runtime_error("TS tables are not initialized");
  }
  if (!bInbred && (nSamples % iPloidy != 0)) {
    throw std::runtime_error("Sample count is not divisible by ploidy in outbred TS mode");
  }
  
  sampleNodeIds.clear();
  sampleNodeIndividualMap.clear();
  sampleNodeIds.reserve(nSamples);
  
  const unsigned int nIndividuals = bInbred ? nSamples : nSamples / iPloidy;
  std::vector<tsk_id_t> individualIds;
  individualIds.reserve(nIndividuals);
  for (unsigned int i = 0; i < nIndividuals; ++i) {
    tsk_id_t indivId = tsk_individual_table_add_row(&pTables->individuals,
                                                    0,
                                                    nullptr,
                                                    0,
                                                    nullptr,
                                                    0,
                                                    nullptr,
                                                    0);
    checkTsk(static_cast<int>(indivId), "Failed to add individual row");
    individualIds.push_back(indivId);
  }
  
  for (unsigned int i = 0; i < nSamples; ++i) {
    NodePtr & node = pSampleNodes[i];
    if (!node) {
      continue;
    }
    const unsigned long long nodeKey = node->getId();
    const unsigned int individualIndex = bInbred ? i : i / iPloidy;
    sampleNodeIndividualMap[nodeKey] = individualIds[individualIndex];
    // Preserve runMacs sample order in TS node ids by creating sample rows first.
    sampleNodeIds.push_back(getOrCreateNode(node));
  }
  bSamplesPreRegistered = true;
}

double TsRecorder::toTsPosition(double posUnit) const {
  double out = positionMode == TsPositionMode::PHYSICAL_BP ?
    posUnit * dSequenceLengthBp : posUnit;
  const double seqLength = positionMode == TsPositionMode::PHYSICAL_BP ?
    dSequenceLengthBp : 1.0;
  if (out < 0.0) {
    out = 0.0;
  } else if (out > seqLength) {
    out = seqLength;
  }
  return out;
}

tsk_id_t TsRecorder::ensurePopulation(short int population) {
  if (population < 0) {
    return TSK_NULL;
  }
  while (pTables->populations.num_rows <= static_cast<tsk_size_t>(population)) {
    tsk_id_t newPopulation = tsk_population_table_add_row(&pTables->populations,
                                                          nullptr, 0);
    checkTsk(static_cast<int>(newPopulation), "Failed to add population row");
  }
  return static_cast<tsk_id_t>(population);
}

tsk_id_t TsRecorder::getOrCreateNode(NodePtr & node) {
  const unsigned long long nodeKey = node->getId();
  const auto it = nodeIdMap.find(nodeKey);
  if (it != nodeIdMap.end()) {
    return it->second;
  }
  
  tsk_flags_t flags = 0;
  if (node->getType() == Node::SAMPLE) {
    flags |= TSK_NODE_IS_SAMPLE;
  }
  const tsk_id_t population = ensurePopulation(node->getPopulation());
  tsk_id_t individual = TSK_NULL;
  if (node->getType() == Node::SAMPLE) {
    const auto sampleIt = sampleNodeIndividualMap.find(nodeKey);
    if (sampleIt != sampleNodeIndividualMap.end()) {
      individual = sampleIt->second;
    }
  }
  tsk_id_t nodeId = tsk_node_table_add_row(&pTables->nodes,
                                           flags,
                                           node->getHeight(),
                                           population,
                                           individual,
                                           nullptr,
                                           0);
  checkTsk(static_cast<int>(nodeId), "Failed to add node row");
  nodeIdMap.insert(std::make_pair(nodeKey, nodeId));
  return nodeId;
}

void TsRecorder::recordTreeInterval(const EdgePtrVector & treeEdges,
                                    unsigned int iTotalTreeEdges,
                                    double leftPosUnit, double rightPosUnit) {
  const double left = toTsPosition(leftPosUnit);
  const double right = toTsPosition(rightPosUnit);
  if (right <= left) {
    return;
  }
  
  for (unsigned int i = 0; i < iTotalTreeEdges; ++i) {
    EdgePtr edge = treeEdges[i];
    if (edge->bDeleted) {
      continue;
    }
    NodePtr & parentNode = edge->getTopNodeRef();
    NodePtr & childNode = edge->getBottomNodeRef();
    const tsk_id_t parent = getOrCreateNode(parentNode);
    const tsk_id_t child = getOrCreateNode(childNode);
    tsk_id_t edgeId = tsk_edge_table_add_row(&pTables->edges,
                                             left,
                                             right,
                                             parent,
                                             child,
                                             nullptr,
                                             0);
    checkTsk(static_cast<int>(edgeId), "Failed to add edge row");
  }
}

void TsRecorder::recordMutation(double mutationPosUnit, EdgePtr & selectedEdge,
                                double mutationTime) {
  double position = toTsPosition(mutationPosUnit);
  const double seqLength = positionMode == TsPositionMode::PHYSICAL_BP ?
    dSequenceLengthBp : 1.0;
  if (position >= seqLength) {
    // Enforce [0, sequence_length) constraint used by tskit for site positions.
    position = std::nextafter(seqLength, 0.0);
  }
  
  static const char ancestralState[] = "0";
  static const char derivedState[] = "1";
  tsk_id_t siteId = tsk_site_table_add_row(&pTables->sites,
                                           position,
                                           ancestralState,
                                           1,
                                           nullptr,
                                           0);
  checkTsk(static_cast<int>(siteId), "Failed to add site row");
  
  NodePtr & childNode = selectedEdge->getBottomNodeRef();
  const tsk_id_t nodeId = getOrCreateNode(childNode);
  tsk_id_t mutationId = tsk_mutation_table_add_row(&pTables->mutations,
                                                   siteId,
                                                   nodeId,
                                                   TSK_NULL,
                                                   mutationTime,
                                                   derivedState,
                                                   1,
                                                   nullptr,
                                                   0);
  checkTsk(static_cast<int>(mutationId), "Failed to add mutation row");
}

void TsRecorder::simplify() {
  if (pTables == nullptr || bSimplified) {
    return;
  }
  if (sampleNodeIds.empty()) {
    bSimplified = true;
    return;
  }
  checkTsk(tsk_table_collection_sort(pTables, nullptr, 0),
           "Failed to sort table collection before simplify");
  checkTsk(tsk_table_collection_simplify(pTables,
                                         sampleNodeIds.data(),
                                         sampleNodeIds.size(),
                                         0,
                                         nullptr),
           "Failed to simplify table collection");
  bSimplified = true;
}

void TsRecorder::expandInbredSamples() {
  if (!bInbred || iPloidy <= 1) {
    return;
  }
  tsPost::expandInbredSamplesInPlace(pTables, iPloidy);
}

tsk_table_collection_t * TsRecorder::release(double timeScale, bool expandInbred) {
  if (pTables == nullptr) {
    return nullptr;
  }
  if (expandInbred) {
    expandInbredSamples();
  }
  rescaleTableCollectionTimes(pTables, timeScale);
  checkTsk(tsk_table_collection_sort(pTables, nullptr, 0),
           "Failed to sort table collection");
  checkTsk(tsk_table_collection_build_index(pTables, 0),
           "Failed to build table collection index");
  tsk_table_collection_t * out = pTables;
  pTables = nullptr;
  return out;
}
