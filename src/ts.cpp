#include "alphasimr.h"
#include <RcppTskit.hpp>
#include "postTS.h"
#include <cmath>
#include <limits>
#include <string>
#include <vector>
// [[Rcpp::depends(RcppTskit)]]
// [[Rcpp::plugins(RcppTskit)]]

namespace {

inline SEXP asTableCollectionXptr(const SEXP tc) {
  if (TYPEOF(tc) == EXTPTRSXP) {
    return tc;
  }
  if (TYPEOF(tc) == ENVSXP) {
    Rcpp::Environment env(tc);
    if (env.exists("xptr")) {
      return env["xptr"];
    }
  }
  if (TYPEOF(tc) == VECSXP) {
    Rcpp::List lst(tc);
    if (lst.containsElementNamed("xptr")) {
      return lst["xptr"];
    }
  }
  Rcpp::stop("tc must be a table-collection external pointer or an object with `$xptr`");
}

inline void stopIfTskError(const tsk_id_t id, const char * context) {
  if (id < 0) {
    Rcpp::stop("%s: %s", context, tsk_strerror(static_cast<int>(id)));
  }
}

} // namespace

// TODO: This is just an example - we will replace it later with other
//       functions working with tree sequences. For example to obtain haplotypes
//       using https://tskit.dev/tskit/docs/stable/c-api.html#decoding-genotypes
//
//       See also ts: Create haplotypes from a tree sequence for downstream
//       AlphaSimR work (but subsample sites!) #4
//       https://github.com/HighlanderLab/AlphaSimR/issues/4
//
//' @title Summarise `tskit` table collection
//' @param tc an external pointer to a \code{tsk_table_collection_t} object.
//' @return A list.
//' @examples
//' ts_file <- system.file("examples", "test.trees", package = "RcppTskit")
//' tc <- RcppTskit:::tc_load(ts_file)
//' RcppTskit:::rtsk_table_collection_summary(tc$xptr)
//' rtsk_table_collection_summary2(tc$xptr)
//' @export
// [[Rcpp::export]]
Rcpp::List rtsk_table_collection_summary2(const SEXP tc) {
  rtsk_table_collection_t tc_xptr(tc);
  const tsk_table_collection_t *tables = tc_xptr;
  return Rcpp::List::create(
      Rcpp::_["num_provenances"] = tables->provenances.num_rows,
      Rcpp::_["num_populations"] = tables->populations.num_rows,
      Rcpp::_["num_migrations"] = tables->migrations.num_rows,
      Rcpp::_["num_individuals"] = tables->individuals.num_rows,
      Rcpp::_["num_nodes"] = tables->nodes.num_rows,
      Rcpp::_["num_edges"] = tables->edges.num_rows,
      Rcpp::_["num_sites"] = tables->sites.num_rows,
      Rcpp::_["num_mutations"] = tables->mutations.num_rows,
      Rcpp::_["sequence_length"] = tables->sequence_length,
      Rcpp::_["has_reference_sequence"] =
          rtsk_table_collection_has_reference_sequence(tc),
      Rcpp::_["time_units"] = rtsk_table_collection_get_time_units(tc),
      Rcpp::_["file_uuid"] = rtsk_table_collection_get_file_uuid(tc),
      Rcpp::_["has_index"] = rtsk_table_collection_has_index(tc));
}

//' @title Get number of individuals in tree sequence
//' @param ts an external pointer to a \code{tsk_treeseq_t} object.
//' @return integer number of individuals.
//' @examples
//' ts_file <- system.file("examples", "test.trees", package = "RcppTskit")
//' ts <- RcppTskit::ts_load(ts_file)
//' ts$num_individuals()
//' RcppTskit:::rtsk_treeseq_get_num_individuals(ts$xptr)
//' rtsk_treeseq_get_num_individuals2(ts$xptr)
//' @export
// [[Rcpp::export]]
int rtsk_treeseq_get_num_individuals2(const SEXP ts) {
  rtsk_treeseq_t ts_xptr(ts);
  return static_cast<int>(tsk_treeseq_get_num_individuals(ts_xptr));
}

// [[Rcpp::export]]
void tsMutateTableCollection(const SEXP tc, const double theta,
                             const uint64_t seed) {
  rtsk_table_collection_t tc_xptr(asTableCollectionXptr(tc));
  tsk_table_collection_t *tables = tc_xptr;
  tsPost::mutateTablesInPlace(tables, theta, seed);
}

// [[Rcpp::export]]
void tsFinalizeInbredTableCollection(const SEXP tc, const int ploidy) {
  if (ploidy <= 1) {
    return;
  }
  rtsk_table_collection_t tc_xptr(asTableCollectionXptr(tc));
  tsk_table_collection_t *tables = tc_xptr;
  if (tables == nullptr) {
    Rcpp::stop("Table collection pointer is null");
  }
  
  tsPost::expandInbredSamplesInPlace(tables, static_cast<unsigned int>(ploidy));
  (void)tsk_table_collection_drop_index(tables, 0);
  tsPost::checkTsk(tsk_table_collection_sort(tables, nullptr, 0),
                   "Failed to sort table collection after inbred finalization");
  tsPost::checkTsk(tsk_table_collection_build_index(tables, 0),
                   "Failed to build index after inbred finalization");
}

// [[Rcpp::export]]
SEXP vcfFounderTableCollection(const Rcpp::IntegerMatrix haplo,
                               const Rcpp::NumericVector pos,
                               const double seqLen,
                               const int ploidy,
                               const bool addTsMut) {
  const int nHaplo = haplo.nrow();
  const int nSites = haplo.ncol();
  if (ploidy <= 0) {
    Rcpp::stop("ploidy must be a positive integer");
  }
  if (nHaplo <= 0 || nSites <= 0) {
    Rcpp::stop("haplo must have at least one row and one column");
  }
  if (nHaplo % ploidy != 0) {
    Rcpp::stop("number of haplotypes must be divisible by ploidy");
  }
  if (pos.size() != nSites) {
    Rcpp::stop("pos length must match ncol(haplo)");
  }
  if (!std::isfinite(seqLen) || seqLen <= 0.0) {
    Rcpp::stop("seqLen must be a finite positive value");
  }

  tsk_table_collection_t *tables = new tsk_table_collection_t;
  int ret = tsk_table_collection_init(tables, 0);
  if (ret < 0) {
    delete tables;
    Rcpp::stop("Failed to initialise tsk tables: %s", tsk_strerror(ret));
  }
  tables->sequence_length = seqLen;

  try {
    static const char generationUnits[] = "generations";
    ret = tsk_table_collection_set_time_units(
      tables, generationUnits, sizeof(generationUnits) - 1);
    if (ret < 0) {
      Rcpp::stop("Failed to set TS time_units: %s", tsk_strerror(ret));
    }

    const int nInd = nHaplo / ploidy;
    std::vector<tsk_id_t> individualIds;
    individualIds.reserve(static_cast<std::size_t>(nInd));
    for (int ind = 0; ind < nInd; ++ind) {
      const tsk_id_t individualId = tsk_individual_table_add_row(
        &tables->individuals,
        0,
        nullptr, 0,
        nullptr, 0,
        nullptr, 0);
      stopIfTskError(individualId, "Failed to add individual row");
      individualIds.push_back(individualId);
    }

    std::vector<tsk_id_t> sampleNodeIds;
    sampleNodeIds.reserve(static_cast<std::size_t>(nHaplo));
    for (int h = 0; h < nHaplo; ++h) {
      const int ind = h / ploidy;
      const tsk_id_t nodeId = tsk_node_table_add_row(
        &tables->nodes,
        TSK_NODE_IS_SAMPLE,
        0.0,
        TSK_NULL,
        individualIds[static_cast<std::size_t>(ind)],
        nullptr, 0);
      stopIfTskError(nodeId, "Failed to add sample node row");
      sampleNodeIds.push_back(nodeId);
    }

    if (addTsMut) {
      double lastPos = -std::numeric_limits<double>::infinity();
      for (int site = 0; site < nSites; ++site) {
        double sitePos = pos[site];
        if (!std::isfinite(sitePos)) {
          Rcpp::stop("site positions must be finite");
        }
        if (sitePos < 0.0 || sitePos > seqLen) {
          Rcpp::stop("site position outside [0, seqLen]");
        }
        if (sitePos >= seqLen) {
          sitePos = std::nextafter(seqLen, 0.0);
        }
        if (sitePos <= lastPos) {
          Rcpp::stop("site positions must be strictly increasing");
        }
        lastPos = sitePos;

        static const char ancestralState[] = "0";
        const tsk_id_t siteId = tsk_site_table_add_row(
          &tables->sites,
          sitePos,
          ancestralState,
          1,
          nullptr, 0);
        stopIfTskError(siteId, "Failed to add site row");

        static const char derivedState[] = "1";
        for (int h = 0; h < nHaplo; ++h) {
          const int allele = haplo(h, site);
          if (allele == NA_INTEGER || (allele != 0 && allele != 1)) {
            Rcpp::stop("haplo entries must be 0 or 1");
          }
          if (allele == 1) {
            const tsk_id_t mutationId = tsk_mutation_table_add_row(
              &tables->mutations,
              siteId,
              sampleNodeIds[static_cast<std::size_t>(h)],
              TSK_NULL,
              TSK_UNKNOWN_TIME,
              derivedState,
              1,
              nullptr, 0);
            stopIfTskError(mutationId, "Failed to add mutation row");
          }
        }
      }
    }

    ret = tsk_table_collection_sort(tables, nullptr, 0);
    if (ret < 0) {
      Rcpp::stop("Failed to sort table collection: %s", tsk_strerror(ret));
    }
    ret = tsk_table_collection_build_index(tables, 0);
    if (ret < 0) {
      Rcpp::stop("Failed to build table collection index: %s", tsk_strerror(ret));
    }
  } catch (...) {
    tsk_table_collection_free(tables);
    delete tables;
    throw;
  }

  rtsk_table_collection_t out(tables, true);
  return out;
}

// [[Rcpp::export]]
Rcpp::IntegerVector tsForwardNodeTableAddRows(const SEXP tc,
                                              const Rcpp::IntegerVector flags,
                                              const Rcpp::NumericVector time,
                                              const Rcpp::IntegerVector population,
                                              const Rcpp::IntegerVector individual) {
  const R_xlen_t n = time.size();
  if (flags.size() != n || population.size() != n || individual.size() != n) {
    Rcpp::stop("flags/time/population/individual must have identical lengths");
  }
  rtsk_table_collection_t tc_xptr(asTableCollectionXptr(tc));
  tsk_table_collection_t *tables = tc_xptr;
  if (tables == nullptr) {
    Rcpp::stop("Table collection pointer is null");
  }

  Rcpp::IntegerVector out(n);
  if (n == 0) {
    return out;
  }
  std::vector<tsk_flags_t> flagsData;
  std::vector<tsk_id_t> populationData;
  std::vector<tsk_id_t> individualData;
  flagsData.reserve(static_cast<std::size_t>(n));
  populationData.reserve(static_cast<std::size_t>(n));
  individualData.reserve(static_cast<std::size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) {
    const double t = time[i];
    if (!R_finite(t)) {
      Rcpp::stop("node time must be finite at row %d", static_cast<int>(i + 1));
    }
    flagsData.push_back(static_cast<tsk_flags_t>(flags[i]));
    populationData.push_back(static_cast<tsk_id_t>(population[i]));
    individualData.push_back(static_cast<tsk_id_t>(individual[i]));
  }

  const tsk_size_t start = tables->nodes.num_rows;
  const int ret = tsk_node_table_append_columns(
    &tables->nodes,
    static_cast<tsk_size_t>(n),
    flagsData.data(),
    REAL(time),
    populationData.data(),
    individualData.data(),
    nullptr,
    nullptr
  );
  stopIfTskError(static_cast<tsk_id_t>(ret), "Failed to append node rows");

  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = static_cast<int>(start + static_cast<tsk_size_t>(i));
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::IntegerVector tsForwardNodeTableAddRowsWithMetadata(
    const SEXP tc,
    const Rcpp::IntegerVector flags,
    const Rcpp::NumericVector time,
    const Rcpp::IntegerVector population,
    const Rcpp::IntegerVector individual,
    const Rcpp::CharacterVector nodeKey) {
  const R_xlen_t n = time.size();
  if (flags.size() != n || population.size() != n || individual.size() != n ||
      nodeKey.size() != n) {
    Rcpp::stop("flags/time/population/individual/nodeKey must have identical lengths");
  }
  rtsk_table_collection_t tc_xptr(asTableCollectionXptr(tc));
  tsk_table_collection_t *tables = tc_xptr;
  if (tables == nullptr) {
    Rcpp::stop("Table collection pointer is null");
  }

  Rcpp::IntegerVector out(n);
  if (n == 0) {
    return out;
  }

  std::vector<tsk_flags_t> flagsData;
  std::vector<tsk_id_t> populationData;
  std::vector<tsk_id_t> individualData;
  std::vector<char> metadataData;
  std::vector<tsk_size_t> metadataOffset;
  flagsData.reserve(static_cast<std::size_t>(n));
  populationData.reserve(static_cast<std::size_t>(n));
  individualData.reserve(static_cast<std::size_t>(n));
  metadataOffset.reserve(static_cast<std::size_t>(n) + 1);
  metadataOffset.push_back(0);

  for (R_xlen_t i = 0; i < n; ++i) {
    const double t = time[i];
    if (!R_finite(t)) {
      Rcpp::stop("node time must be finite at row %d", static_cast<int>(i + 1));
    }
    if (Rcpp::CharacterVector::is_na(nodeKey[i])) {
      Rcpp::stop("nodeKey must not be NA at row %d", static_cast<int>(i + 1));
    }

    flagsData.push_back(static_cast<tsk_flags_t>(flags[i]));
    populationData.push_back(static_cast<tsk_id_t>(population[i]));
    individualData.push_back(static_cast<tsk_id_t>(individual[i]));

    const std::string key = Rcpp::as<std::string>(nodeKey[i]);
    const std::string metadata = "{\"alphaSimR\":{\"id\":\"" + key + "\"}}";
    metadataData.insert(metadataData.end(), metadata.begin(), metadata.end());
    metadataOffset.push_back(static_cast<tsk_size_t>(metadataData.size()));
  }

  const tsk_size_t start = tables->nodes.num_rows;
  const int ret = tsk_node_table_append_columns(
    &tables->nodes,
    static_cast<tsk_size_t>(n),
    flagsData.data(),
    REAL(time),
    populationData.data(),
    individualData.data(),
    metadataData.data(),
    metadataOffset.data()
  );
  stopIfTskError(static_cast<tsk_id_t>(ret), "Failed to append node rows with metadata");

  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = static_cast<int>(start + static_cast<tsk_size_t>(i));
  }
  return out;
}

// [[Rcpp::export]]
void tsForwardEdgeTableAddRows(const SEXP tc,
                               const Rcpp::NumericVector left,
                               const Rcpp::NumericVector right,
                               const Rcpp::IntegerVector parent,
                               const Rcpp::IntegerVector child) {
  const R_xlen_t n = left.size();
  if (right.size() != n || parent.size() != n || child.size() != n) {
    Rcpp::stop("left/right/parent/child must have identical lengths");
  }
  rtsk_table_collection_t tc_xptr(asTableCollectionXptr(tc));
  tsk_table_collection_t *tables = tc_xptr;
  if (tables == nullptr) {
    Rcpp::stop("Table collection pointer is null");
  }

  std::vector<double> leftData;
  std::vector<double> rightData;
  std::vector<tsk_id_t> parentData;
  std::vector<tsk_id_t> childData;
  leftData.reserve(static_cast<std::size_t>(n));
  rightData.reserve(static_cast<std::size_t>(n));
  parentData.reserve(static_cast<std::size_t>(n));
  childData.reserve(static_cast<std::size_t>(n));

  for (R_xlen_t i = 0; i < n; ++i) {
    const double l = left[i];
    const double r = right[i];
    if (!R_finite(l) || !R_finite(r) || r <= l) {
      continue;
    }
    leftData.push_back(l);
    rightData.push_back(r);
    parentData.push_back(static_cast<tsk_id_t>(parent[i]));
    childData.push_back(static_cast<tsk_id_t>(child[i]));
  }

  if (leftData.empty()) {
    return;
  }

  const int ret = tsk_edge_table_append_columns(
    &tables->edges,
    static_cast<tsk_size_t>(leftData.size()),
    leftData.data(),
    rightData.data(),
    parentData.data(),
    childData.data(),
    nullptr,
    nullptr
  );
  stopIfTskError(static_cast<tsk_id_t>(ret), "Failed to append edge rows");
}

// [[Rcpp::export]]
void tsForwardSetSampleFlags(const SEXP tc,
                             const Rcpp::IntegerVector samples,
                             const bool clearExisting = true) {
  rtsk_table_collection_t tc_xptr(asTableCollectionXptr(tc));
  tsk_table_collection_t *tables = tc_xptr;
  if (tables == nullptr) {
    Rcpp::stop("Table collection pointer is null");
  }

  const tsk_size_t nNodes = tables->nodes.num_rows;
  if (clearExisting) {
    for (tsk_size_t i = 0; i < nNodes; ++i) {
      tables->nodes.flags[i] &= ~static_cast<tsk_flags_t>(TSK_NODE_IS_SAMPLE);
    }
  }

  for (R_xlen_t i = 0; i < samples.size(); ++i) {
    const int node = samples[i];
    if (node == NA_INTEGER) {
      continue;
    }
    if (node < 0 || static_cast<tsk_size_t>(node) >= nNodes) {
      Rcpp::stop("sample node id out of bounds at index %d", static_cast<int>(i + 1));
    }
    tables->nodes.flags[node] |= static_cast<tsk_flags_t>(TSK_NODE_IS_SAMPLE);
  }
}
