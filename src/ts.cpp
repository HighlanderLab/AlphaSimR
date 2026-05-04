#include "alphasimr.h"
#include <RcppTskit.hpp>
// [[Rcpp::depends(RcppTskit)]]
// [[Rcpp::plugins(RcppTskit)]]

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
