#ifndef ALPHASIMR_POST_TS_H
#define ALPHASIMR_POST_TS_H

#include <cstdint>

#include "tskit.h"

namespace tsPost {

void checkTsk(int status, const char *context);

void expandInbredSamplesInPlace(tsk_table_collection_t *tables,
                                unsigned int ploidy);

void mutateTablesInPlace(tsk_table_collection_t *tables,
                         double theta,
                         uint64_t seed);

} // namespace tsPost

#endif
