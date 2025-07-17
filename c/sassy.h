#ifndef SASSY_H
#define SASSY_H

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#define sassy_Searcher_CHECK_AT_LEAST_ROWS 8

typedef struct sassy_SearcherType sassy_SearcherType;

typedef struct sassy_CMatch {
  int32_t pattern_start;
  int32_t text_start;
  int32_t pattern_end;
  int32_t text_end;
  int32_t cost;
  uint8_t strand;
} sassy_CMatch;

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * Create a new `Searcher` instance.
 *
 * * `alphabet` – one of "ascii", "dna", "iupac" (case-insensitive).
 * * `rc`       – whether to also search the reverse-complement (- strand).
 * * `alpha`    – overhang parameter. Pass `NAN` if you do not want to set it.
 *
 * Returns a pointer to an opaque `Searcher` or NULL on error.
 */
struct sassy_SearcherType *sassy_searcher(const char *alphabet, bool rc, float alpha);

/**
 * Free a `Searcher` previously created with `sassy_searcher`.
 */
void sassy_searcher_free(struct sassy_SearcherType *ptr);

/**
 * Search for `pattern` inside `text` allowing up to `k` edits.
 *
 * `out_matches` will point to a newly allocated array of `CMatch` results. The
 * function returns the number of matches found. The caller takes ownership of
 * the array and must later free it using `sassy_matches_free`.
 */
uintptr_t search(struct sassy_SearcherType *ptr,
                                const uint8_t *pattern,
                                uintptr_t pattern_len,
                                const uint8_t *text,
                                uintptr_t text_len,
                                uintptr_t k,
                                struct sassy_CMatch **out_matches);

/**
 * Free a match array previously obtained from `sassy_searcher_search`.
 */
void sassy_matches_free(struct sassy_CMatch *ptr, uintptr_t len);

#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus

#endif  /* SASSY_H */
