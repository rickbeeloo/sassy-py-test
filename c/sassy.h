#ifndef SASSY_H
#define SASSY_H

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct sassy_SearcherType sassy_SearcherType;

typedef struct sassy_Match {
  uintptr_t text_start;
  uintptr_t text_end;
  uintptr_t pattern_start;
  uintptr_t pattern_end;
  int32_t cost;
  /**
   * 0 = Fwd, 1 = Rc
   */
  uint8_t strand;
} sassy_Match;



#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * Create a new `Searcher` instance.
 *
 * `alphabet`: one of "ascii", "dna", "iupac" (case-insensitive).
 * `rc`: whether to also search the reverse-complement strand.
 * `alpha`: overhang parameter. Pass `NAN` to disable.
 *
 * Returns a pointer to an opaque `Searcher` object, or panics on error.
 */
struct sassy_SearcherType *sassy_searcher(const char *alphabet, bool rc, float alpha);

/**
 * Free a `Searcher` previously created with `sassy_searcher`.
 */
void sassy_searcher_free(struct sassy_SearcherType *ptr);

/**
 * Search for `pattern` in `text` allowing up to `k` edits.
 *
 * `out_matches` will point to a newly allocated rust `Vec` of `Match` results. The
 * function returns the number of matches found.
 * Matches should be freed using `sassy_matches_free`.
 */
uintptr_t search(struct sassy_SearcherType *searcher,
                 const uint8_t *pattern,
                 uintptr_t pattern_len,
                 const uint8_t *text,
                 uintptr_t text_len,
                 uintptr_t k,
                 struct sassy_Match **out_matches);

/**
 * Free a match array previously obtained from `sassy_searcher_search`.
 */
void sassy_matches_free(struct sassy_Match *ptr, uintptr_t len);

#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus

#endif  /* SASSY_H */
