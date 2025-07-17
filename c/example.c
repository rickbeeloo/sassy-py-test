
#include "sassy.h"
#include <math.h>     
#include <stdio.h>
#include <string.h>   

static void
print_match(const sassy_CMatch *m, size_t idx)
{
    printf("#%zu  pat[%d-%d]  txt[%d-%d]  cost=%d  strand=%c\n",
           idx,
           m->pattern_start,
           m->pattern_end,
           m->text_start,
           m->text_end,
           m->cost,
           m->strand == 0 ? '+' : '-');
}

int
main(void)
{
    /* ---------------------------------------------------------------------
     * 1.  Create a searcher (DNA alphabet, search reverse-complement too)
     * ------------------------------------------------------------------ */
    sassy_SearcherType *s =
        sassy_searcher(/* alphabet */ "dna",
                           /* rc       */ true,
                           /* alpha    */ NAN);      /* no overhang */
    if (!s) {
        fprintf(stderr, "Couldn’t create searcher – bad alphabet?\n");
        return 1;
    }

    /* ---------------------------------------------------------------------
     * 2.  Define pattern / text and run the search
     * ------------------------------------------------------------------ */
    const char  *pattern = "AAGGGGA";
    const char  *text    = "CCCCCCCCCAAGGGGACCCCCAAGGCGACCCCCCCCC";
    const size_t k       = 1;      /* maximum edit distance */

    sassy_CMatch *matches = NULL;
    size_t n =
        search(s,
                              (const uint8_t *)pattern, strlen(pattern),
                              (const uint8_t *)text,    strlen(text),
                              k,
                              &matches);

    printf("Found %zu match(es):\n", n);
    for (size_t i = 0; i < n; i++)
        print_match(&matches[i], i);

    /* ---------------------------------------------------------------------
     * 3.  Clean up
     * ------------------------------------------------------------------ */
    sassy_matches_free(matches, n);
    sassy_searcher_free(s);

    return 0;
}