/**
   Find longest common substring of character(n) using suffix array
*/

#include <string.h>             // strcmp

#include <vector>
#include <algorithm>
#include <numeric>

inline int lcp(const char *a, const char *b) {
    int i = 0;
    while (a[i] == b[i] && a[i] != '\0' && b[i] != '\0')
        i += 1;
    return i;
}

inline bool cmp_less_than(const char *a, const char *b) {
    return strcmp(a, b) < 0;
}

inline bool cmp_equal(const char *a, const char *b) {
    return strcmp(a, b) == 0;
}

#include <Rinternals.h>

SEXP longest_common_substring(SEXP x_sexp)
{
    if (!Rf_isString(x_sexp))
        Rf_error("'x' must be a character vector");

    /* Create suffix array */
    R_len_t n = Rf_length(x_sexp);
    size_t tot_len = 0;
    std::vector<int> len(n, 0);
    for (R_len_t i = 0; i < n; ++i) {
        tot_len += Rf_length(STRING_ELT(x_sexp, i)) + 1; // null-terminate
        len[i] = tot_len;
    }

    std::vector<char> cstr(tot_len, '\0');
    for (R_len_t i = 0; i < n; ++i) {
        size_t off = i == 0 ? 0 : len[i - 1], n_elt = len[i] - off - 1;
        memcpy(&cstr[off], CHAR(STRING_ELT(x_sexp, i)), n_elt * sizeof(char));
    }

    std::vector<const char *> array(tot_len - n, NULL);
    std::vector<int> interval(tot_len, 0);
    int index = 0;
    for (size_t i = 0; i < tot_len; ++i) {
        if (i == len[index] - 1) { // ignore null termination character
            index += 1;
            continue;
        }
        array[i - index] = &cstr[i];
        interval[i] = index;
    }
    std::sort(array.begin(), array.end(), cmp_less_than);

    /* Find longest prefix */
    std::vector<const char *> found;
    int lcs_len = 0;
    for (size_t i = 0; i + 1 < array.size(); ++i) {
        int i0 = interval[array[i] - &cstr[0]],
            i1 = interval[array[i + 1] - &cstr[0]];
        if (i0 == i1)           // same character element
            continue;
        int lcp_len = lcp(array[i], array[i + 1]);
        if (lcp_len > lcs_len) { // longer prefix found
            found.clear();
            lcs_len = lcp_len;
        }
        if (lcp_len && (lcp_len == lcs_len))
            found.push_back(array[i]);
    }

    /* Post-processing / return */
    std::sort(found.begin(), found.end(), cmp_less_than);
    found.erase(unique(found.begin(), found.end(), cmp_equal), found.end());

    SEXP result = PROTECT(Rf_allocVector(STRSXP, found.size()));
    for (size_t i = 0; i < found.size(); ++i)
        SET_STRING_ELT(result, i, mkCharLen(found[i], lcs_len));

    UNPROTECT(1);
    return result;
}

#include <R_ext/Rdynload.h>

extern "C" {

    static const R_CallMethodDef callMethods[] = {
        {".longest_common_substring", (DL_FUNC) & longest_common_substring, 1},
        {NULL, NULL, 0}
    };

    void R_init_microRNA(DllInfo * info)
    {
        R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    }

}
