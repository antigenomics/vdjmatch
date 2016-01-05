package com.antigenomics.vdjdb.fixer

class KmerScanner {
    final int minHitSize
    final Map<String, Integer> kmers = new HashMap<>();

    KmerScanner(String seq, int minHitSize = 2) {
        this.minHitSize = minHitSize

        for (int i = minHitSize; i < seq.length(); i++) {
            for (int j = 0; j < seq.length() - i; j++) {
                String kmer = seq.substring(j, j + i)
                kmers.put(kmer, j)
            }
        }
    }

    SearchResult scan(String seq) {
        // iterate from largest window to smallest one
        for (int i = seq.length(); i >= minHitSize; i--) {
            // sliding window scan
            for (int j = 0; j < seq.length() - i; j++) {
                String kmer = seq.substring(j, j + i)

                def hit = kmers[kmer]

                if (hit != null) {
                    return new SearchResult(hit,
                            j,
                            kmer.length()
                    )
                }
            }
        }

        null
    }
}
