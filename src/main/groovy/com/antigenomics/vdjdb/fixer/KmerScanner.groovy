package com.antigenomics.vdjdb.fixer

class KmerScanner {
    final int minHitSize, refLength
    final boolean reverse
    final Map<String, Integer> kmers = new HashMap<>();

    KmerScanner(String seq, boolean reverse, int minHitSize = 2) {
        this.reverse = reverse
        this.minHitSize = minHitSize
        this.refLength = seq.length()

        seq = reverse ? seq.reverse() : seq

        for (int i = minHitSize; i < seq.length(); i++) {
            for (int j = 0; j < seq.length() - i; j++) {
                String kmer = seq.substring(j, j + i)
                kmers.put(kmer, j)
            }
        }
    }

    SearchResult scan(String seq) {
        seq = reverse ? seq.reverse() : seq

        // iterate from largest window to smallest one
        for (int i = seq.length(); i >= minHitSize; i--) {
            // sliding window scan
            for (int j = 0; j < seq.length() - i; j++) {
                String kmer = seq.substring(j, j + i)

                def hit = kmers[kmer]

                if (hit) {
                    return new SearchResult(reverse ? refLength - hit - 1 : hit,
                            reverse ? seq.length() - i - 1 : i,
                            kmer.length()
                    )
                }
            }
        }

        null
    }
}
