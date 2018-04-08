package com.antigenomics.vdjdb;

import java.util.HashSet;
import java.util.Set;
import java.util.function.Function;

public class HashSetStringGen implements Function<String, Set<String>> {
    public static final HashSetStringGen INSTANCE = new HashSetStringGen();
    private HashSetStringGen() {

    }

    @Override
    public Set<String> apply(String s) {
        return new HashSet<>();
    }
}
