package com.antigenomics.vdjdb.misc;

public interface DistanceEdge {
    int getFrom();
    int getTo();
    double getDissimilarity();
}
