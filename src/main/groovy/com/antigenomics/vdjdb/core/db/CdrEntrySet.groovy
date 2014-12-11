package com.antigenomics.vdjdb.core.db

import groovy.transform.PackageScope


class CdrEntrySet implements Iterable<CdrEntry> {
    private final CdrDatabase parent
    private final List<CdrEntry> entries = new LinkedList<>()
    private final String cdr3aa

    public CdrEntrySet(CdrDatabase parent, String cdr3aa) {
        this.parent = parent
        this.cdr3aa = cdr3aa
    }

    @PackageScope
    void addEntry(String v, String j, String[] annotation) {
        if (parent.annotationHeader.length != annotation.length)
            throw new IndexOutOfBoundsException("Annotation header length don't match provided annotation row")
        entries.add(new CdrEntry(v, j, annotation, this))
    }

    String getCdr3aa() {
        return cdr3aa
    }

    public CdrDatabase getParent() {
        return parent
    }

    @Override
    public Iterator iterator() {
        entries.iterator()
    }

    @Override
    public boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        CdrEntrySet entries = (CdrEntrySet) o

        cdr3aa == entries.cdr3aa && parent == entries.parent
    }

    @Override
    public int hashCode() {
        31 * parent.hashCode() + cdr3aa.hashCode()
    }
}
