package com.antigenomics.vdjdb.core.db


class CdrEntry {
    public final String v, j, key
    private final String[] annotation
    public final CdrEntrySet parent

    public CdrEntry(String v, String j, String[] annotation, CdrEntrySet parent) {
        this.v = v
        this.j = j
        this.annotation = annotation
        this.parent = parent
        this.key = v + "\t" + j + "\t" + annotation.collect().join("\t")
    }

    public String getAt(String field) {
        int index = parent.parent.getIndex(field)
        if (index < 0)
            throw new Exception("Field $field not found within annotations. " +
                    "Allowed fileds: ${parent.parent.annotationHeader.join(", ")}")
        getAt(index)
    }

    public String getAt(int index) {
        annotation[index]
    }

    public List<String> getAnnotation() {
        annotation.collect()
    }

    public String getCdr3aa() {
        parent.cdr3aa
    }

    public CdrEntrySet getParent() {
        parent
    }

    @Override
    public boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        CdrEntry cdrEntry = (CdrEntry) o

        key == cdrEntry.key && parent == cdrEntry.parent
    }

    @Override
    public int hashCode() {
        31 * key.hashCode() + parent.hashCode()
    }

    @Override
    public String toString() {
        key
    }
}
