package com.antigenomics.vdjdb.core2.text

import com.antigenomics.vdjdb.core2.db.Column
import com.antigenomics.vdjdb.core2.db.ColumnType
import com.antigenomics.vdjdb.core2.db.Entry

class TextColumn extends Column {
    private final Map<String, List<Entry>> map = new HashMap<>()

    TextColumn(String name, List<String> metadata) {
        super(name, metadata, ColumnType.Text)
    }

    @Override
    void append(Entry entry) {
        if (entry.value.length() > 0) {
            def entries = map[entry.value]
            if (entries == null) {
                map.put(entry.value, entries = new ArrayList<Entry>())
            }
            entries.add(entry)
        }
    }
}
