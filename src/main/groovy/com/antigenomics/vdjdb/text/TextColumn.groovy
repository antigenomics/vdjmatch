package com.antigenomics.vdjdb.text

import com.antigenomics.vdjdb.db.Column
import com.antigenomics.vdjdb.db.ColumnType
import com.antigenomics.vdjdb.db.Entry

class TextColumn extends Column {
    private final Map<String, List<Entry>> map = new HashMap<>()

    TextColumn(String name, Map<String, String> metadata = [:]) {
        super(name, metadata, ColumnType.Text)
    }

    Set<String> getValues() {
        map.keySet()
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
