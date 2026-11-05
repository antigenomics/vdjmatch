package com.antigenomics.vdjdb

import org.junit.Test
import java.nio.file.Files
import java.util.zip.GZIPOutputStream
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream

class ReleaseTransportTest {
    @Test
    void historicalAndVersionedGzipDownloadsAreAtomicAndReadable() {
        def root = Files.createTempDirectory('vdjmatch-transport-').toFile()
        try {
            def pointer = new File(root, 'pointer.txt')
            def home = new File(root, 'installed')
            ['flat', 'dated'].each { layout ->
                def archive = new File(root, layout + '.zip')
                archive.withOutputStream { out ->
                    new ZipOutputStream(out).withCloseable { zip ->
                        ['vdjdb.txt', 'vdjdb.slim.txt', 'vdjdb.meta.txt', 'vdjdb.slim.meta.txt'].each { name ->
                            zip.putNextEntry(new ZipEntry((layout == 'dated' ? 'vdjdb-legacy-2026-11-06/' : '') + name + (layout == 'dated' ? '.gz' : '')))
                            if (layout == 'dated') {
                                def bytes = new ByteArrayOutputStream()
                                new GZIPOutputStream(bytes).withCloseable { it.write((name + '\n').getBytes('UTF-8')) }
                                zip.write(bytes.toByteArray())
                            } else zip.write((name + '\n').getBytes('UTF-8'))
                            zip.closeEntry()
                        }
                    }
                }
                pointer.text = archive.toURI().toString() + '\n'
                assert Util.checkDatabase(home, true, pointer.toURI().toURL())
                assert !Util.checkDatabase(home, true, pointer.toURI().toURL())
                def installed = Util.databaseDirectory(home)
                assert new File(installed, 'vdjdb.slim.txt' + (layout == 'dated' ? '.gz' : '')).delete()
                assert Util.checkDatabase(home, true, pointer.toURI().toURL())
                Util.tableAsStream(new File(Util.databaseDirectory(home), 'vdjdb.txt')).withCloseable {
                    assert it.getText('UTF-8') == 'vdjdb.txt\n'
                }
            }
            def previous = new File(home, 'latest-version.txt').text
            ['../vdjdb.txt', 'prefix/../../vdjdb.txt', '/vdjdb.txt', 'prefix\\vdjdb.txt', 'incomplete.txt'].each { name ->
                def broken = new File(root, 'broken.zip')
                broken.withOutputStream { out ->
                    new ZipOutputStream(out).withCloseable { zip ->
                        zip.putNextEntry(new ZipEntry(name)); zip.write('bad'.bytes); zip.closeEntry()
                    }
                }
                pointer.text = broken.toURI().toString() + '\n'
                try {
                    Util.checkDatabase(home, true, pointer.toURI().toURL())
                    assert false: 'unsafe/incomplete download accepted'
                } catch (IOException expected) { }
                assert new File(home, 'latest-version.txt').text == previous
                Util.tableAsStream(new File(Util.databaseDirectory(home), 'vdjdb.slim.txt')).withCloseable {
                    assert it.getText('UTF-8') == 'vdjdb.slim.txt\n'
                }
            }
        } finally { root.deleteDir() }
    }
}
