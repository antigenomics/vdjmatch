/*
 * Copyright 2015 Mikhail Shugay
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


package com.antigenomics.vdjdb

import com.milaboratory.core.sequence.AminoAcidSequence
import java.nio.file.Files
import java.nio.file.StandardCopyOption
import java.security.MessageDigest
import java.util.zip.GZIPInputStream
import java.util.zip.ZipInputStream

class Util {
    /**
     * Parent directory of executable JAR
     */
    static
    final String HOME_DIR = new File(Util.class.protectionDomain.codeSource.location.path).parent.replaceAll("%20", " ")

    /**
     * Check if a local copy of the database exists and downloads one if required.
     * @param updateIfNewer perform update if the local copy is older than the one available online
     * @return true if a new database was downloaded
     */
    static boolean checkDatabase(boolean updateIfNewer = false) {
        checkDatabase(new File(HOME_DIR), updateIfNewer,
                new URL("https://raw.githubusercontent.com/antigenomics/vdjdb-db/master/latest-version.txt"))
    }

    static boolean checkDatabase(File home, boolean updateIfNewer, URL pointer) {
        def versionFile = new File(home, "latest-version.txt")
        def version = versionFile.isFile() ? versionFile.readLines("UTF-8").find { it.trim() } : null
        if (version && !updateIfNewer) return false
        def latest = pointer.withInputStream { it.readLines("UTF-8").find { line -> line.trim() } }
        if (!latest) throw new IOException("Empty database release pointer")
        if (version == latest && DATABASE_FILES.every { tableFile(new File(databaseDirectory(home), it)).isFile() })
            return false

        home.mkdirs()
        def staging = Files.createTempDirectory(home.toPath(), ".vdjdb-install-").toFile()
        try {
            def names = [] as Set
            def parents = [] as Set
            new URL(latest).withInputStream { raw ->
                new ZipInputStream(new BufferedInputStream(raw)).withCloseable { zip ->
                    def entry
                    while ((entry = zip.nextEntry) != null) {
                        def parts = entry.name.split("/", -1)
                        if (entry.name.startsWith("/") || entry.name.contains("\\") ||
                                parts.any { it == ".." || it == "." } || !names.add(entry.name))
                            throw new IOException("Unsafe or duplicate release member: " + entry.name)
                        def name = parts[-1]
                        def logical = name.endsWith(".gz") ? name[0..<-3] : name
                        if (!entry.directory && RELEASE_FILES.contains(logical)) {
                            if (tableFile(new File(staging, logical)).exists())
                                throw new IOException("Ambiguous release table: " + logical)
                            parents.add(parts.size() > 1 ? parts[0..<-1].join("/") : "")
                            new File(staging, name).withOutputStream { out -> out << zip }
                        }
                        zip.closeEntry() // validates the ZIP member CRC
                    }
                }
            }
            if (parents.size() != 1 || !DATABASE_FILES.every { tableFile(new File(staging, it)).isFile() })
                throw new IOException("Release does not contain one complete legacy database")
            def target = new File(home, releaseDirectory(latest))
            if (target.exists() && !DATABASE_FILES.every { tableFile(new File(target, it)).isFile() }) {
                def backup = new File(home, ".vdjdb-repair-" + UUID.randomUUID())
                Files.move(target.toPath(), backup.toPath(), StandardCopyOption.ATOMIC_MOVE)
                try {
                    Files.move(staging.toPath(), target.toPath(), StandardCopyOption.ATOMIC_MOVE)
                } catch (Exception error) {
                    Files.move(backup.toPath(), target.toPath(), StandardCopyOption.ATOMIC_MOVE)
                    throw error
                }
                backup.deleteDir()
            } else if (!target.exists()) {
                Files.move(staging.toPath(), target.toPath(), StandardCopyOption.ATOMIC_MOVE)
            }
            def marker = Files.createTempFile(home.toPath(), ".vdjdb-version-", ".txt")
            try {
                marker.toFile().setText(latest + "\n", "UTF-8")
                Files.move(marker, versionFile.toPath(), StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING)
            } finally { Files.deleteIfExists(marker) }
            System.err.println("[VDJDB-update] Done, using " + latest)
            true
        } finally { staging.deleteDir() }
    }

    private static final List<String> DATABASE_FILES =
            ["vdjdb.txt", "vdjdb.meta.txt", "vdjdb.slim.txt", "vdjdb.slim.meta.txt"]

    private static final List<String> RELEASE_FILES = DATABASE_FILES +
            ["vdjdb_full.txt", "cluster_members.txt", "motif_pwms.txt",
             "cluster_members_tcremp.txt", "motif_pwms_tcremp.txt", "vdjdb_summary_embed.html", "LICENSE"]

    private static String releaseDirectory(String url) {
        "vdjdb-" + MessageDigest.getInstance("SHA-256").digest(url.getBytes("UTF-8")).encodeHex().toString()
    }

    static File databaseDirectory(File home = new File(HOME_DIR)) {
        def marker = new File(home, "latest-version.txt")
        def url = marker.isFile() ? marker.readLines("UTF-8").find { it.trim() } : null
        def installed = url ? new File(home, releaseDirectory(url)) : home
        installed.isDirectory() ? installed : home // historical flat installations
    }

    private static File tableFile(File file) {
        def gzip = new File(file.path + ".gz")
        if (file.isFile() && gzip.isFile()) throw new IOException("Ambiguous table: " + file)
        !file.isFile() && gzip.isFile() ? gzip : file
    }

    static InputStream tableAsStream(File file) {
        def resolved = tableFile(file)
        def stream = new FileInputStream(resolved)
        if (!resolved.name.endsWith(".gz")) return stream
        try { new GZIPInputStream(stream) }
        catch (Exception error) { stream.close(); throw error }
    }

    /**
     * Gets an input stream for a given resource
     * @param resourceName resource file name
     * @return input stream
     */
    static InputStream resourceAsStream(String resourceName) {
        Util.class.classLoader.getResourceAsStream(resourceName)
    }

    /**
     * Converts a given string to amino acid sequence
     * @param aaSeq amino acid string
     * @return binary AminoAcidSequence object or null if failed to convert
     */
    static AminoAcidSequence convert(String aaSeq) {
        aaSeq = aaSeq.trim()
        if (aaSeq.length() > 0 && aaSeq =~ /^[FLSYCWPHQRIMTNKVADEG]+$/)
            return new AminoAcidSequence(aaSeq)

        System.err.println("Error converting '$aaSeq' to amino acid sequences, entry will be skipped from search")
        null
    }


    static void sout(String header, String message) {
        println "[${new Date()} $header] $message"
    }

    static void error(String message) {
        sout("ERROR", message)
        System.exit(1)
    }

    static void info(String message) {
        sout("INFO", message)
    }
}
