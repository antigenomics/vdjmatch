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
        def version

        def versionFile = new File(HOME_DIR + "/latest-version.txt")
        if (!versionFile.exists()) {
            System.err.println("[VDJDB-update] No database is present (running for the first time, I assume), downloading it. " +
                    "NOTE: No automatic update will be performed next time VDJdb is run - use Update command to check for database updates.")
            version = ""
        } else {
            version = versionFile.readLines()[0]
        }

        if (!versionFile.exists() || updateIfNewer) {
            def latestVersion = new URL("https://raw.githubusercontent.com/antigenomics/vdjdb-db/master/latest-version.txt").openStream().readLines()[0]

            if (version != latestVersion) {
                versionFile << latestVersion
            }

            def dbZip = new File(HOME_DIR + "/vdjdb.zip")
            def out = new BufferedOutputStream(new FileOutputStream(dbZip))
            out << new URL(latestVersion).openStream()
            out.close()

            def ant = new AntBuilder()

            ant.unzip(src: dbZip.absolutePath,
                    dest: HOME_DIR,
                    overwrite: "true")

            System.err.println("[VDJDB-update] Done, you are now using VDJDB-V${latestVersion.split("/")[-2]}")

            return true
        }

        false
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
