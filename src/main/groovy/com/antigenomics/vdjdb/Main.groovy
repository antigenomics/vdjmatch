/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */


package com.antigenomics.vdjdb

import com.antigenomics.vdjdb.impl.ClonotypeDatabase
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.ExecUtil

import static com.antigenomics.vdjdb.Util.resourceAsStream

def DEFAULT_PARAMERES = "2,1,1,2"
def cli = new CliBuilder(usage: "vdjdb [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.D(longOpt: "database", argName: "string", args: 1, "Path to an external database file.")
cli._(longOpt: "filter", argName: "logical expression(__field__,...)", args: 1,
        "Logical filter evaluated for database columns. Supports Regex, .contains(), .startsWith(), etc.")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.p(longOpt: "parameters", argName: "s,i,d,t", args: 1,
        "CDR3 sequence search parameters: " +
                "allowed number of substitutions (s), insertions (i), deletions (d) and total number of mutations. " +
                "[default=$DEFAULT_PARAMERES]")
cli.v(longOpt: "v-match", "Require V segment matching.")
cli.j(longOpt: "j-match", "Require J segment matching.")
cli.c("Compressed output")

def opt = cli.parse(args)

if (opt == null || opt.h || opt.arguments().size() == 0) {
    cli.usage()
    System.exit(1)
}

// Check if metadata is provided

def metadataFileName = opt.m

if (metadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 2) {
    if (metadataFileName)
        println "[ERROR] Only output prefix should be provided in case of -m"
    else
        println "[ERROR] At least 1 sample files should be provided if not using -m"
    cli.usage()
    System.exit(1)
}

// Remaining arguments

def dbPrefix = (String) (opt.D ?: null), p = (opt.p ?: DEFAULT_PARAMERES).split(",").collect { it.toInteger() },
    compress = (boolean) opt.c,
    vMatch = (boolean) opt."v-match", jMatch = (boolean) opt."j-match",
    filter = (String) (opt.'filter' ?: null),
    outputFileName = opt.arguments()[-1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

println "[${new Date()} $scriptName] Loading database..."

ClonotypeDatabase database

def metaStream = dbPrefix ? new FileInputStream("${dbPrefix}.meta") : resourceAsStream("vdjdb_legacy.meta"),
    dataStream = dbPrefix ? new FileInputStream("${dbPrefix}.txt") : resourceAsStream("vdjdb_legacy.txt")

database = new ClonotypeDatabase(metaStream, vMatch, jMatch, p[0], p[1], p[2], p[3])
database.addEntries(dataStream, filter)

println "[${new Date()} $scriptName] Finished.\n $database"

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading sample(s)..."

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-2])

println "[${new Date()} $scriptName] ${sampleCollection.size()} sample(s) to process."

//
// Main loop
//

println "[${new Date()} $scriptName] Annotating sample(s) & writing results."

def sw = new SampleWriter(compress)

//new File(ExecUtil.formOutputPath(outputFileName, "annot", "summary")).withPrintWriter { pwSummary ->
//    def header = "$MetadataTable.SAMPLE_ID_COLUMN\t" +
//            sampleCollection.metadataTable.columnHeader + "\tdatabase\tfilter"
//    pwSummary.println(header)

sampleCollection.eachWithIndex { Sample sample, int ind ->
    def sampleId = sample.sampleMetadata.sampleId

    def results = database.search(sample)

    def writer = sw.getWriter(ExecUtil.formOutputPath(outputFileName, sampleId, "annot"))

    writer.println(sw.header + "\tpenalty\t" + database.header)
    results.each { result ->
        result.value.each {
            writer.println(sw.getClonotypeString(result.key) + "\t" +
                    it.result.penalty + "\t" + it.row.toString())
        }
    }

    writer.close()

    println "[${new Date()} $scriptName] ${ind + 1} sample(s) done."
}
//}

println "[${new Date()} $scriptName] Finished."