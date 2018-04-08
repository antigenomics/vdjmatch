package com.antigenomics.vdjdb

import com.antigenomics.vdjtools.annotate.Annotate
import com.antigenomics.vdjtools.misc.ExecUtil

import java.util.jar.JarFile

def getVersion = {
    (getClass().classLoader.findResource(JarFile.MANIFEST_NAME).text =~
            /Implementation-Version: (.+)/)[0][1]
}

def printHelp = {
    def version = getVersion()

    println "VDJmatch V$version"
    println ""
    println "Run as \$java -jar vdjmatch-${version}.jar ROUTINE_NAME arguments"
    println ""
    println "Supported ROUTINE_NAME:"
    println "match   - matches against a TCR specificity database"
    println "cluster - cluster TCR clonotypes in a sample based on sequence similarity"
    println "update  - update you local copy of VDJdb specificity database to latest version"
}

def getScript = { String scriptName ->
    switch (scriptName.toUpperCase()) {
        case "MATCH":
            return new RunMatch()
        case "CLUSTER":
            return new RunCluster()
        case "UPDATE":
            return new RunUpdate()
        case "ANNOTATE":
            return new Annotate()

        case "-H":
        case "H":
        case "-HELP":
        case "HELP":
        case "":
            printHelp()
            println ""
            System.exit(0)
            break

        case "DUMMY":
            System.exit(0)
            break

        default:
            printHelp()
            println ""
            println "Unknown routine name $scriptName"
            System.exit(0)
    }
}

if (args.length == 0) {
    printHelp()
} else {
    def script = getScript(args[0])
    try {
        ExecUtil.run(script, args.size() > 1 ? args[1..-1] : [""])
    } catch (Exception e) {
        def version = getVersion()
        println "[CRITICAL ERROR] ${e.toString()}, see _vdjmatch_error.log for details"
        new File("_vdjmatch_error.log").withWriterAppend { writer ->
            writer.println("[${new Date()} BEGIN]")
            writer.println("[Script]")
            writer.println(args[0])
            writer.println("[CommandLine]")
            writer.println("executing vdjmatch-${version}.jar ${args.join(" ")}")
            writer.println("[Message]")
            writer.println(e.toString())
            writer.println("[StackTrace-Short]")
            writer.println(e.stackTrace.findAll { it.toString().contains("com.antigenomics.vdjdb") }.join("\n"))
            writer.println("[StackTrace-Full]")
            e.printStackTrace(new PrintWriter(writer))
            writer.println("[END]")
        }
        System.exit(1)
    }
}