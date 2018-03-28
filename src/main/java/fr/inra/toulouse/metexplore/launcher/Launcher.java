package fr.inra.toulouse.metexplore.launcher;

import fr.inra.toulouse.metexplore.io.WritingComportment;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.File;

import static java.lang.System.exit;

public abstract class Launcher implements WritingComportment{
    @Option(name = "-h", aliases = "-help", usage = "Prints this help.")
    protected boolean phelp = false;

    @Option(name = "-v", aliases = "-version", usage = "Prints the current version of the program.")
    protected boolean version;

    @Option(name = "-gal", aliases = "-galaxy", usage = "For galaxy compliance: formatting pathway output and creating a new one containing log information.")
    protected String galaxyFile = "";

    protected String mappingWarnings = "[WARNING] By default, a mapping has been set with the name and the SBML id respectively on the 1st and the 2nd column of your dataset.\n" +
            "[WARNING] Other mapping available: ChEBI, InChI, InChIKey, SMILES, CSID, PubChem and HMDB (check -help).\n";

    protected static File logFile;
    protected String logContent;

    protected static long startTime;

    public static void timeCalculation(long elapsedTime) {
        long min = elapsedTime / 60000000000L;
        long sec = elapsedTime / 1000000000L - (min * 60L);
        System.out.println("Time to run the process : " + min + "min " + sec + "s");

    }

    public void printError(CmdLineParser parser, CmdLineException e){
        System.err.println(e.getMessage());
        System.err.println("Options:");
        parser.printUsage(System.err);
        exit(1);
    }

    public void printInfo(CmdLineParser parser){

        if (this.phelp) {
            System.out.println("Options:");
            parser.printUsage(System.out);
            exit(0);
        }

        if (this.version) {
            System.out.println("Version: 2.0");
            exit(0);
        }
    }

    public void writeLog(String warning){
        this.logContent = writeLog(this.logContent,warning);
    }
}