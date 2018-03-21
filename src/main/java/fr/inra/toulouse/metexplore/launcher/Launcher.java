package fr.inra.toulouse.metexplore.launcher;

import fr.inra.toulouse.metexplore.WritingComportment;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import static java.lang.System.exit;

public abstract class Launcher {
    @Option(name = "-h", aliases = "--help", usage = "Prints this help.")
    protected boolean phelp = false;

    @Option(name = "-v", aliases = "--version", usage = "Prints the current version of the program.")
    protected boolean version;

    @Option(name = "-gal", aliases = "--galaxy", usage = "For galaxy compliance: formatting pathway output and creating a new one containing log information.")
    protected String galaxy = "";

    protected String mappingWarnings = "[WARNING] By default, a mapping has been set with the name and the SBML id respectively on the 1st and the 2nd column of your dataset.\n" +
            "[WARNING] Other mapping available: ChEBI, InChI, InChIKey, SMILES, CSID, PubChem and HMDB (check --help).\n";

    protected WritingComportment write = new WritingComportment(this.galaxy);

    protected static long startTime = System.nanoTime();

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
            System.out.println("Version: 1.1");
            exit(0);
        }
    }

    /* @SuppressWarnings("deprecation")
   public static void exec(String[] args) {

        long startTime = System.nanoTime();
        Launcher launch = new Launcher();
        CmdLineParser parser = new CmdLineParser(launch);

        try {
            parser.parseArgument(args);

        launch.timeCalculation(System.nanoTime() - startTime);
    }*/


}