package fr.inra.toulouse.metexplore;

import java.io.IOException;
import java.util.regex.Pattern;

import parsebionet.biodata.BioNetwork;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import static java.lang.System.exit;
import org.kohsuke.args4j.CmdLineException;

public class Launcher_MetExplore4Galaxy {

    @Option(name="-h", usage="Prints this help.")
    public boolean phelp = false;

    @Option(name="-o1", usage="Output file name for mapping result (by default: mapping.tsv).")
    String outFile1 = "mapping.tsv";

    @Option(name="-o2", usage="Output file name for pathway enrichment result (by default: pathwayEnrichment.tsv).")
    String outFile2 = "pathwayEnrichment.tsv";

    @Option(name="-o3", usage="Output file name for general information resulting from mapping and pathway enrichment results (by default: info.txt).")
    String outFile3 = "info.txt";

    @Option(name="-s", usage="Sbml file name.")
    public String sbml = "data/recon2.v03_ext_noCompartment_noTransport_v2.xml";

    @Option(name="-i", usage="[REQUIRED] Input file in tsv file format.")
    public String inFile ;

    @Option(name="-f", usage="Number of the filtered column (by default: 0 for none)")
    public int colFiltered = -1;

    @Option(name="-inchi", usage="Number of the file's column containing the InChI data (by default: 5; 0 for none).")
    public int inchiColumn = -1;

    @Option(name="-chebi", usage="Number of the file's column containing the CHEBI data (by default: 0 for none).")
    public int chebiColumn = -1;

    @Option(name="-id", usage="Number of the file's column containing the metabolite identifier in the SBML file (by default: 0 for none).")
    public int idSBMLColumn = 2;

    @Option(name="-l", usage="List containing the number - separated by comma without blank spaces - of the InChi's layer concerned by the mapping (by default: c,h; for a mapping including all the layers, enter c,h,q,p,b,t,i,f,r; for a mapping on formula layer only, enter the -l option with no parameter).")
    public String inchiLayers = "c,h";


    @SuppressWarnings("deprecation")
    public static void main(String[] args) {

        long startTime = System.nanoTime();
        Launcher_MetExplore4Galaxy launch = new Launcher_MetExplore4Galaxy();
        CmdLineParser parser = new CmdLineParser(launch);

        try {
            parser.parseArgument(args);

            //Print help
            if (launch.phelp) {
                System.out.println("Options:");
                parser.printUsage(System.out);
                exit(0);
            }

            //Error messages for bad parameters
            if(launch.inFile==null){
                throw new CmdLineException("-i parameter required");
            }

            if (launch.chebiColumn < 1 && launch.inchiColumn < 1 && launch.idSBMLColumn < 1) {
                throw new CmdLineException("-chebi, -inchi and idSBML parameters cannot be all setted at < 1. Choose at less one mapping criterion.");
            }

            if (!Pattern.matches("([chqpbtifr],)*[chqpbtifr]", launch.inchiLayers)) {
                throw new CmdLineException("-l parameter badly formatted");
            }

        //Personalised error print for help
        } catch (CmdLineException e) {
            if(e.getMessage().equals("Option \"-l\" takes an operand")){
                launch.inchiLayers="";
            }else {
                System.err.println(e.getMessage());
                System.err.println("Options:");
                parser.printUsage(System.err);
                exit(1);
            }
        }

        //Regex for inchiLayers parameter
        String[] inchiLayers = launch.inchiLayers.replaceAll(" ","").split(",");
        //Extract SBML
        BioNetwork bionet = (new JSBML2Bionetwork4Galaxy(launch.sbml)).getBioNetwork();

        try{
            MetExplore4Galaxy met = new MetExplore4Galaxy(bionet,launch.inFile,launch.outFile1,launch.outFile2,launch.outFile3,(launch.chebiColumn -1),(launch.inchiColumn -1),(launch.idSBMLColumn -1),(launch.colFiltered -1),inchiLayers);
            met.exec(startTime);
        }
        catch (IOException e){
            e.printStackTrace();
        }
    }
}