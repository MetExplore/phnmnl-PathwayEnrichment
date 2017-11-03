package fr.inra.toulouse.metexplore;

import java.io.IOException;
import java.util.Set;
import java.util.HashMap;
import java.util.regex.Pattern;

import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPhysicalEntity;
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
    String outFile3 = "";

    @Option(name="-s", usage="Sbml file name.")
    public String sbml = "data/recon2.v03_ext_noCompartment_noTransport_v2.xml";

    @Option(name="-i", usage="[REQUIRED] Input file in tsv file format.")
    public String inFile ;

    @Option(name="-f", usage="Number of the filtered column (by default: 0 for none)")
    public int colFiltered = -1;

    @Option(name="-inchi", usage="Number of the file's column containing the InChI data (by default: 5; 0 for none).")
    public int inchiColumn = 5;

    @Option(name="-chebi", usage="Number of the file's column containing the chebi data (by default: 0 for none).")
    public int chebiColumn = -1;

    @Option(name="-l", usage="List containing the number - separated by comma without blank spaces - of the InChi's layer concerned by the mapping (by default: c,h; for a mapping including all the layers, enter c,h,q,p,b,t,i,f,r; for a mapping on formula layer only, enter the -l option with no parameter).")
    public String inchiLayers = "c,h";


    @SuppressWarnings("deprecation")
    public static void main(String[] args) {

        //Parameters
        long startTime = System.nanoTime();
        Launcher_MetExplore4Galaxy launch = new Launcher_MetExplore4Galaxy();
        MetExplore4Galaxy met = new MetExplore4Galaxy(launch.outFile3);
        CmdLineParser parser = new CmdLineParser(launch);

        //CmdParsing
        try {
            parser.parseArgument(args);

            if (launch.phelp) {
                System.out.println("Options:");
                parser.printUsage(System.out);
                exit(0);
            }

            if(launch.inFile==null){
                throw new CmdLineException("-i parameter required");
            }

            if (launch.chebiColumn < 1 && launch.inchiColumn < 1) {
                throw new CmdLineException("-chebi and -inchi parameters cannot be both setted at < 1");
            }

            if (!Pattern.matches("([chqpbtifr],)*[chqpbtifr]", launch.inchiLayers)) {
                throw new CmdLineException("-l parameter badly formatted");
            }

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
        String[] inchiLayers = launch.inchiLayers.replaceAll(" ","").split(",");
        BioNetwork bionet = (new JSBML2Bionetwork4Galaxy(launch.sbml)).getBioNetwork();

        try{

            //Pathway Enrichment
            HashMap <String, String[]> parsedFile = met.extractData(launch.inFile, (launch.colFiltered -1));
            Set<BioPhysicalEntity> map = met.mapping(bionet, launch.outFile1,launch.outFile3, parsedFile, (launch.chebiColumn -1), (launch.inchiColumn -1), inchiLayers);
            met.writeOutput(met.pathwayEnrichment (bionet, launch.outFile3, map), map, launch.outFile2);

        } catch (IOException e2){
            e2.printStackTrace();
        }

        met.runTime(System.nanoTime() - startTime);
    }
}