package fr.inra.toulouse.metexplore;

import java.io.IOException;
import java.util.Set;
import java.util.HashMap;
import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPhysicalEntity;
import parsebionet.io.JSBMLToBionetwork;
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

    @Option(name="-s", usage="Sbml file name.")
    public String sbml = "data/recon2.v03_ext_noCompartment_noTransport_v2.xml";

    @Option(name="-i", usage="[Required] Input file in tsv file format.")
    public String inFile ;

    @Option(name="-f", usage="Number of the filtered column (by default: 0 for none)")
    public int colFiltered = -1;

    @Option(name="-inchi", usage="Number of the file's column containing the InChI data (by default: 5).")
    public int inchiColumn = 5;

    @Option(name="-chebi", usage="Number of the file's column containing the chebi data (by default: 2).")
    public int chebiColumn = 2;

    @Option(name="-l", usage="List containing the number - separated by comma without blank spaces - of the InChi's layer concerned by the mapping (by default: c,h; for all layers selection, enter c,h,q,p,b,t,i,f,r).")
    public String inchiLayers = "c,h";

    @Option(name="-c", usage="Test correction (by default: Bonferoni; other options: enter 2 for Benjamini-Hochberg correction and 3 for Holm-Bonferroni one).")
    public int correction = 1;

    @SuppressWarnings("deprecation")
    public static void main(String[] args) {

        //Parameters
        long startTime = System.nanoTime();
        Launcher_MetExplore4Galaxy launch = new Launcher_MetExplore4Galaxy();
        MetExplore4Galaxy met = new MetExplore4Galaxy();
        CmdLineParser parser = new CmdLineParser(launch);
        String[] inchiLayers = launch.inchiLayers.split(",");


        //CmdParsing
        try {
            parser.parseArgument(args);

            if(launch.phelp){
                System.out.println("Options:");
                parser.printUsage(System.out);
                exit(0);
            }

        } catch (CmdLineException e) {

            System.out.println("Options:");
            parser.printUsage(System.out);
            exit(0);

        }

        BioNetwork bionet = (new JSBMLToBionetwork(launch.sbml)).getBioNetwork();

        try{

            //Pathway Enrichment
            HashMap <String, String[]> parsedFile = met.extractData(launch.inFile, (launch.colFiltered -1));
            Set<BioPhysicalEntity> map = met.mapping(bionet, launch.outFile1, parsedFile, (launch.chebiColumn -1), (launch.inchiColumn -1), inchiLayers);
            met.writeOutput(met.pathwayEnrichment (bionet, map, launch.correction), map, launch.outFile2);

        } catch (IOException e2){
            e2.printStackTrace();
        }

        met.runTime(System.nanoTime() - startTime);
    }
}

