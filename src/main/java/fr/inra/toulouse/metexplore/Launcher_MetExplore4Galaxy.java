package fr.inra.toulouse.metexplore;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPathway;
import parsebionet.biodata.BioPhysicalEntity;
import parsebionet.io.JSBMLToBionetwork;
import parsebionet.statistics.PathwayEnrichment;

import java.io.IOException;

import java.util.ArrayList;
import java.util.Set;
import java.util.HashMap;

public class Launcher_MetExplore4Galaxy {

    @Option(name="-h", usage="Prints this help.")
    public boolean phelp = false;

    @Option(name="-o", usage="Output file name  (by default: output.tsv).")
    String outFile = "output.tsv";

    @Option(name="-sbml", usage="Sbml file name.")
    public String sbml = "recon2.v03_ext_noCompartment_noTransport_v2.xml";

    @Option(name="-i", usage="[Required] Input file in tsv file format.")
    //public String inFile = "Galaxy15-[Biosigner_Multivariate_Univariate_Multivariate_variableMetadata.tsv].tabular";
    public String inFile = "sacurineVariableMetadataEnhanced.tsv";
    //public String inFile ;

    @Option(name="-f", usage="Choose to filter some columns  (by default: false).")
    public Boolean filter = false;

    @Option(name="-fcol", usage="Number of the filtered column (if the filter option is choosen)")
    public int colFiltered = 0;

    @Option(name="-inchi", usage="Number of the file's column containing the InChI data (first case: 0; by default: 4).")
    public int inchiColumn = 5;

    @Option(name="-form", usage="List containing the number (separated by comma; ex: 0,1,2,3,4; first case: 0) of the InChi's layer concerned by the mapping (by default: all of them).")//TODO: implements this functionality
    public String formulaColumn = "1,2,3,4,5,6,7";//TODO: catch the elements in an array of Integer with formulaColumn.split(",")

    @SuppressWarnings("deprecation")
    public static void main(String[] args) {

        //Parameters
        long startTime = System.nanoTime();
        Launcher_MetExplore4Galaxy launch = new Launcher_MetExplore4Galaxy();
        MetExplore4Galaxy met = new MetExplore4Galaxy();
        CmdLineParser parser = new CmdLineParser(launch);


        //CmdParsing
        try {
            parser.parseArgument(args);

            if(launch.phelp){
                System.out.println("Application Usage:");
                parser.printUsage(System.out);
                System.exit(0);
            }

        } catch (CmdLineException e) {

            System.out.println("Application Usage:");
            parser.printUsage(System.out);
            System.exit(0);

        }


        //SBML parsing
        JSBMLToBionetwork jsbml = new JSBMLToBionetwork(launch.sbml);//TODO: test function SBML
        BioNetwork bionet = jsbml.getBioNetwork();

        try{

            //Mapping
            HashMap <String, String[]> parsedFile = met.extractData(launch.inFile, launch.filter, (launch.colFiltered -1));
            Set<BioPhysicalEntity> map = met.mapping(bionet, parsedFile, (launch.inchiColumn -1));

            //PathwayEnrichment
            PathwayEnrichment enr = new PathwayEnrichment(bionet, map);
            ArrayList<HashMap<BioPathway, Double>> resultList = new ArrayList<HashMap<BioPathway, Double> >();

            resultList.add(enr.computeEnrichment());
            for (int i = 0; i < 3; i++) {
                resultList.add(enr.computeEnrichment(i));
            }

            System.err.println(resultList.get(0).size() + " pathways are concerned among the network (on a total of " + bionet.getPathwayList().size() + ").");
            met.writeOutput(resultList, map, launch.outFile);
       }
        catch (IOException e2){
            e2.printStackTrace();
        }


        long elapsedTime = System.nanoTime() - startTime;
        if (elapsedTime / 1000000000L == 0L) {
            System.err.println("Time to run the process in miliseconde :" + elapsedTime / 1000000L);
        } else {
            System.err.println("Time to run the process in seconde :" + elapsedTime / 1000000000L);
        }

    }
}

