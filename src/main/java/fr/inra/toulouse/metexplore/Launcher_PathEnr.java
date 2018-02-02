package fr.inra.toulouse.metexplore;

import java.io.IOException;
import java.util.regex.Pattern;

import parsebionet.biodata.BioNetwork;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import static java.lang.System.exit;
import org.kohsuke.args4j.CmdLineException;

public class Launcher_PathEnr {

    @Option(name="-h", usage="Prints this help.")
    public boolean phelp = false;

    @Option(name="-o1", usage="Output file name for mapping result (by default: NONE).")
    String outFileMapping = "";

    @Option(name="-gal", usage="Output file name for general information resulting from mapping and pathway enrichment results (by default: NONE).")
    String galaxy;

    @Option(name="-o2", usage="Output file name for pathway enrichment result (by default: pathwayEnrichment.tsv).")
    String outFilePathEnr = "pathwayEnrichment.tsv";

    @Option(name="-s", usage="SBML file name.")
    public String sbml = "data/recon2.v03_ext_noCompartment_noTransport_v2.xml";

    @Option(name="-i", usage="[REQUIRED] Input file containing a fingerprint in tsv file format.")
    public String inFileFingerprint ;

    @Option(name="-name", usage="Number of the file's column containing the metabolite name.")
    public int nameColumn = 1;

    @Option(name="-f", usage="Number of the filtered column (by default: 0 for none)")
    public int colFiltered = -1;

    @Option(name="-inchi", usage="Number of the file's column containing the InChI data (by default: 5; 0 for none).")
    public int inchiColumn = -1;

    @Option(name="-chebi", usage="Number of the file's column containing the CHEBI data (by default: 0 for none).")
    public int chebiColumn = -1;

    @Option(name="-id", usage="Number of the file's column containing the metabolite identifier (by default: 0 for none).")
    public int idSBMLColumn = 2;

    @Option(name="-smiles", usage="Number of the file's column containing the SMILES data (by default: 0 for none).")
    public int smilesColumn = -1;

    @Option(name="-pubchem", usage="Number of the file's column containing the PubChem identifier (by default: 0 for none).")
    public int pubchemColum = -1;

    @Option(name="-inchikey", usage="Number of the file's column containing the InChIKey (by default: 0 for none).")
    public int inchikeysColumn = -1;

    @Option(name="-kegg", usage="Number of the file's column containing the KEGG identifier (by default: 0 for none).")
    public int keggColumn = -1;

    @Option(name="-hmdb", usage="Number of the file's column containing the HMD identifier (by default: 0 for none).")
    public int hmdColumn = -1;

    @Option(name="-chemspider", usage="Number of the file's column containing the ChemSpider identifier (by default: 0 for none).")
    public int chemspiderColumn = -1;

    @Option(name="-weight", usage="Number of the file's column containing the weigth of the metabolites (by default: 0 for none).")
    public int weightColumn = -1;

    @Option(name="-l", usage="List containing the number - separated by comma without blank spaces - of the InChi's layer concerned by the mapping (by default: c,h; for a mapping including all the layers, enter c,h,q,p,b,t,i,f,r; for a mapping on formula layer only, enter the -l option with no parameter).")
    public String inchiLayers = "c,h";

    @Option(name="h=F", usage="Activate this option if the fingerprint dataset contains no header.")
    public String header;

    @Option(name="-sep", usage="Activate this option if the fingerprint dataset contains no header.")
    public String separator = "\t";

    public void timeCalculation(long elapsedTime){
        long min = elapsedTime / 60000000000L;
        long sec = elapsedTime / 1000000000L - (min * 60L);
        System.out.println("Time to run the process : " + min + "min " + sec + "s");

    }

    @SuppressWarnings("deprecation")
    public static void main(String[] args) {

        long startTime = System.nanoTime();
        Boolean ifHeader = true;//Take account of the header
        Boolean ifGalaxy = false;//Galaxy compliance
        Launcher_PathEnr launch = new Launcher_PathEnr();
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
            if(launch.inFileFingerprint==null){
                throw new CmdLineException("-i parameter required");
            }

            if (launch.chebiColumn < 1 && launch.inchiColumn < 1 && launch.idSBMLColumn < 1 &&
                    launch.smilesColumn < 1 && launch.pubchemColum < 1 && launch.inchikeysColumn < 1
                    && launch.keggColumn < 1 && launch.hmdColumn < 1 && launch.chemspiderColumn < 1
                    && launch.weightColumn < 1){
                throw new CmdLineException("Mapping parameters cannot be all set at < 1. Choose at less one criterion.");
            }

            if (!Pattern.matches("([chqpbtifr],)*[chqpbtifr]", launch.inchiLayers)) {
                throw new CmdLineException("-l parameter badly formatted");
            }

            //Personalised error print for help
        } catch (CmdLineException e) {
            if(e.getMessage().equals("No argument is allowed: h=F")) {
                ifHeader = false;
            }else if(e.getMessage().equals("Option \"-gal\" takes an operand")) {
                ifGalaxy = true;
            }else if(e.getMessage().equals("Option \"-l\" takes an operand")){
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
        BioNetwork network = (new JSBML2Bionetwork4Galaxy(launch.sbml)).getBioNetwork();
        int[] mappingColumns = {(launch.idSBMLColumn-1), (launch.inchiColumn-1), (launch.chebiColumn-1),
                (launch.smilesColumn-1), (launch.pubchemColum-1), (launch.inchikeysColumn-1),
                (launch.keggColumn-1), (launch.hmdColumn-1), (launch.chemspiderColumn-1), (launch.weightColumn-1)};

        try{
            Fingerprint fingerprint = new Fingerprint(launch.inFileFingerprint,ifHeader, launch.separator, (launch.nameColumn-1),
                    mappingColumns, (launch.colFiltered-1));
            Mapping mapping = new Mapping(network, fingerprint.list_metabolites, inchiLayers,
                    launch.outFileMapping, ifGalaxy);
            PathwayEnrichment pathEnr = new PathwayEnrichment(network, mapping.list_mappedMetabolites,
                    launch.outFilePathEnr,ifGalaxy);
        }
        catch (IOException e){
            e.printStackTrace();
        }
        launch.timeCalculation(System.nanoTime() - startTime);
    }
}
