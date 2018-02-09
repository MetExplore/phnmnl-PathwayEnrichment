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
    protected boolean phelp = false;

    /******FILES PARAMETERS*****/

    @Option(name="-i", usage="[REQUIRED] Input file containing a fingerprint in tsv file format.")
    protected String inFileFingerprint ;

    @Option(name="-o1", usage="Output file name for mapping result (by default: disabled).")
    protected String outFileMapping = "mapping.tsv";

    @Option(name="-o2", usage="Output file name for pathway enrichment result (by default: pathwayEnrichment.tsv).")
    protected String outFilePathEnr = "pathwayEnrichment.tsv";

    @Option(name="-s", usage="SBML file name (by default: Recon 2v02).")
    protected String sbml = "data/recon2.02_without_compartment.xml.xml";

    /******PARSING PARAMETERS*****/

    @Option(name="-header", usage="Activate this option if the fingerprint dataset contains no header (by default: disabled).")
    protected boolean ifNoHeader = false;

    @Option(name="-sep", usage="Character used as separator in the dataset (by default: \\t for tab).")
    protected String separator = "\t";

    @Option(name="-f", usage="Number of the filtered column (by default: disabled)")
    protected int colFiltered = -1;

    @Option(name="-gal", usage="Formating output in a galaxy compliant way (by default: disabled).")
    protected boolean ifGalaxy = false;

    /*****MAPPING PARAMETERS*****/

    @Option(name = "-t", usage = "Type of biological object selected : 1 for metabolites, 2 for reactions, 3 for pathway, 4 for enzyme, 5 for protein, 6 for gene (by default: metabolites).")
    protected int bioEntityType = 1;

    @Option(name="-name", usage="Number of the file's column containing the metabolite name (by default: 1st column).")
    protected int nameColumn = 1;

    @Option(name="-l", usage="List containing the number - separated by comma without blank spaces - of the InChi's layer concerned by the mapping" +
            " (by default: c,h; for a mapping including all the layers, enter c,h,q,p,b,t,i,f,r; for a mapping on formula layer only, enter the -l option with no parameter).")
    protected String inchiLayers = "c,h";

    @Option(name="-inchi", usage="Number of the file's column containing the InChI data (by default: 5th column).")
    protected int inchiColumn = 5;

    @Option(name="-chebi", usage="Number of the file's column containing the CHEBI data (by default: disabled).")
    protected int chebiColumn = -1;

    @Option(name="-id", usage="Number of the file's column containing the metabolite identifier (by default: disabled).")
    protected int idSBMLColumn = -1;

    @Option(name="-smiles", usage="Number of the file's column containing the SMILES data (by default: disabled).")
    protected int smilesColumn = -1;

    @Option(name="-pubchem", usage="Number of the file's column containing the PubChem identifier (by default: disabled).")
    protected int pubchemColum = -1;

    @Option(name="-inchikey", usage="Number of the file's column containing the InChIKey (by default: disabled).")
    protected int inchikeysColumn = -1;

    @Option(name="-kegg", usage="Number of the file's column containing the KEGG identifier (by default: disabled).")
    protected int keggColumn = -1;

    @Option(name="-hmdb", usage="Number of the file's column containing the HMD identifier (by default: disabled).")
    protected int hmdColumn = -1;

    @Option(name="-chemspider", usage="Number of the file's column containing the ChemSpider identifier (by default: disabled).")
    protected int chemspiderColumn = -1;

    @Option(name="-weight", usage="Number of the file's column containing the weigth of the metabolites (by default: disabled).")
    protected int weightColumn = -1;


    public void timeCalculation(long elapsedTime){
        long min = elapsedTime / 60000000000L;
        long sec = elapsedTime / 1000000000L - (min * 60L);
        System.out.println("Time to run the process : " + min + "min " + sec + "s");

    }

    @SuppressWarnings("deprecation")
    public static void main(String[] args) {

        long startTime = System.nanoTime();
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

            if (launch.nameColumn < 1 && launch.chebiColumn < 1 && launch.inchiColumn < 1 && launch.idSBMLColumn < 1 &&
                    launch.smilesColumn < 1 && launch.pubchemColum < 1 && launch.inchikeysColumn < 1
                    && launch.keggColumn < 1 && launch.hmdColumn < 1 && launch.chemspiderColumn < 1
                    && launch.weightColumn < 1){
                throw new CmdLineException("Mapping parameters cannot be all set at < 1. Choose at less one criterion.");
            }

            if (!Pattern.matches("([chqpbtifr],)*[chqpbtifr]", launch.inchiLayers)) {
                throw new CmdLineException("-l parameter badly formatted");
            }

            if (launch.bioEntityType < 1 || launch.bioEntityType > 6) {
                throw new CmdLineException("Type of biological object must be between 1 and 6.");
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
        BioNetwork network = (new JSBML2Bionetwork4Galaxy(launch.sbml)).getBioNetwork();
        int[] mappingColumns = {(launch.idSBMLColumn-1), (launch.inchiColumn-1), (launch.chebiColumn-1),
                (launch.smilesColumn-1), (launch.pubchemColum-1), (launch.inchikeysColumn-1),
                (launch.keggColumn-1), (launch.hmdColumn-1), (launch.chemspiderColumn-1), (launch.weightColumn-1)};

        try{
            Fingerprint fingerprint = new Fingerprint(launch.inFileFingerprint,launch.ifNoHeader, launch.separator, (launch.nameColumn-1),
                    mappingColumns, (launch.colFiltered-1));
            Mapping mapping = new Mapping(network, fingerprint.list_entities, inchiLayers,
                    launch.outFileMapping, launch.ifGalaxy, launch.bioEntityType);
            PathwayEnrichment pathEnr = new PathwayEnrichment(network, fingerprint.list_entities, mapping.list_mappedEntities,
                    launch.outFilePathEnr,launch.ifGalaxy, launch.bioEntityType);
        }
        catch (IOException e){
            e.printStackTrace();
        }
        launch.timeCalculation(System.nanoTime() - startTime);
    }
}