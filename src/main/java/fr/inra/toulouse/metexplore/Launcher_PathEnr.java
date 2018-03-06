package fr.inra.toulouse.metexplore;

import java.io.IOException;
import java.util.regex.Pattern;

import parsebionet.biodata.BioNetwork;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import static java.lang.System.exit;
import org.kohsuke.args4j.CmdLineException;

public class Launcher_PathEnr {

    @Option(name="-h", aliases="--help", usage="Prints this help.")
    protected boolean phelp = false;

    @Option(name="-v", aliases="--version", usage="Prints the current version of the program.")
    protected boolean version ;

    /******FILES PARAMETERS*****/

    @Option(name="-i", aliases="--inFile", usage="[REQUIRED] Input file containing a fingerprint (in tsv file format).")
    protected String inFileFingerprint ;

    @Option(name="-o1", aliases="--outMap", usage="Output file name for mapping result.")
    protected String outFileMapping = "mapping.tsv";

    @Option(name="-o2", aliases="--outPath", usage="Output file name for pathway enrichment result (by default: pathwayEnrichment.tsv).")
    protected String outFilePathEnr = "pathwayEnrichment.tsv";

    @Option(name="-s", aliases="--sbml", usage="SBML file name (by default: Recon 2v02).")
    protected String sbml = "data/recon2.02_without_compartment.xml";

    /******PARSING PARAMETERS*****/

    @Option(name="--header", usage="Activate this option if the fingerprint dataset contains no header.")
    protected boolean ifNoHeader = false;

    @Option(name="-sep", aliases="--separator", usage="Character used as separator in the dataset (by default: \\t for tab).")
    protected String separator = "\t";

    @Option(name="-f", aliases="--filter", usage="Number of the filtered column")
    protected int colFiltered = -1;

    @Option(name="-gal", aliases="--galaxy", usage="For galaxy compliance: formatting pathway output and creating a new one containing log information.")
    protected String galaxy = "";

    /*****MAPPING PARAMETERS*****/

    @Option(name = "-t", aliases="--bioType", usage = "Type of biological object selected : 1 for metabolites or 2 for reactions (by default: metabolites).")
    protected int bioEntityType = 1;

    @Option(name="-name", usage="Number of the file's column containing the metabolite name (by default: 1st column).")
    protected int nameColumn = 1;

    @Option(name="-l", aliases="--layers", usage="List containing the number - separated by comma without blank spaces - of the InChi's layer concerned by the mapping" +
            " (by default: c,h; for a mapping including all the layers, enter c,h,q,p,b,t,i,f,r; for a mapping on formula layer only, enter the -l option with no parameter).")
    protected String inchiLayers = "c,h";

    @Option(name="-inchi", usage="Number of the file's column containing the InChI data.")
    protected int inchiColumn = -1;

    @Option(name="-chebi", usage="Number of the file's column containing the ChEBI data.")
    protected int chebiColumn = -1;

    @Option(name="-idSBML", usage="Number of the file's column containing the metabolite identifier (by default: 2nd column).")
    protected int idSBMLColumn = -1;

    @Option(name="-smiles", usage="Number of the file's column containing the SMILES data.")
    protected int smilesColumn = -1;

    @Option(name="-pubchem", usage="Number of the file's column containing the PubChem identifier.")
    protected int pubchemColumn = -1;

    @Option(name="-inchikey", usage="Number of the file's column containing the InChIKey.")
    protected int inchikeyColumn = -1;

    @Option(name="-kegg", usage="Number of the file's column containing the KEGG identifier.")
    protected int keggColumn = -1;

    @Option(name="-hmdb", usage="Number of the file's column containing the HMDB identifier.")
    protected int hmdbColumn = -1;

    @Option(name="-csid", aliases="--chemspider", usage="Number of the file's column containing the ChemSpider identifier.")
    protected int csidColumn = -1;

    @Option(name="-weight", usage="Number of the file's column containing the weight of the metabolites.")
    protected int weightColumn = -1;


    public void timeCalculation(long elapsedTime){
        long min = elapsedTime / 60000000000L;
        long sec = elapsedTime / 1000000000L - (min * 60L);
        System.out.println("Time to run the process : " + min + "min " + sec + "s");

    }

    public static Boolean testInchiParameter(String[] args){

        for (String arg2 : args) {
            if (Pattern.matches("-inchi[ ]*", arg2)) {
                return true;
            }
        }
        return false;
    }

    @SuppressWarnings("deprecation")
    public static void main(String[] args) {

        long startTime = System.nanoTime();
        Launcher_PathEnr launch = new Launcher_PathEnr();
        CmdLineParser parser = new CmdLineParser(launch);
        String mappingWarnings = "#Warning: By default, a mapping has been set with the name and the SBML id respectively on the 1st and the 2nd column of your dataset.\n" +
                "#Warning: Other mapping available: ChEBI, InChI, InChIKey, SMILES, CSID, PubChem and HMDB (check --help).";

        try {
            parser.parseArgument(args);

            //Print help
            if (launch.phelp) {
                System.out.println("Options:");
                parser.printUsage(System.out);
                exit(0);
            }

            if (launch.version) {
                System.out.println("Version: 1.1");
                exit(0);
            }

            //Error messages for bad parameters
            if(launch.inFileFingerprint==null){
                throw new CmdLineException("-i parameter required");
            }

            if (!Pattern.matches("([chqpbtifr],)*[chqpbtifr]", launch.inchiLayers)) {
                throw new CmdLineException("-l parameter badly formatted");
            }

            if (launch.bioEntityType < 1 || launch.bioEntityType > 6) {
                throw new CmdLineException("Type of biological object must be between 1 and 6.");
            }

            Boolean ifLayerMappingParameter = false, ifInchiMappingParameter = false;
            for (String arg : args) {
                if (Pattern.matches("-l[ ]*", arg)) {
                    ifLayerMappingParameter = true;
                    ifInchiMappingParameter = testInchiParameter(args);
                }
            }
            if(ifLayerMappingParameter && !ifInchiMappingParameter){
                launch.inchiColumn = 2;
                System.out.println("#Warning: InChI layers parameters set without having specified the InChI column (-inchi).\n" +
                        "#Warning: By default, the column used for InChI mapping is the 2nd of your dataset.");
            }else{

                //The user have use any mapping parameters
                int i=0;
                Boolean ifMappingParameter = false;
                for (String arg : args) {
                    if(Pattern.matches("-(name|chebi|inchi|idSBML|smiles|pubchem|inchikey|kegg|hmdb|csid|weight)[ ]*", arg)) {
                        ifMappingParameter = true;
                        break;
                    }
                    i++;
                }
                if(!ifMappingParameter){
                    launch.idSBMLColumn = 2;
                    System.out.println("#Warning: No mapping parameters have been chosen.\n" + mappingWarnings);
                }

                //All mapping parameters are disabled
                i=0;
                if (launch.nameColumn < 1 && launch.chebiColumn < 1 && launch.inchiColumn < 1 && launch.idSBMLColumn  < 1 &&
                        launch.smilesColumn < 1 && launch.pubchemColumn < 1 && launch.inchikeyColumn < 1
                        && launch.keggColumn < 1 && launch.hmdbColumn < 1 && launch.csidColumn < 1
                        && launch.weightColumn < 1) {
                    launch.nameColumn = 1;
                    launch.idSBMLColumn = 2;
                    System.out.println("#Warning: All parameters for mapping your dataset on the SBML are disabled.\n" + mappingWarnings);
                }else {
                    for (String arg : args) {
                        if (Pattern.matches("-name[ ]*", arg) && Pattern.matches("-1[ ]*", args[i + 1])) {
                            System.out.println("#Warning: By disabling the name parameters, name of the entities will not appear.");
                            break;
                        } else {
                            i++;
                        }
                    }
                }
            }

        //Personalised error print for help
        } catch (CmdLineException e) {
            if(e.getMessage().equals("Option \"-l (--layers)\" takes an operand")){
                launch.inchiLayers="";
                Boolean ifInchiMappingParameter = testInchiParameter(args);
                if(!ifInchiMappingParameter) {
                    launch.inchiColumn = 2;
                    System.out.println("#Warning: InChI layers parameters set without having specified the InChI column (-inchi).\n" +
                            "#Warning: By default, the column used for InChI mapping is the 2nd of your dataset.");
                }
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
                (launch.smilesColumn-1), (launch.pubchemColumn-1), (launch.inchikeyColumn-1),
                (launch.keggColumn-1), (launch.hmdbColumn-1), (launch.csidColumn-1), (launch.weightColumn-1)};

        try{
            Fingerprint fingerprint = new Fingerprint(launch.inFileFingerprint,launch.ifNoHeader, launch.separator, (launch.nameColumn-1),
                    mappingColumns, (launch.colFiltered-1));
            Mapping mapping = new Mapping(network, fingerprint.list_entities, inchiLayers,
                    launch.outFileMapping, launch.galaxy, launch.bioEntityType);
            PathwayEnrichment pathEnr = new PathwayEnrichment(network, fingerprint.list_entities, mapping.list_mappedEntities,
                    launch.outFilePathEnr,launch.galaxy, launch.bioEntityType);
        }
        catch (IOException e){
            e.printStackTrace();
        }
        launch.timeCalculation(System.nanoTime() - startTime);
    }
}