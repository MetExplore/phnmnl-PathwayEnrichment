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

    @Option(name="-o1", usage="Output file name for checking format process (by default: disabled).")
    protected String checkingFile = "";

    @Option(name="-o2", aliases="--outMap", usage="Output file name for mapping result (by default: mapping.tsv).")
    protected String outFileMapping = "mapping.tsv";

    @Option(name="-o3", aliases="--outPath", usage="Output file name for pathway enrichment result (by default: pathwayEnrichment.tsv).")
    protected String outFilePathEnr = "pathwayEnrichment.tsv";

    @Option(name="-s", aliases="--sbml", usage="SBML file name (by default: Recon v2.02).")
    protected String sbml = "data/recon2.02_without_compartment.xml";

    /******PARSING PARAMETERS*****/

    @Option(name="--check", usage="Activate this option to check database identifier format.")
    protected boolean noFormatCheck = false;

    @Option(name="--header", usage="Activate this option if the fingerprint dataset contains no header.")
    protected boolean ifNoHeader = false;

    @Option(name="-sep", aliases="--separator", usage="Character used as separator for columns in the dataset (by default: tabulation).")
    protected String columnSeparator = "\t";

    @Option(name="-sepID", aliases="--separatorID", usage="Character used as separator if there are multiple values in a database column of the dataset (by default: ,).")
    protected String IDSeparator = ";";
    //TODO: warning if ',' separator because of the InChI
    //TODO: remove this separator in Fingerprint class for InChI only

    @Option(name="-f", aliases="--filter", usage="Number of the filtered column")
    protected int colFiltered = -1;

    @Option(name="-gal", aliases="--galaxy", usage="For galaxy compliance: formatting pathway output and creating a new one containing log information.")
    protected String galaxy = "";

    @Option(name="-lWarn", aliases="-layersWarning", usage="When format checking option is activated, points out badly formatted InChI layers")
    protected Boolean layerWarning = false;

    /*****MAPPING PARAMETERS*****/

    @Option(name = "-t", aliases="--type", usage = "1 for metabolites, 2 for reactions, 3 for pathway, 4 for enzyme, 5 for protein, 6 for gene (by default: metabolites).")
    protected int entityType2Map = 1;

    @Option(name = "-tEnr", aliases="--typeEnr", usage = "1 for metabolites, 2 for reactions, 3 for pathway, 4 for enzyme, 5 for protein, 6 for gene (by default: pathways).")
    protected int entityType2Enrich = 3;

    @Option(name="-name", usage="Activate this option for a name mapping .")
    protected Boolean nameMapping = false;

    @Option(name="-nameCol", usage="Number of the file's column containing the bio-entity name (by default: 1st column).")
    protected int nameColumn = 1;

    @Option(name="-l", aliases="--layers", usage="List containing the number - separated by comma without blank spaces - of the InChi's layer concerned by the mapping" +
            " (by default: c,h; for a mapping including all the layers, enter c,h,q,p,b,t,i,f,r; for a mapping on formula layer only, enter the -l option with no parameter).")
    protected String inchiLayers = "c,h";

    @Option(name="-inchi", usage="Number of the file's column containing the InChI data.")
    protected int inchiColumn = -1;

    @Option(name="-chebi", usage="Number of the file's column containing the ChEBI data.")
    protected int chebiColumn = -1;

    @Option(name="-idSBML", usage="Number of the file's column containing the SBML identifier (by default: 2nd column).")
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

    //TODO: --adv without interface
    //TODO: Check separator

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
        String mappingWarnings = "[WARNING] By default, a mapping has been set with the name and the SBML id respectively on the 1st and the 2nd column of your dataset.\n" +
                "[WARNING] Other mapping available: ChEBI, InChI, InChIKey, SMILES, CSID, PubChem and HMDB (check --help).\n";

        try {
            parser.parseArgument(args);
            WritingComportment write = new WritingComportment(launch.galaxy);

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

            if (launch.entityType2Map < 1 || launch.entityType2Map > 6 ||
                    launch.entityType2Enrich < 1 || launch.entityType2Enrich > 6) {
                throw new CmdLineException("Type of biological object must be between 1 and 6.");
            }

           if(!launch.noFormatCheck && !launch.inchiLayers.equals("c,h")){
                launch.layerWarning = true;
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
                write.writeLog("[WARNING] InChI layers parameters set without having specified the InChI column (-inchi).\n" +
                        "[WARNING] By default, the column used for InChI mapping is the 2nd of your dataset.\n");
            }else{

                //The user have use any mapping parameters
                int i=0;
                Boolean ifMappingParameter = false;
                for (String arg : args) {
                    if(Pattern.matches("-(name|chebi|inchi|idSBML|smiles|pubchem|inchikey|kegg|hmdb|csid|weight)", arg)) {
                        ifMappingParameter = true;
                        break;
                    }
                    i++;
                }
                if(!ifMappingParameter){
                    launch.idSBMLColumn = 2;
                    write.writeLog("[WARNING] No mapping parameters have been chosen.\n" + mappingWarnings);
                }

                //All mapping parameters are disabled
                i=0;
                if (launch.nameColumn < 1 && launch.chebiColumn < 1 && launch.inchiColumn < 1 && launch.idSBMLColumn  < 1 &&
                        launch.smilesColumn < 1 && launch.pubchemColumn < 1 && launch.inchikeyColumn < 1
                        && launch.keggColumn < 1 && launch.hmdbColumn < 1 && launch.csidColumn < 1
                        && launch.weightColumn < 1) {
                    launch.nameColumn = 1;
                    launch.idSBMLColumn = 2;
                    write.writeLog("[WARNING] All parameters for mapping your dataset on the SBML are disabled.\n" + mappingWarnings);
                }else {
                    for (String arg : args) {
                        if (Pattern.matches("-nameCol", arg) && Pattern.matches("-[ ]*", args[i + 1])) {
                            throw new CmdLineException("Name column must be positive.\n");
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
                    WritingComportment write = new WritingComportment(launch.galaxy) ;
                    System.out.println("[WARNING] InChI layers parameters set without having specified the InChI column (-inchi).\n" +
                            "[WARNING] By default, the column used for InChI mapping is the 2nd of your dataset.\n");
                }
            }else {
                System.err.println(e.getMessage());
                System.err.println("Options:");
                parser.printUsage(System.err);
                exit(1);
            }
        }

        WritingComportment write = new WritingComportment(launch.galaxy) ;
        //Regex for inchiLayers parameter
        String[] inchiLayers = launch.inchiLayers.replaceAll(" ","").split(",");
        //Extract SBML
        //BioNetwork network = (new JSBML2Bionetwork4Galaxy(launch.sbml)).getBioNetwork();
        int[] mappingColumns = {(launch.idSBMLColumn-1), (launch.inchiColumn-1), (launch.chebiColumn-1),
                (launch.smilesColumn-1), (launch.pubchemColumn-1), (launch.inchikeyColumn-1),
                (launch.keggColumn-1), (launch.hmdbColumn-1), (launch.csidColumn-1), (launch.weightColumn-1)};

        try{
            Fingerprint fingerprint = new Fingerprint(launch.layerWarning,launch.noFormatCheck,launch.checkingFile,launch.inFileFingerprint,launch.ifNoHeader, launch.columnSeparator,
                    launch.IDSeparator,(launch.nameColumn-1),mappingColumns, inchiLayers,(launch.colFiltered-1));
            /*Mapping mapping = new Mapping(network, fingerprint.list_entities, inchiLayers, launch.nameMapping,
                    launch.outFileMapping, launch.galaxy, launch.entityType2Map);
            PathwayEnrichment pathEnr = new PathwayEnrichment(network, fingerprint.list_entities, mapping.list_mappedEntities,
                    launch.outFilePathEnr,launch.galaxy, launch.entityType2Map, launch.entityType2Enrich);
            write.writeOutputInfo();*/
        }
        catch (IOException e){
            e.printStackTrace();
        }
        launch.timeCalculation(System.nanoTime() - startTime);
    }
}