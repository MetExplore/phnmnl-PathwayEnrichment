package fr.inra.toulouse.metexplore.launcher;

import fr.inra.toulouse.metexplore.Fingerprint;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.IOException;
import java.util.regex.Pattern;

public class Launcher_Fingerprint extends Launcher{

    /******FILES PARAMETERS*****/

    @Option(name="-i", aliases="--inFile", usage="[REQUIRED] Input file containing a fingerprint (in tsv file format).")
    protected String inFileFingerprint ;

    @Option(name="-o1", usage="Output file name for checking format process (by default: disabled).")
    protected String checkingFile = "";

    /******PARSING PARAMETERS*****/

    @Option(name="--check", usage="Activate this option to check database identifier format.")
    protected boolean noFormatCheck = false;

    @Option(name="--header", usage="Activate this option if the fingerprint dataset contains no header.")
    protected boolean ifNoHeader = false;

    @Option(name="-sep", aliases="--separator", usage="Character used as separator for columns in the dataset (by default: tabulation).")
    protected String columnSeparator = "\t";
    //TODO: --adv without interface
    //TODO: Check separator

    @Option(name="-sepID", aliases="--separatorID", usage="Character used as separator if there are multiple values in a database column of the dataset (by default: ,).")
    protected String IDSeparator = ";";
    //TODO: warning if ',' separator because of the InChI
    //TODO: remove this separator in Fingerprint class for InChI only

    @Option(name="-f", aliases="--filter", usage="Number of the filtered column")
    protected int colFiltered = -1;

    @Option(name="-lWarn", aliases="-layersWarning", usage="When format checking option is activated, points out badly formatted InChI layers")
    protected Boolean layerWarning = false;

    /*****MAPPING PARAMETERS*****/

    @Option(name = "-nameCol", usage = "Number of the file's column containing the bio-entity name (by default: 1st column).")
    protected int nameColumn = 1;

    @Option(name = "-l", aliases = "--layers", usage = "List containing the number - separated by comma without blank spaces - of the InChi's layer concerned by the mapping" +
            " (by default: c,h; for a mapping including all the layers, enter c,h,q,p,b,t,i,f,r; for a mapping on formula layer only, enter the -l option with no parameter).")
    protected String inchiLayers = "c,h";

    protected String[] tab_inchiLayers;

    @Option(name = "-inchi", usage = "Number of the file's column containing the InChI data.")
    protected int inchiColumn = -1;

    @Option(name = "-chebi", usage = "Number of the file's column containing the ChEBI data.")
    protected int chebiColumn = -1;

    @Option(name = "-idSBML", usage = "Number of the file's column containing the SBML identifier (by default: 2nd column).")
    protected int idSBMLColumn = -1;

    @Option(name = "-smiles", usage = "Number of the file's column containing the SMILES data.")
    protected int smilesColumn = -1;

    @Option(name = "-pubchem", usage = "Number of the file's column containing the PubChem identifier.")
    protected int pubchemColumn = -1;

    @Option(name = "-inchikey", usage = "Number of the file's column containing the InChIKey.")
    protected int inchikeyColumn = -1;

    @Option(name = "-kegg", usage = "Number of the file's column containing the KEGG identifier.")
    protected int keggColumn = -1;

    @Option(name = "-hmdb", usage = "Number of the file's column containing the HMDB identifier.")
    protected int hmdbColumn = -1;

    @Option(name = "-csid", aliases = "--chemspider", usage = "Number of the file's column containing the ChemSpider identifier.")
    protected int csidColumn = -1;

    @Option(name = "-weight", usage = "Number of the file's column containing the weight of the metabolites.")
    protected int weightColumn = -1;

    public static Boolean testInchiParameter(String[] args) {

        for (String arg2 : args) {
            if (Pattern.matches("-inchi[ ]*", arg2)) {
                return true;
            }
        }
        return false;
    }

    @SuppressWarnings("deprecation")
    public void printInfo(CmdLineParser parser, String[] args) throws CmdLineException{

        super.printInfo(parser);

        if(!this.noFormatCheck && !this.inchiLayers.equals("c,h")){
            this.layerWarning = true;
        }

        //Error messages for bad parameters
        if(this.inFileFingerprint==null){
            throw new CmdLineException("-i parameter required");
        }
        if (!this.noFormatCheck && !this.inchiLayers.equals("c,h")) {
            this.layerWarning = true;

        }

        if (!Pattern.matches("([chqpbtifr],)*[chqpbtifr]", this.inchiLayers)) {
            throw new CmdLineException("-l parameter badly formatted");
        }

        Boolean ifLayerMappingParameter = false, ifInchiMappingParameter = false;
        for (String arg : args) {
            if (Pattern.matches("-l[ ]*", arg)) {
                ifLayerMappingParameter = true;
                ifInchiMappingParameter = testInchiParameter(args);
            }
        }
        if (ifLayerMappingParameter && !ifInchiMappingParameter) {
            this.inchiColumn = 2;
            write.writeLog("[WARNING] InChI layers parameters set without having specified the InChI column (-inchi).\n" +
                    "[WARNING] By default, the column used for InChI mapping is the 2nd of your dataset.\n");
        } else {

            //The user have use any mapping parameters
            int i = 0;
            Boolean ifMappingParameter = false;
            for (String arg : args) {
                if (Pattern.matches("-(name|chebi|inchi|idSBML|smiles|pubchem|inchikey|kegg|hmdb|csid|weight)", arg)) {
                    ifMappingParameter = true;
                    break;
                }
                i++;
            }
            if (!ifMappingParameter) {
                this.idSBMLColumn = 2;
                write.writeLog("[WARNING] No mapping parameters have been chosen.\n" + mappingWarnings);
            }

            //All mapping parameters are disabled
            i = 0;
            if (this.nameColumn < 1 && this.chebiColumn < 1 && this.inchiColumn < 1 && this.idSBMLColumn < 1 &&
                    this.smilesColumn < 1 && this.pubchemColumn < 1 && this.inchikeyColumn < 1
                    && this.keggColumn < 1 && this.hmdbColumn < 1 && this.csidColumn < 1
                    && this.weightColumn < 1) {
                this.nameColumn = 1;
                this.idSBMLColumn = 2;
                write.writeLog("[WARNING] All parameters for mapping your dataset on the SBML are disabled.\n" + mappingWarnings);
            } else {
                for (String arg : args) {
                    if (Pattern.matches("-nameCol", arg) && Pattern.matches("-[ ]*", args[i + 1])) {
                        throw new CmdLineException("Name column must be positive.\n");
                    } else {
                        i++;
                    }
                }
            }
        }
    }

    public void printError(CmdLineParser parser, CmdLineException e, String[] args) {
        super.printError(parser, e);
        if (e.getMessage().equals("Option \"-l (--layers)\" takes an operand")) {
            this.inchiLayers = "";
            Boolean ifInchiMappingParameter = testInchiParameter(args);
            if (!ifInchiMappingParameter) {
                this.inchiColumn = 2;
                System.out.println("[WARNING] InChI layers parameters set without having specified the InChI column (-inchi).\n" +
                        "[WARNING] By default, the column used for InChI mapping is the 2nd of your dataset.\n");
            }
        } else {
            super.printError(parser, e);
        }
    }

    public Object exec(CmdLineParser parser, String[] args) throws IOException {
        Fingerprint fingerprint = null;

        try {
            parser.parseArgument(args);
            this.printInfo(parser, args);
        } catch (CmdLineException e) {
            this.printError(parser, e, args);
        }

        //Regex for inchiLayers parameter

        this.tab_inchiLayers = this.inchiLayers.replaceAll(" ", "").split(",");
        int[] mappingColumns = {(this.idSBMLColumn - 1), (this.inchiColumn - 1), (this.chebiColumn - 1),
                (this.smilesColumn - 1), (this.pubchemColumn - 1), (this.inchikeyColumn - 1),
                (this.keggColumn - 1), (this.hmdbColumn - 1), (this.csidColumn - 1), (this.weightColumn - 1)};

        try {
            fingerprint = new Fingerprint(this.layerWarning, this.noFormatCheck, this.checkingFile,
                    this.inFileFingerprint, this.ifNoHeader, this.columnSeparator,
                    this.IDSeparator, (this.nameColumn - 1), mappingColumns, tab_inchiLayers,
                    (this.colFiltered - 1));
        } catch (IOException e) {
            e.printStackTrace();
        }

        return fingerprint;
    }

    @SuppressWarnings("deprecation")
    public static void main(String[] args) {
        Launcher_Fingerprint launch = new Launcher_Fingerprint();
        try {
            launch.exec(new CmdLineParser(launch), args);
        }catch (IOException e){}
        timeCalculation(System.nanoTime() - startTime);
    }

}
