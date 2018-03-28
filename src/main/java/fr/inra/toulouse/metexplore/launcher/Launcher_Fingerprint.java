package fr.inra.toulouse.metexplore.launcher;

import fr.inra.toulouse.metexplore.io.Fingerprint;
import fr.inra.toulouse.metexplore.io.WritingComportment;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.IOException;
import java.util.regex.Pattern;

public class Launcher_Fingerprint extends Launcher implements WritingComportment{

    /******FILES PARAMETERS*****/

    @Option(name="-i", aliases="-inFile", usage="[REQUIRED] Input file containing a fingerprint (in tsv file format).")
    protected String inFileFingerprint ;

    @Option(name="-o1", aliases = "--outCheck", usage="Output file name for checking format process (by default: disabled).")
    protected String checkingFile = "";

    /******PARSING PARAMETERS*****/

    @Option(name="-noCheck", usage="Activate this option to check database identifier format.")
    protected boolean noFormatCheck = false;

    @Option(name="-header", usage="Activate this option if the fingerprint dataset contains no header.")
    protected boolean ifNoHeader = false;

    @Option(name="-sep", aliases="-separator", usage="Character used as separator for columns in the dataset (by default: tabulation).")
    protected String columnSeparator = "\t";
    //TODO: --adv without interface
    //TODO: Check separator

    @Option(name="-sepID", aliases="-separatorID", usage="Character used as separator if there are multiple values in a database column of the dataset (by default: ,).")
    protected String IDSeparator = ";";
    //TODO: warning if ',' separator because of the InChI
    //TODO: remove this separator in Fingerprint class for InChI only

    @Option(name="-f", aliases="-filter", usage="Number of the filtered column")
    protected int colFiltered = -1;

    @Option(name="-lWarn", aliases="-layersWarning", usage="When format checking option is activated, points out badly formatted InChI layers")
    protected Boolean layerWarning = false;

    /*****MAPPING PARAMETERS*****/

    @Option(name = "-nameCol", usage = "Number of the file's column containing the bio-entity name (by default: 1st column).")
    protected int nameColumn = 1;

    @Option(name = "-l", aliases = "-layers", usage = "List containing the number - separated by comma without blank spaces - of the InChi's layer concerned by the mapping" +
            " (by default: c,h; for a mapping including all the layers, enter c,h,q,p,b,t,i,f,r; for a mapping on formula layer only, enter the -l option with no parameter).")
    protected String inchiLayers = "c,h";

    protected String[] tab_inchiLayers;

    @Option(name = "-inchi", aliases = "-InChI", usage = "Number of the file's column containing the InChI data.")
    protected int inchiColumn = -1;

    @Option(name = "-chebi",  aliases = {"-ChEBI", "-CHEBI"}, usage = "Number of the file's column containing the ChEBI data.")
    protected int chebiColumn = -1;

    @Option(name = "-idSBML", aliases = {"-id", "-idSbml", "-idsbml"}, usage = "Number of the file's column containing the SBML identifier (by default: 2nd column).")
    protected int idSBMLColumn = -1;

    @Option(name = "-smiles", aliases = "-SMILES", usage = "Number of the file's column containing the SMILES data.")
    protected int smilesColumn = -1;

    @Option(name = "-pubchem", aliases = {"-pub", "-pubChem"}, usage = "Number of the file's column containing the PubChem identifier.")
    protected int pubchemColumn = -1;

    @Option(name = "-inchikey", aliases = {"-inchiKey", "-key"}, usage = "Number of the file's column containing the InChIKey.")
    protected int inchikeyColumn = -1;

    @Option(name = "-kegg", aliases = "-KEGG", usage = "Number of the file's column containing the KEGG identifier.")
    protected int keggColumn = -1;

    @Option(name = "-hmdb", aliases = "-HMDB", usage = "Number of the file's column containing the HMDB identifier.")
    protected int hmdbColumn = -1;

    @Option(name = "-csid", aliases = {"-chemspider", "-chemSpider"}, usage = "Number of the file's column containing the ChemSpider identifier.")
    protected int csidColumn = -1;

    @Option(name = "-mass", aliases = "-weight", usage = "Number of the file's column containing the isotopic mass of the metabolites.")
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

        //required=true in @option required -i with -h and -v !
        if(this.inFileFingerprint==null){
            throw new CmdLineException("-i parameter required");
        }

        if(!this.noFormatCheck && !this.inchiLayers.equals("c,h")){
            this.layerWarning = true;
        }

        //Error messages for bad parameters
        if (!this.noFormatCheck){
            if(!this.inchiLayers.equals("c,h")) {
                this.layerWarning = true;
            }
        }else {
            if(!this.checkingFile.equals("")){
                this.logContent = writeLog(logContent,"[WARNING] Checking format option has been disabled.\n" +
                        "[WARNING] To prevent checking file to be empty, it has been activated by default.\n");
                this.noFormatCheck = false;
                if(!this.inchiLayers.equals("c,h")) {
                    this.layerWarning = true;
                }
            }
            if(this.layerWarning && this.noFormatCheck){
                this.logContent = writeLog(logContent,"[WARNING] Checking format option has been disabled.\n" +
                        "[WARNING] Without checking, layer warnings option will be useless.\n");
            }
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
            this.logContent = writeLog(logContent,"[WARNING] InChI layers parameters set without having specified the InChI column (-inchi).\n" +
                    "[WARNING] By default, the column used for InChI mapping is the 2nd of your dataset.\n");
        } else {

            //The user have use any mapping parameters
            int i = 0;
            Boolean ifMappingParameter = false;
            for (String arg : args) {
                if (Pattern.matches("-(nameCol|name|chebi|inchi|idSBML|smiles|pubchem|inchikey|kegg|hmdb|csid|mass)$", arg)) {
                    ifMappingParameter = true;
                    break;
                }
                i++;
            }
            if (!ifMappingParameter && weightColumn < 1) {
                this.idSBMLColumn = 2;
                this.logContent = writeLog(logContent,"[WARNING] No mapping parameters have been chosen.\n" + mappingWarnings);
            }

            //All mapping parameters are disabled
            ifMappingParameter = false;
            i = 0;
            if (this.nameColumn < 1 && this.chebiColumn < 1 && this.inchiColumn < 1 && this.idSBMLColumn < 1 &&
                    this.smilesColumn < 1 && this.pubchemColumn < 1 && this.inchikeyColumn < 1
                    && this.keggColumn < 1 && this.hmdbColumn < 1 && this.csidColumn < 1
                    && this.weightColumn < 1) {
                //case for name mapping in corresponding Launcher and avoid a -name option in Fingerprint launcher
                for (String arg : args) {
                    if (Pattern.matches("-name$", arg)) {
                        ifMappingParameter = true;
                        if (Pattern.matches("^-.*", args[i + 1])) {
                            this.nameColumn = 1;
                            this.idSBMLColumn = 2;
                            this.logContent = writeLog(logContent, "[WARNING] All parameters for mapping your dataset on the SBML are disabled.\n" + mappingWarnings);
                        }
                    }else {
                        i++;
                    }
                }
                if (!ifMappingParameter){
                    this.nameColumn = 1;
                    this.idSBMLColumn = 2;
                    this.logContent = writeLog(logContent, "[WARNING] All parameters for mapping your dataset on the SBML are disabled.\n" + mappingWarnings);
                }
            } else {
                i = 0;
                for (String arg : args) {
                    //this.nameColumn < 0 : case for name mapping in corresponding Launcher
                    if (this.nameColumn < 1 && Pattern.matches("-nameCol", arg) && Pattern.matches("-.*", args[i + 1])) {
                        throw new CmdLineException("Name column must be positive.\n");
                    } else {
                        i++;
                    }
                }
            }
        }
    }

    public void printError(CmdLineParser parser, CmdLineException e, String[] args) {
        if (e.getMessage().equals("Option \"-l (-layers)\" takes an operand")) {
            this.inchiLayers = "";
            Boolean ifInchiMappingParameter = testInchiParameter(args);
            if (!ifInchiMappingParameter) {
                this.inchiColumn = 2;
                this.logContent = writeLog(logContent,"[WARNING] InChI layers parameters set without having specified the InChI column (-inchi).\n" +
                        "[WARNING] By default, the column used for InChI mapping is the 2nd of your dataset.\n");
            }
        } else {
            super.printError(parser, e);
        }
    }

    public Object analyse(CmdLineParser parser, String[] args) throws IOException {
        Fingerprint fingerprint = null;

        //Regex for inchiLayers parameter
        this.tab_inchiLayers = this.inchiLayers.replaceAll(" ", "").split(",");
        int[] mappingColumns = {(this.idSBMLColumn - 1), (this.inchiColumn - 1), (this.chebiColumn - 1),
                (this.smilesColumn - 1), (this.pubchemColumn - 1), (this.inchikeyColumn - 1),
                (this.keggColumn - 1), (this.hmdbColumn - 1), (this.csidColumn - 1), (this.weightColumn - 1)};

        try {
            fingerprint = new Fingerprint(this.logContent, this.layerWarning, this.noFormatCheck, this.checkingFile,
                    this.inFileFingerprint, this.ifNoHeader, this.columnSeparator,
                    this.IDSeparator, (this.nameColumn - 1), mappingColumns, tab_inchiLayers,
                    (this.colFiltered - 1));
            this.logContent = fingerprint.getLogContent();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return fingerprint;
    }

    public static void exec(Launcher_Fingerprint launch, String[] args) {
        startTime = System.nanoTime();
        CmdLineParser parser = new CmdLineParser(launch);

        try {
            parser.parseArgument(args);
            launch.printInfo(parser, args);
        } catch (CmdLineException e) {
            launch.printError(parser, e, args);
        }

        if (!launch.galaxyFile.equals("")) {
            logFile = launch.createFile(launch.galaxyFile);
            launch.logContent = "";
        }

        try {
            launch.analyse(parser, args);
            if (!launch.galaxyFile.equals("")) launch.writeOutput(launch.logContent, logFile);
        }catch (IOException e){
            e.printStackTrace();
        }
        timeCalculation(System.nanoTime() - startTime);
    }

    @SuppressWarnings("deprecation")
    public static void main(String[] args) {
        exec(new Launcher_Fingerprint(), args);
    }

}
