/*******************************************************************************
 * Copyright INRA
 *
 *  Contact: ludovic.cottret@toulouse.inra.fr
 *
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *  In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *  The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *******************************************************************************/

package fr.inra.toulouse.metexplore.launcher;

import fr.inra.toulouse.metexplore.io.Fingerprint;
import fr.inra.toulouse.metexplore.io.WritingComportment;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.IOException;
import java.util.regex.Pattern;

public class Launcher_Fingerprint extends Launcher implements WritingComportment {

    /******FILES PARAMETERS*****/

    @Option(name = "-i", aliases = "-inFile", usage = "[REQUIRED] Input file containing a fingerprint (in tsv file format).")
    protected String inFileFingerprint;

    @Option(name = "-o1", aliases = "--outCheck", usage = "Output file name for checking format process (by default: disabled).")
    protected String checkingFile = "checking_format.tsv";

    /******PARSING PARAMETERS*****/

    @Option(name = "-noCheck", usage = "Activate this option to check database identifier format.")
    protected boolean noFormatCheck = false;

    @Option(name = "-header", usage = "Activate this option if the fingerprint dataset contains no header.")
    protected boolean ifNoHeader = false;

    @Option(name = "-sep", aliases = "-separator", usage = "Character used as separator for columns in the dataset (by default: tabulation).")
    protected String columnSeparator = "\t";
    //TODO: --adv without interface
    //TODO: Check separator

    @Option(name = "-sepID", aliases = "-separatorID", usage = "Character used as separator if there are multiple values in a database column of the dataset (by default: ,).")
    protected String IDSeparator = ";";
    //TODO: warning if ',' separator because of the InChI
    //TODO: remove this separator in Fingerprint class for InChI only

    @Option(name = "-f", aliases = "-filter", usage = "Number of the filtered column")
    protected int colFiltered = -1;

    @Option(name = "-lWarn", aliases = "-layersWarning", usage = "When format checking option is activated, points out badly formatted InChI layers")
    protected Boolean layerWarning = false;

    /*****MAPPING PARAMETERS*****/

    @Option(name = "-nameCol", usage = "Number of the file's column containing the bio-entity name (by default: 1st column).")
    protected int nameColumn = -1;

    @Option(name = "-l", aliases = "-layers", usage = "List containing the number - separated by comma without blank spaces - of the InChi's layer concerned by the mapping" +
            " (by default: c,h; for a mapping including all the layers, enter c,h,q,p,b,t,i,f,r; for a mapping on formula layer only, enter the -l option with no parameter).")
    protected String inchiLayers = "c,h";

    protected String[] tab_inchiLayers;

    @Option(name = "-inchi", aliases = "-InChI", usage = "Number of the file's column containing the InChI data.")
    protected int inchiColumn = -1;

    @Option(name = "-chebi", aliases = {"-ChEBI", "-CHEBI"}, usage = "Number of the file's column containing the ChEBI data.")
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
    public void printInfo(CmdLineParser parser, String[] args) throws CmdLineException {

        super.printInfo(parser);

        //required=true in @option required -i with -h and -v !
        if (this.inFileFingerprint == null) {
            throw new CmdLineException("-i parameter required");
        }

        //Check inchi Layers format
        if (!Pattern.matches("([chqpbtifr],)*[chqpbtifr]", this.inchiLayers)) {
            throw new CmdLineException("-l parameter badly formatted: it must be a list containing the number - separated by comma without blank spaces - of the InChi's layer concerned by the mapping " +
                    "(by default: c,h; for a mapping including all the layers, enter c,h,q,p,b,t,i,f,r; for a mapping on formula layer only, enter the -l option with no parameter)");
        }


        //Case for layers settings without InChI
        Boolean ifLayerMappingParameter = false, ifInchiMappingParameter = false;
        for (String arg : args) {
            if (Pattern.matches("-l[ ]*", arg)) {
                ifLayerMappingParameter = true;
                ifInchiMappingParameter = testInchiParameter(args);
            }
        }
        if (ifLayerMappingParameter && !ifInchiMappingParameter) {
            this.inchiColumn = 2;
            writeLog("[WARNING] InChI layers parameters set without having specified the InChI column (-inchi).\n" +
                    "[WARNING] By default, the column used for InChI mapping is the 2nd of your dataset.\n");
        }


        //Mapping parameters
        int i = 0;
        Boolean ifMappingParameter = false;
        String mappingWarnings = " By default, it was set on the SBML identifiers at the 2nd column of your dataset." +
                " Other mapping available: ChEBI, InChI, InChIKey, SMILES, CSID, PubChem, isotopic mass and HMDB (check -help).\n";

        for (String arg : args) {
            if (Pattern.matches("-(name|chebi|inchi|idsbml|id|smiles|pubchem|inchikey|key|kegg|hmdb|csid|mass|weight)$", arg.toLowerCase())) {
                ifMappingParameter = true;
                if (Pattern.matches("(^-.*|0)", args[i + 1]) && !arg.equals("-name")) {
                    writeLog("[WARNING] " + arg + " column parameter must be positive.\n");
                }
            }
            i++;
        }
        //also catch "-prec" only default setting with weight column
        if (this.weightColumn < 1 && !ifMappingParameter) {
            this.idSBMLColumn = 2;
            writeLog("[WARNING] No mapping parameters has been chosen." + mappingWarnings);
        }

        //All mapping parameters are disabled (case with all set to negative values)
        ifMappingParameter = false;
        i = 0;
        String mappingNegativeWarnings = "[WARNING] All your mapping parameters have negative column." + mappingWarnings;
        if (this.chebiColumn < 1 && this.inchiColumn < 1 && this.idSBMLColumn < 1 &&
                this.smilesColumn < 1 && this.pubchemColumn < 1 && this.inchikeyColumn < 1
                && this.keggColumn < 1 && this.hmdbColumn < 1 && this.csidColumn < 1
                && this.weightColumn < 1) {
            //case for name mapping in corresponding Launcher and avoid a -name option in Fingerprint launcher
            for (String arg : args) {
                if (Pattern.matches("-name$", arg)) {
                    ifMappingParameter = true;
                    if (Pattern.matches("(^-.*|0)", args[i + 1])) {
                        this.idSBMLColumn = 2;
                        writeLog(mappingNegativeWarnings);
                    }
                } else {
                    i++;
                }
            }
            if (!ifMappingParameter) {
                this.idSBMLColumn = 2;
                writeLog(mappingNegativeWarnings);
            }
        }

        testParameters(args);

    }

    public void testParameters(String[] args) throws CmdLineException2{

        //check name column setting
        int i = 0;
        Boolean ifNameColumn = false;
        // check if name parameters have been called with negative values
        for (String arg : args) {
            if (Pattern.matches("-name.*", arg)) {
                ifNameColumn = true;
                /*if(Pattern.matches("-.*", args[i + 1])) {
                    writeLog("[WARNING] "+  arg + " column parameter must be positive.\n");
                }*/
            } else {
                i++;
            }
        }
        String nameWarning = "; by default it was set to the 1rst column.\n";
        if (!ifNameColumn) {
            //no name parameters have been called
            this.nameColumn = 1;
            writeLog("[WARNING] No column number has been chosen for the name of the chemicals" + nameWarning);
        } else {
            //name parameters have been called but all with negative column
            i = 0;
            ifNameColumn = false;
            String nameNegativeWarning = "[WARNING] Your column number for the name of the chemicals is negative" + nameWarning;
            for (String arg : args) {
                if (Pattern.matches("-name$", arg)) {
                    ifNameColumn = true;
                    if (Pattern.matches("^-\\d", args[i + 1]) && this.nameColumn < 1) {
                        this.nameColumn = 1;
                        writeLog(nameNegativeWarning);
                    }
                } else {
                    i++;
                }
            }
            //no nameMapping parameters have been called but only nameColum with negative values
            if (!ifNameColumn && this.nameColumn < 1) {
                this.nameColumn = 1;
                writeLog(nameNegativeWarning);
            }
        }

        if (!this.noFormatCheck) {
            //The user has set inchi layers and has not disabled the checking step
            if (!this.inchiLayers.equals("c,h")) {
                this.layerWarning = true;
            }
        } else if (this.layerWarning && this.noFormatCheck) {
            writeLog("[WARNING] Checking format option has been disabled.\n" +
                    "[WARNING] Without checking, layer warnings option will be useless.\n");
        }
    }

    class CmdLineException2 extends CmdLineException{

        public CmdLineException2(String message){
            super(message);
            System.out.println("Vous essayez d'instancier une classe Ville avec un nombre d'habitants négatif !");

        }

    }

    public void printError(CmdLineParser parser, CmdLineException e, String[] args) {
        if (!e.getMessage().equals("Option \"-l (-layers)\" takes an operand")) {
            super.printError(parser, e);
        } else {
            this.inchiLayers = "";
            Boolean ifInchiMappingParameter = testInchiParameter(args);
            if (!ifInchiMappingParameter) {
                this.inchiColumn = 2;
                writeLog("[WARNING] InChI layers parameters set without having specified the InChI column (-inchi).\n" +
                        "[WARNING] By default, the column used for InChI mapping is the 2nd of your dataset.\n");
            }
            try {
                testParameters(args);
            }catch (CmdLineException2 e2){
                    super.printError(parser, e2);
            }
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
            fingerprint = new Fingerprint(this.logContent, this.galaxyFile, this.layerWarning, this.noFormatCheck,
                    this.inFileFingerprint, this.ifNoHeader, this.columnSeparator,
                    this.IDSeparator, (this.nameColumn - 1), mappingColumns, tab_inchiLayers,
                    (this.colFiltered - 1), this.checkingFile);
            this.logContent = fingerprint.getLogContent();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return fingerprint;
    }

    public void setFile(){
        if (!this.galaxyFile.equals("")) {
            this.logContent = "";
            logFile = this.createFile(this.galaxyFile);
        }
    }

    public static void exec(Launcher_Fingerprint launch, String[] args) {
        startTime = System.nanoTime();
        CmdLineParser parser = new CmdLineParser(launch);

        try {
            parser.parseArgument(args);
            launch.setFile();
            launch.printInfo(parser, args);
        } catch (CmdLineException e) {
            //setFile here because of -l parameter catch
            launch.setFile();
            launch.printError(parser, e, args);
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
