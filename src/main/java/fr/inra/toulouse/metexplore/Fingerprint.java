package fr.inra.toulouse.metexplore;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.regex.Pattern;

import static java.lang.System.exit;

public class Fingerprint {

    protected int nameColumn, filteredColumn, nbLine = 3;
    protected int[] columnNumbers, nbWarningPerDatabases;
    protected String separator, IDSeparator, warnings ="";
    protected String[] inchiLayers;
    protected Boolean ifNoHeader, noFormatCheck, layerWarning;
    protected ArrayList<String[]> list_entities = new ArrayList<String[]>(); //input file after formatting and filtering
    protected String[] databases = {"InChI","ChEBI","SMILES","PubChem ID","InChIKey","KEGG ID","HMDB ID","ChemSpider ID", "weight"};
    protected BufferedReader fileBuffer;

    //TODO: excel parsing

    public Fingerprint (Boolean layerWarning, Boolean noFormatCheck, String inFileFingerprint, Boolean ifNoHeader, String separator, String IDSeparator, int nameColumn,
                        int[] mappingColumns, String[] inchiLayers, int filteredColumn) throws IOException {
        this.fileBuffer=new BufferedReader(new FileReader(new File(inFileFingerprint)));;
        this.layerWarning=layerWarning;
        this.noFormatCheck=noFormatCheck;
        this.ifNoHeader=ifNoHeader;
        this.separator=separator;
        this.IDSeparator=IDSeparator;
        this.nameColumn=nameColumn;
        this.nbWarningPerDatabases=new int[databases.length];
        for (int i = 0; i< databases.length; i++){
            this.nbWarningPerDatabases[i] = 0;
        }
        this.columnNumbers=new int[mappingColumns.length+1];
        this.columnNumbers[0]=this.nameColumn;
        for (int i = 0; i< mappingColumns.length; i++){
            this.columnNumbers[i+1] = mappingColumns[i];
        }
        this.inchiLayers=inchiLayers;
        this.filteredColumn=filteredColumn;
        extractData();
        //TODO: check e.g. inchi lines selected begin with INCHI=
    }

    public void  extractData() throws IOException {

        Boolean isFiltered = (this.filteredColumn >= 0) ? true : false;
        String line;

        //if (debug) System.out.println(Arrays.toString(columnNumbers));

        if(!this.ifNoHeader){
            fileBuffer.readLine(); //skip the header
            this.nbLine = 2;
        }

        //Loop on each lines from the input file
        while ((line = fileBuffer.readLine()) != null) {

            String[] lineInFile = line.replaceAll("\"", "").split(this.separator);//splitting by tabulation
            String[] lineFormatted = new String[columnNumbers.length];
            //if (debug)
            //System.out.println(Arrays.toString(lineInFile));

            for (int i = 0; i < columnNumbers.length; i++) {
                putValueIfExists(lineFormatted, lineInFile, i, columnNumbers[i]);
            }
            try {
                if (isFiltered == false || lineInFile[this.filteredColumn] != "") { //optional filtering on a specified column
                   // if (debug)
                   //     System.out.println(Arrays.toString(lineFormatted));
                    this.list_entities.add(lineFormatted);//add to hashmap
                }
            } catch (ArrayIndexOutOfBoundsException e) {
                //avoid errors with filtering functionality containing empty values
            }
            nbLine++;
        }
        if (fileBuffer != null) fileBuffer.close();
        if (this.list_entities.size() < 1) {//no extraction = error generation
            System.err.println("File badly formatted");
            exit(1);
        }else if(!noFormatCheck && !testWrongColumn()){
            if (warnings.isEmpty()) System.out.println("All your databases identifiers seems valid.");
            else System.out.println(warnings);
        }

        /*for (String[] lineInFile : list_entities) {
            System.out.println(Arrays.toString(lineInFile));
        }*/
    }

    public Boolean testWrongColumn() {
        Boolean ifWrongCol = false;
        for (int i = 0; i < databases.length; i++) {
            if (nbWarningPerDatabases[i] >= nbLine) {
                System.out.println("[WARNING] For " + databases[i] + " values, all the lines are badly formatted. Check your column number for this parameter.\n");
                ifWrongCol = true;
            }
        }
        return ifWrongCol;
    }

    public void putValueIfExists (String[] lineFormatted, String[] lineInFile, int columnInTable, int columnInFile){
        if (columnInFile >= 0) {
            try {
                String[] tab_ids = lineInFile[columnInFile].split(IDSeparator);
                ArrayList <String> ids = new ArrayList <String>();
                for (String id : tab_ids) {

                    id = id.replace("\\s$", "").replace("^\\s","").replaceAll("\\s{2,}", "");
                    if(id.toUpperCase().equals("NA")) id = "";
                    if (columnInTable > 1){
                        //avoid to replace space for example for pathway name (could be refactored)
                        id = id.replaceAll("\\s", "");
                        if(!this.noFormatCheck && !id.isEmpty()){
                            checkIDFormat(id, lineInFile,columnInTable);
                        }else{
                            nbWarningPerDatabases[columnInTable-2] = nbWarningPerDatabases[columnInTable-2] + 1;
                        }
                    }
                    ids.add(id);
                }
                lineFormatted[columnInTable] = String.join(";", ids);
            } catch (ArrayIndexOutOfBoundsException e) {
                //if a column contains some blank values
                lineFormatted[columnInTable] =  "";
            }
        }else {
            lineFormatted[columnInTable] =  "";
        }
    }

    public void checkIDFormat(String id, String[] lineInFile, int columnInTable){
        String[] patterns = {"^CHEBI:[0-9]+$",".*","^[0-9]*$","^[A-Z]{14}-[A-Z]{10}-[A-Z]$","^[A-Z]{1,2}[0-9]{5}$","^HMDB[0-9]{5}$","^[0-9]*$"};
        String[] warning = {"[WARNING] For " + lineInFile[this.nameColumn] + " (line n°" + this.nbLine + "), ", " is badly formatted: " + id + "\n"};
        String[] databases_subset = new String[databases.length-2];
        System.arraycopy(databases, 1, databases_subset, 0, databases.length-2);

        if (columnInTable == 10){
            try {
                Float.parseFloat(id);
            } catch (NumberFormatException e) {
                setWarnings(warning,"weight",columnInTable);
            }
        } else if (columnInTable > 2){
            if(!Pattern.matches(patterns[columnInTable-3], id)) {
                setWarnings(warning,databases_subset[columnInTable - 3],columnInTable);
            }
        }else if (columnInTable == 2){
            InChI4Galaxy inchi = new InChI4Galaxy(id, this.inchiLayers);
            if(!inchi.validity) {
                String message = "InChI";
                if (this.layerWarning) message+= " (wrong layer: " + inchi.wrongLayer +")";
                setWarnings(warning,message,columnInTable);
            }
        }
    }

    public void setWarnings(String[] warning, String database, int columnInTable) {
        warnings += warning[0] + database + warning[1];
        nbWarningPerDatabases[columnInTable-2] = nbWarningPerDatabases[columnInTable-2] + 1;
    }

}