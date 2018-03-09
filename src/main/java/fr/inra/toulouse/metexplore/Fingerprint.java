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

    protected int nameColumn, chebiColumn, inchiColumn, idSBMLColumn, smilesColumn, pubchemColum;
    protected int inchikeysColumn, keggColumn, hmdColumn, chemspiderColumn, weightColumn,  filteredColumn;
    protected int nbLine = 3;
    protected String separator, IDSeparator, inFileFingerprint;
    protected String[] inchiLayers;
    protected Boolean ifNoHeader;
    protected ArrayList<String[]> list_entities = new ArrayList<String[]>(); //input file after formatting and filtering

    //TODO: excel parsing

    public Fingerprint (String inFileFingerprint, Boolean ifNoHeader, String separator, String IDSeparator, int nameColumn,
                        int[] mappingColumns, String[] inchiLayers, int filteredColumn) throws IOException {
        this.inFileFingerprint=inFileFingerprint;
        this.ifNoHeader=ifNoHeader;
        this.separator=separator;
        this.IDSeparator=IDSeparator;
        this.nameColumn=nameColumn;
        this.idSBMLColumn=mappingColumns[0];
        this.inchiColumn=mappingColumns[1];
        this.chebiColumn=mappingColumns[2];
        this.smilesColumn=mappingColumns[3];
        this.pubchemColum=mappingColumns[4];
        this.inchikeysColumn=mappingColumns[5];
        this.keggColumn=mappingColumns[6];
        this.hmdColumn=mappingColumns[7];
        this.chemspiderColumn=mappingColumns[8];
        this.weightColumn=mappingColumns[9];
        this.inchiLayers=inchiLayers;
        this.filteredColumn=filteredColumn;

        extractData();
        //TODO: check e.g. inchi lines selected begin with INCHI=
    }

    public void  extractData() throws IOException {

        Boolean isFiltered = (this.filteredColumn >= 0) ? true : false;
        BufferedReader fileBuffer=new BufferedReader(new FileReader(new File(inFileFingerprint)));
        String line;
        int[] columnNumbers = {this.nameColumn, this.idSBMLColumn, this.inchiColumn, this.chebiColumn,
                this.smilesColumn, this.pubchemColum, this.inchikeysColumn, this.keggColumn, this.hmdColumn,
                this.chemspiderColumn, this.weightColumn};
        //if (debug) System.out.println(Arrays.toString(columnNumbers));

        if(!this.ifNoHeader){
            fileBuffer.readLine(); //skip the header
            nbLine = 2;
        }

        //Loop on each lines from the input file
        while ((line = fileBuffer.readLine()) != null) {

            String[] lineInFile = line.replaceAll("\"", "").split(this.separator);//splitting by tabulation
            String[] lineFormatted = new String[11];
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
        }else {
            //TODO: good format
        }

        /*for (String[] lineInFile : list_entities) {
            System.out.println(Arrays.toString(lineInFile));
        }*/
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
                        if(!id.isEmpty()) checkIDFormat(id, lineInFile,columnInTable);
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
        String[] databases = {"ChEBI","SMILES","PubChem ID","InChIKey","KEGG ID","HMDB ID","ChemSpider ID"};
        String[] warning = {"[WARNING] For " + lineInFile[this.nameColumn] + " (line n°" + this.nbLine + "), ", " is badly formatted: " + id};

        if (columnInTable == 10){
            try {
                Float.parseFloat(id);
            } catch (NumberFormatException e) {
                System.out.println(warning[0] + "weight" + warning[1]);
            }
        } else if (columnInTable > 2){
            if(!Pattern.matches(patterns[columnInTable-3], id)) {
                System.out.println(warning[0] + databases[columnInTable - 3] + warning[1]);
            }
        }
        else if (columnInTable == 2){
            InChI4Galaxy inchi = new InChI4Galaxy(id, inchiLayers);
            /*System.out.println(inchi.validity);
            System.out.println(inchi.connectivity);
            System.out.println(inchi.hLayer);
            System.out.println(inchi.protonationLayer);
            System.out.println(inchi.dbStereoLayer);
            System.out.println(inchi.tetraStereoLayer);
            System.out.println(inchi.isotopicLayer);
            System.out.println(inchi.fixedLayer);
            System.out.println(inchi.reconnectedLayer);*/
            if(!inchi.validity) {
                System.out.println(warning[0] + "InChI" + warning[1]);
            }
        }
    }

}