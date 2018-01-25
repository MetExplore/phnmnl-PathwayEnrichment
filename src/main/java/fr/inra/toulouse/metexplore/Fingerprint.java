package fr.inra.toulouse.metexplore;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import static java.lang.System.exit;

public class Fingerprint {
    BufferedReader fileBuffer;
    public int nameColumn, chebiColumn, inchiColumn, idSBMLColumn, smilesColumn, pubchemColum,
            inchikeysColumn, keggColumn, filteredColumn;
    public String separator;
    public String[] lineInFile, lineFormatted = new String[3];
    public Boolean ifHeader;
    public HashMap<String, String[]> list_metabolites = new HashMap<String, String[]>(); //input file after formating and filtering

    //TODO: excel parsing

    public Fingerprint (String inFileFingerprint, Boolean ifHeader, String separator, int nameColumn, int[] mappingColumns,
                        int filteredColumn) throws IOException {
        this.fileBuffer=new BufferedReader(new FileReader(new File(inFileFingerprint)));
        this.ifHeader=ifHeader;
        this.separator=separator;
        this.nameColumn=nameColumn;
        this.idSBMLColumn=mappingColumns[0];
        this.inchiColumn=mappingColumns[1];
        this.chebiColumn=mappingColumns[2];
        this.smilesColumn=mappingColumns[3];
        this.pubchemColum=mappingColumns[4];
        this.inchikeysColumn=mappingColumns[5];
        this.keggColumn=mappingColumns[6];
        this.filteredColumn=filteredColumn;
    }

    public void  extractData() throws IOException {

        Boolean isFiltered = (filteredColumn >= 0) ? true : false;

        String line;
        int id = 1;

        if(ifHeader) this.fileBuffer.readLine(); //skip the header

        //Loop on each lines from the input file
        while ((line = this.fileBuffer.readLine()) != null) {
            this.lineInFile = line.replaceAll("\"", "").split(this.separator);//splitting by tabulation

            int[] columnNumbers = {this.nameColumn, this.idSBMLColumn, this.inchiColumn, this.chebiColumn,
            this.smilesColumn, this.pubchemColum, this.inchikeysColumn, this.keggColumn};
            for (int i = 0; i < columnNumbers.length; i++) {
                putValueIfExists(i, columnNumbers[i]);
            }

            try {
                if (isFiltered == false || this.lineInFile[this.filteredColumn] != "") { //optional filtering on a specified column
                    this.list_metabolites.put("" + id, this.lineFormatted);//add to hashmap
                    id++;
                }
            } catch (ArrayIndexOutOfBoundsException e) {//avoid errors with filtering functionality containing empty values
            }
        }
        if (fileBuffer != null) fileBuffer.close();
        if (this.list_metabolites.size() < 1) {//no extraction = error generation
            System.err.println("File badly formatted");
            exit(1);
        }
    }

    public void putValueIfExists (int columnInTable, int columnInFile){
        this.lineFormatted[columnInTable] = (columnInFile >= 0) ? this.lineInFile[columnInFile] : "";
    }

}