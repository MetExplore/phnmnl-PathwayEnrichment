package fr.inra.toulouse.metexplore;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

import static java.lang.System.exit;

public class Fingerprint {

    public int nameColumn, chebiColumn, inchiColumn, idSBMLColumn, smilesColumn, pubchemColum;
    public int inchikeysColumn, keggColumn, hmdColumn, chemspiderColumn, weightColumn,  filteredColumn;
    public String separator, inFileFingerprint;
    public Boolean ifHeader;
    public HashMap<String, String[]> list_metabolites = new HashMap<String, String[]>(); //input file after formatting and filtering

    //TODO: excel parsing

    public Fingerprint (String inFileFingerprint, Boolean ifHeader, String separator, int nameColumn, int[] mappingColumns,
                        int filteredColumn) throws IOException {
        this.inFileFingerprint=inFileFingerprint;
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
        this.hmdColumn=mappingColumns[7];
        this.chemspiderColumn=mappingColumns[8];
        this.weightColumn=mappingColumns[9];
        this.filteredColumn=filteredColumn;

        extractData();
        //TODO: check e.g. inchi lines selected begin with INCHI=
    }

    public void  extractData() throws IOException {

        Boolean isFiltered = (this.filteredColumn >= 0) ? true : false;
        BufferedReader fileBuffer=new BufferedReader(new FileReader(new File(inFileFingerprint)));
        String line;
        int id = 1;
        int[] columnNumbers = {this.nameColumn, this.idSBMLColumn, this.inchiColumn, this.chebiColumn,
                this.smilesColumn, this.pubchemColum, this.inchikeysColumn, this.keggColumn, this.hmdColumn,
                this.chemspiderColumn, this.weightColumn};
        //if (verbose=true) System.out.println(Arrays.toString(columnNumbers));

        if(this.ifHeader) fileBuffer.readLine(); //skip the header

        //Loop on each lines from the input file
        while ((line = fileBuffer.readLine()) != null) {

            String[] lineInFile = line.replaceAll("\"", "").split(this.separator);//splitting by tabulation
            String[] lineFormatted = new String[11];
            //if (verbose=true) System.out.println(Arrays.toString(lineInFile));            

            for (int i = 0; i < columnNumbers.length; i++) {
                putValueIfExists(lineFormatted, lineInFile, i, columnNumbers[i]);
            }
            try {
                if (isFiltered == false || lineInFile[this.filteredColumn] != "") { //optional filtering on a specified column
                    //if (verbose=true) System.out.println(Arrays.toString(lineFormatted));
                    this.list_metabolites.put(""+ id, lineFormatted);//add to hashmap
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

        /*for (String[] lineInFile : list_metabolites.values()) {
            System.out.println(Arrays.toString(lineInFile));
        }*/
    }

    public void putValueIfExists (String[] lineFormatted, String[] lineInFile, int columnInTable, int columnInFile){
        if (columnInFile >= 0) {
            try {
                lineFormatted[columnInTable] = lineInFile[columnInFile];
            } catch (ArrayIndexOutOfBoundsException e) {
                //if a column contains some blank values
                lineFormatted[columnInTable] =  "";
            }
        }else {
            lineFormatted[columnInTable] =  "";
        }
    }

}