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

package fr.inra.toulouse.metexplore.io;

import fr.inra.toulouse.metexplore.biodata.InChI;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Pattern;

import static java.lang.System.exit;

public class Fingerprint implements WritingComportment {

    protected int nameColumn, filteredColumn, nbLine = 2, nbWarningPerLine, nbSelectedDatabases;
    protected int[] columnNumbers, nbWarningPerDatabases;
    protected String galaxy, separator, IDSeparator, warnings ="", logContent, checkingFile;
    protected String[] inchiLayers;
    protected Boolean ifNoHeader, noFormatCheck, layerWarning;
    protected ArrayList<String[]> list_entities = new ArrayList<String[]>(); //input file after formatting and filtering
    protected String[] databases = {"InChI","ChEBI","SMILES","PubChem ID","InChIKey","KEGG ID","HMDB ID","ChemSpider ID", "mass"};
    protected BufferedReader inBuffer;
    protected BufferedWriter outBuffer;

    public ArrayList<String[]> getEntityList() {
        return list_entities;
    }

    public String getLogContent() {
        return logContent;
    }

    //TODO: parsing for Excel files

    public Fingerprint (String logContent, String galaxy, Boolean layerWarning, Boolean noFormatCheck, String inFileFingerprint, Boolean ifNoHeader, String separator, String IDSeparator, int nameColumn,
                        int[] mappingColumns, String[] inchiLayers, int filteredColumn, String checkingFile) throws IOException {

        this.logContent = logContent;
        this.galaxy = galaxy;
        this.checkingFile = checkingFile;
        this.inBuffer=new BufferedReader(new FileReader(new File(inFileFingerprint)));
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

    }
    
    public void setFileChecking(){
        File f = new File(this.checkingFile);
        if (f.isFile()) {
            //write a new file if already exists
            f.delete();
        }
        try {
            f.createNewFile();

            this.outBuffer = new BufferedWriter(new FileWriter(f));
            if (layerWarning) {
                this.outBuffer.write("Line number\tBioentity name\tDatabase\tInChI layer\tWrong value\n");
            } else {
                this.outBuffer.write("Line number\tBioentity name\tDatabase\tWrong value\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void  extractData() throws IOException {

        Boolean isFiltered = this.filteredColumn >= 0;
        String line;

        //System.out.println(Arrays.toString(columnNumbers));
        //get the number of positive values (i.e., the number of used database)
        int[] datasesColumnNumbers = Arrays.copyOfRange(columnNumbers, 2, columnNumbers.length);
        nbSelectedDatabases = datasesColumnNumbers.length - Arrays.toString(datasesColumnNumbers).replaceAll("[^-]+", "").length();
        if(!this.ifNoHeader){
            inBuffer.readLine(); //skip the header
            this.nbLine--;
        }

        //Loop on each lines from the input file
        while ((line = inBuffer.readLine()) != null) {

            String[] lineInFile = line.replaceAll("\"", "").split(this.separator);//splitting by tabulation
            String[] lineFormatted = new String[columnNumbers.length];
            //if (debug)
            //System.out.println(Arrays.toString(lineInFile));

            for (int i = 0; i < columnNumbers.length; i++) {
                parseIDValueIfExists(lineFormatted, lineInFile, i, columnNumbers[i]);
            }
            try {
                if (!isFiltered || !lineInFile[this.filteredColumn].equals("")) { //optional filtering on a specified column
                   // if (debug)
                   //     System.out.println(Arrays.toString(lineFormatted));
                    this.list_entities.add(lineFormatted);//add to hashmap
                }
            } catch (ArrayIndexOutOfBoundsException e) {
                //avoid errors with filtering functionality containing empty values
            }
            //System.out.println(nbLine);
            nbLine++;
        }
        if (inBuffer != null) inBuffer.close();
        if (this.list_entities.size() < 1) {//no extraction = error generation
            String message="[FATAL] File badly formatted. Please, check the column number of your -f parameter if enable. Otherwise, check your column separator (default: tabulated file, set with -sep parameter).";
            sysExit(this.logContent,message,this.galaxy,10);
        }else if(!noFormatCheck && !ifWrongDatabase()){
            if (warnings.isEmpty()) this.logContent = writeLog(this.logContent,"All your databases identifiers seem valid.\n");
            else{
                setFileChecking();
                this.outBuffer.write(warnings);
                outBuffer.close();
                this.logContent = writeLog(this.logContent,"[WARNING] Some database identifiers are badly formatted, please take a look to \"checking_format.tsv\".\n");
            }
        }

        /*for (String[] lineInFile : list_entities) {
            System.out.println(Arrays.toString(lineInFile));
        }*/

    }

    public Boolean ifWrongDatabase() {
        //tests if all the selected database have all there values wrong
        //else, if they are errors, write into a file

        int nbWrongDatabase = 0;
        String wrongDatabaseWarning = "";

        //check the number of database with all the values which are badly formatted
        for (int i = 0; i < databases.length; i++) {
            if (nbWarningPerDatabases[i] == nbLine-1) {
                wrongDatabaseWarning +="[WARNING] For " + databases[i] + " values, all the lines are badly formatted. Check your column number for this parameter.\n";
                nbWrongDatabase++;
            }
        }

        //if all the selected database have all there values wrong, exit
        // (if they are > 0, meaning they are not included only name or sbml id mapping)
        if (nbSelectedDatabases > 0 && nbWrongDatabase == nbSelectedDatabases) {
            String message="[FATAL] All the values of the selected database(s) are badly formatted. Please check the column number set for these databases.\n" +
                    "-noCheck to ignore the bad format exit and run the analysis anyway.";
            sysExit(this.logContent,message,this.galaxy,11);
        }
        //else write errors
        if(nbWrongDatabase>0) this.logContent = writeLog(this.logContent, wrongDatabaseWarning);

        return (nbWrongDatabase>0);
    }

    public void parseIDValueIfExists (String[] lineFormatted, String[] lineInFile, int columnInTable, int columnInFile){
    //formatting id values if they are blank
    //separate them if they are many among one column
    //check their format if they are not corresponding to name or sbml id

        if (columnInFile >= 0) {
            try {
                //parse values separated by a designed character from a line of the fingerprint dataset
                String[] tab_ids = lineInFile[columnInFile].split(IDSeparator);
                ArrayList <String> ids = new ArrayList<>();
                nbWarningPerLine = 0;
                for (String id : tab_ids) {
                    //catch blank characters
                    id = id.replace("\\s$", "").replace("^\\s","").replaceAll("\\s{2,}", "");
                    if(id.toUpperCase().equals("NA")) id = "";
                    //discard name and sbml id column from format checking process
                    if (columnInTable > 1){
                        //avoid to replace space in name and id parsing step (e.g. for pathway name)
                        id = id.replaceAll("\\s", "");
                        if(!this.noFormatCheck && !id.isEmpty()){
                            checkIDFormat(id, lineInFile,columnInTable);
                        }else{
                            //blank line are included in line count
                            nbWarningPerDatabases[columnInTable-2] = nbWarningPerDatabases[columnInTable-2] + 1;
                        }
                    }
                    ids.add(id);
                }
                //join multiple identifiers from the same database column
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

        //discard InChI and mass
        String[] databases_subset = new String[databases.length-2];
        System.arraycopy(databases, 1, databases_subset, 0, databases.length-2);

        String database ="", layer="";

        if (columnInTable == 10){
            //check mass is a float
            try {
                Float.parseFloat(id);
            } catch (NumberFormatException e) {
                database = "weight"; //set error type
            }
        } else if (columnInTable > 2){
            //check using above regex
            if(!Pattern.matches(patterns[columnInTable-3], id)) {
                database = databases_subset[columnInTable - 3]; //set error type
            }
        }else if (columnInTable == 2){
            //check using InChI object attributes
            InChI inchi = new InChI(id, this.inchiLayers);
            if(!inchi.validity) { //if error in format
                database = "InChI";
                //list of wrong inchi layers
                if (this.layerWarning) layer+= String.join("; ", inchi.wrongLayer);
            }
        }

        if(!database.equals("")) { //if error in format
            setWarnings(id, lineInFile,database,layer);

            //if multiple ID values for one database column in the same line, count only one wrong line
            if(nbWarningPerLine == 0) {
                nbWarningPerDatabases[columnInTable - 2] = nbWarningPerDatabases[columnInTable - 2] + 1;
            }

            //set nb error for case of multiple ID values for one database column in the same line
            nbWarningPerLine++;
        }
    }

    public void setWarnings(String id, String[] lineInFile, String database, String layer) {

        if (layerWarning) {
            this.warnings += this.nbLine + "\t" + lineInFile[this.nameColumn] + "\t" + database + "\t" + layer + "\t" + id + "\n";
        } else {
            this.warnings += this.nbLine + "\t" + lineInFile[this.nameColumn] + "\t" + database + "\t" + id + "\n";
        }

    }

    public ArrayList<String[]> getList_entities() {
        return list_entities;
    }
}