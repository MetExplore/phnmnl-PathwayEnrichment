package fr.inra.toulouse.metexplore;

import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPathway;
import parsebionet.biodata.BioPhysicalEntity;
import parsebionet.biodata.BioRef;
import parsebionet.statistics.PathwayEnrichment;

import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.File;

import java.util.Set;
import java.util.HashSet;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.Collections;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Comparator;

import static java.lang.System.exit;

public class MetExplore4Galaxy {

    public String inputFile, outputFileMapping, outputFilePathEnr, outputFileInfo;
    public int nameColumn, chebiColumn, inchiColumn, idSBMLColumn, filteredColumn;
    public String[] inchiLayers;
    public BioNetwork bioNetwork;
    public BufferedWriter f3;
    public HashMap <String, String[]> list_lineInFile= new HashMap<String, String[]>(); //input file after formating and filtering
    public HashMap<String, String[]> list_unmappedMetabolites = new HashMap<String, String[]>(); //list of non-mapped metabolites
    public Set <BioPhysicalEntity> list_mappingBpe = new HashSet<BioPhysicalEntity>(); //list of mapped metabolites used for analysis
    //Set type is used to avoid metabolites duplicates (no need to used Set now, could be refactored)
    public List <MappingElement> list_mappingElement = new ArrayList<MappingElement>(); //list of mapped metabolites used only for writing mapping output into a file
    public List<HashMap<BioPathway, Double>> list_pathwayEnr = new ArrayList<HashMap<BioPathway, Double>>(); //list of pathway containing mapped metabolites, p-value and corrections values


    public MetExplore4Galaxy(BioNetwork bioNetwork, String inputFile, String outputFileMapping, String outputFilePathEnr, String outputFileInfo, int nameColumn, int chebiColumn, int inchiColumn, int idSBMLColumn, int filteredColumn, String[] inchiLayers) throws IOException{
       this.bioNetwork=bioNetwork;
       this.inputFile=inputFile;
       this.outputFileMapping=outputFileMapping;
       this.outputFilePathEnr=outputFilePathEnr;
       this.outputFileInfo=outputFileInfo;
       this.nameColumn=nameColumn;
       this.chebiColumn=chebiColumn;
       this.inchiColumn=inchiColumn;
       this.idSBMLColumn=idSBMLColumn;
       this.filteredColumn=filteredColumn;
       this.inchiLayers=inchiLayers;
       if (outputFileInfo != "") {
           File fo3 = new File(outputFileInfo);
           fo3.createNewFile();
           this.f3 = new BufferedWriter(new FileWriter(fo3));
       }
    }

    public void  extractData() throws IOException {

        Boolean isFiltered = (filteredColumn >= 0 ) ? true : false;
        BufferedReader f = new BufferedReader(new FileReader(new File(inputFile)));
        String line;
        int i = 0;

        while ((line = f.readLine()) != null) {//Loop on each lines from the input file
            if (i != 0) {//exception on header
                String[] tabuledLine = line.replaceAll("\"","").split("\t");//splitting by tabulation
                try {
                    if (isFiltered ==  false || tabuledLine[filteredColumn] != ""){ //optional filtering on a specified column
                        list_lineInFile.put(""+i, tabuledLine);//add to hashmap
                    }
                } catch (ArrayIndexOutOfBoundsException e) {//avoid errors with filtering functionality containing empty values
                }
            }
            i++;
        }
        if (f != null) f.close();
        if (list_lineInFile.size() < 1 ){//no extraction = error generation
            System.err.println("File badly formatted");
            exit(1);
        }
    }

    public void mapping() throws IOException{

            Boolean isMapped = false, isMappedCurrentBpe = false;
            list_unmappedMetabolites = (HashMap<String, String[]>) list_lineInFile.clone();//will contain non-mapped metabolites
            int mappingOccurences, id = 1;
            BioPhysicalEntity mappedBpe;

            //Loop on each metabolite from the input file
            for (String[] lineInFile : list_lineInFile.values()) {
                isMapped = false;
                mappedBpe = new BioPhysicalEntity();
                mappingOccurences = 0;//identification of multiple mapping

                //Loop for each metabolite from the SBML
                for (BioPhysicalEntity bpe : bioNetwork.getPhysicalEntityList().values()) {
                    isMappedCurrentBpe = false;

                    //Mapping on metabolite identifier associated with a bionetwork
                    try {
                        if ((idSBMLColumn > 0) && (bpe.getId().equals(lineInFile[idSBMLColumn])))
                            isMappedCurrentBpe = mapping4ID_INCHI(lineInFile,"ID",bpe);
                    } catch (ArrayIndexOutOfBoundsException e) {//avoid errors with idSBML column containing empty values
                    }

                    //CHEBI mapping
                    if(!isMappedCurrentBpe){
                        if ( mappingPreRequisite(lineInFile,chebiColumn) && !((lineInFile[chebiColumn]).equals(""))) {

                            //Loop on attribute of the metabolite from the SBML
                            for (Map.Entry<String, Set<BioRef>> key : bpe.getRefs().entrySet()) {
                                if (key.getKey().equals("chebi")) {//researching the one positioned on chebi
                                    //Sometimes different CHEBI can be associated
                                    for (BioRef val : key.getValue()) {
                                        if (lineInFile[chebiColumn].equals(val.id)) {
                                            //add a mappingElement to the List
                                            addMappingElement2List(lineInFile,bpe,chebiColumn,val.id);
                                            isMappedCurrentBpe = true; break;
                                        }
                                    } break;
                                }
                            }
                        }
                    }

                    //InChI mapping
                    if (!isMappedCurrentBpe) {
                        //Call to the "equal" function of the InChI4Galaxy class (allow mapping on selected layers functionality)
                            if ((inchiColumn > 0) && (new InChI4Galaxy(bpe.getInchi(), inchiLayers)).equals(new InChI4Galaxy(lineInFile[inchiColumn], inchiLayers)))
                                isMappedCurrentBpe = mapping4ID_INCHI(lineInFile,"InChI",bpe);
                    }

                    if (isMappedCurrentBpe){
                        isMapped = true; mappingOccurences++; mappedBpe = bpe;
                    }

                }

                if (isMapped) {
                    //If there is no doublets, add the mapped metabolite to the mapped metabolite list (mappingList variable) used for pathway enrichment
                    //Warning: in any case, the multiple matches will be written in the mapping output file (mappingElementList variable) (but not used in the analysis)
                    if(mappingOccurences == 1) list_mappingBpe.add(mappedBpe);
                    else {
                        String nameInInputFile = (nameColumn >= 0 ) ? lineInFile[nameColumn]: 
                                (idSBMLColumn >= 0 ) ? lineInFile[idSBMLColumn] :  
                                                        "line n°" + id + " from the input file";
                        System.out.println("###Warning: There is " + mappingOccurences + " possible matches for " + nameInInputFile + ". Please, check the corresponding lines in the mapping output file. These duplicates will be discarded from the pathway analysis.");
                    }

                    //Remove the mapped metabolite from unmapped list
                    list_unmappedMetabolites.remove(""+id);
                }
                id++;
            }
        writeOutputMapping();
    }

    public void addMappingElement2List(String[] lineInFile, BioPhysicalEntity bpe, int mappingColumnInfile, String matchValueSbml){
        String nameInInputFile = (nameColumn >= 0) ? lineInFile[nameColumn] : "";
        list_mappingElement.add(new MappingElement(true,nameInInputFile,bpe.getName(),bpe.getId(),lineInFile[mappingColumnInfile],matchValueSbml));
    }

    public Boolean mappingPreRequisite(String[] lineInFile, int matchingColumnInFile){
        //Discarding mapping on NA and blank values
        return ((matchingColumnInFile >= 0) && !(lineInFile[matchingColumnInFile]).equals("NA") && !(lineInFile[matchingColumnInFile]).equals(""));
    }

    public Boolean mapping4ID_INCHI (String[] lineInFile, String mappingType, BioPhysicalEntity bpe) {

        String matchValueSbml;
        int mappingColumnInfile;

        if (mappingType =="ID"){
            matchValueSbml = bpe.getId(); mappingColumnInfile = idSBMLColumn;
        }else{
            matchValueSbml = bpe.getInchi(); mappingColumnInfile = inchiColumn;
        }

        if (mappingPreRequisite(lineInFile, mappingColumnInfile)) {
            addMappingElement2List(lineInFile, bpe, mappingColumnInfile, matchValueSbml);
            return true;
        }
        return false;
    }

    public void writeOutputMapping() throws IOException{

        File fo1 = new File(outputFileMapping);
        fo1.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo1));

        f.write("Mapped\tName (Input File)\tName (SBML)\tSBML ID\tMatched value (input File)\tMatched value (SBML)\n");

        String printingMessage = (list_lineInFile.size() - list_unmappedMetabolites.size()) + " metabolites have been mapped (on " + list_lineInFile.size() + ").";
        System.out.println(printingMessage);

        //if "writing console output in a file" functionality is activated
        if (outputFileInfo != "") {
            this.f3.write(printingMessage);
        }

        //Add non-mapped metabolites to the mapping output file
        for (String[] unmappedMetabolites : list_unmappedMetabolites.values()) {
            String nameInInputFile = (nameColumn >= 0) ? unmappedMetabolites[nameColumn] :
                    (idSBMLColumn >= 0 ) ? unmappedMetabolites[idSBMLColumn] :
                            unmappedMetabolites[0];
            list_mappingElement.add(new MappingElement(false, nameInInputFile,"", "", "", ""));
        }

        //Sorting by input file metabolites names (and by true/false mapping) and writing in output file
        Collections.sort(list_mappingElement, new MappedComparator());
        for (int i = 0; i < list_mappingElement.size(); i++) {
            MappingElement m = list_mappingElement.get(i);
            f.write(m.isMapped + "\t" + m.inFileName + "\t" + m.sbmlName + "\t" + m.ID + "\t" + m.inFileVal + "\t" + m.sbmlVal + "\n");
        }
        f.close();
    }

    public void pathwayEnrichment() throws IOException{

        System.out.println("Pathway enrichment in progress...");
        PathwayEnrichment enr = new PathwayEnrichment(bioNetwork, list_mappingBpe);
        HashMap<BioPathway, Double> result = enr.computeEnrichment(); //obtaining p-values for mapped pathway

        list_pathwayEnr.add(result);
        list_pathwayEnr.add(enr.bonferroniCorrection(result)); //obtaining Bonferroni q-values
        list_pathwayEnr.add(enr.benjaminiHochbergCorrection(result));//same for Benjamini Hochberg

        String printingMessage = list_pathwayEnr.get(0).size() + " pathways are concerned among the bioNetwork (on " + bioNetwork.getPathwayList().size() + ").";
        System.out.println(printingMessage);

        //if "writing console output in a file" functionality is activated
        if (outputFileInfo != "") {
            this.f3.write("\n" + printingMessage);
            this.f3.close();
        }
        writeOutputPathEnr();
    }

    public void writeOutputPathEnr() throws IOException{

        File fo2 = new File(outputFilePathEnr);
        fo2.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo2));
        List <PathwayEnrichment4Galaxy> list_pathwayEnr4Galaxy = new ArrayList<PathwayEnrichment4Galaxy>(); //list of pathway enrichment instantiation for sorting
        List <String> listPathwayMetabolites = new ArrayList<String>();

        f.write("Pathway enrichment\tFisher's p-value\tBonferroni correction\tBenjamini-Hochberg correction\tMapped metabolites\tNb of mapped\tCoverage (%)\n");
        HashMap<BioPathway, Double> result = list_pathwayEnr.get(0);//get pathway enrichment without corrections
        Iterator itBonCorr = list_pathwayEnr.get(1).values().iterator(); //Create a loop on Bonferonni values
        Iterator itBenHocCorr = list_pathwayEnr.get(2).values().iterator();//Same for Benjamini Hochberg

        for (Map.Entry<BioPathway, Double> pathEnrEntry : result.entrySet()) {//Loop on pathway enrichment without corrections
            listPathwayMetabolites = new ArrayList<String>();
            BioPathway path = pathEnrEntry.getKey();

            int j = 0; //number of mapped metabolites contained in a BioPathway

            //Extracting metabolites from the mapping list contained in BioPathway
            for (BioPhysicalEntity bpe : list_mappingBpe) {
                if (path.getListOfInvolvedMetabolite().containsValue(bpe)) {
                   listPathwayMetabolites.add(bpe.getName());
                   j++;
                }
            }
            Collections.sort(listPathwayMetabolites);
            String coverage = round((double) j / (double) path.getListOfInvolvedMetabolite().size() * (double) 100);
            list_pathwayEnr4Galaxy.add(new PathwayEnrichment4Galaxy(pathEnrEntry.getKey().getName(),pathEnrEntry.getValue(),(double)itBonCorr.next(),(double)itBenHocCorr.next(),listPathwayMetabolites,j,coverage));
        }
        Collections.sort(list_pathwayEnr4Galaxy);
        for (int i=0;i< list_pathwayEnr4Galaxy.size();i++){
            f.write(list_pathwayEnr4Galaxy.get(i).toString());
        }
        if (f != null) {
            f.close();
        }
    }

    public class MappingElement implements Comparable <MappingElement>{

        public Boolean isMapped;
        public String inFileName, sbmlName, ID, inFileVal, sbmlVal;


        public MappingElement(Boolean isMapped, String inFileName, String sbmlName, String ID, String inFileVal, String sbmlVal){
            this.isMapped=isMapped;
            this.inFileName=inFileName;
            this.sbmlName=sbmlName;
            this.ID=ID;
            this.inFileVal=inFileVal;
            this.sbmlVal=sbmlVal;
        }

        public int compareTo(MappingElement m){
            return (this.inFileName).compareToIgnoreCase(m.inFileName);
        }
    }

    public class MappedComparator implements Comparator <MappingElement> {

        public int compare(MappingElement m1, MappingElement m2) {
            if (m1.isMapped == m2.isMapped) return m1.inFileName.compareToIgnoreCase(m2.inFileName);
            else if (m1.isMapped == true && m2.isMapped == false) return -1;
            return 1;
        }
    }

    public class PathwayEnrichment4Galaxy implements Comparable <PathwayEnrichment4Galaxy>{

        public String pathName, mappedMetabolites, coverage;
        public double p_value, q_value_Bonf, q_value_BenHoc;
        public int nb_mapped;

        public PathwayEnrichment4Galaxy(String pathName, double p_value, double q_value_Bonf, double q_value_BenHoc,
                                        List <String> mappedMetabolites, int nb_mapped, String coverage){
            this.pathName=pathName;
            this.p_value=p_value;
            this.q_value_Bonf=q_value_Bonf;
            this.q_value_BenHoc=q_value_BenHoc;
            this.mappedMetabolites=String.join(";",mappedMetabolites);
            this.nb_mapped=nb_mapped;
            this.coverage=coverage;
        }

        public int compareTo(PathwayEnrichment4Galaxy p) {
            return (this.pathName).compareToIgnoreCase(p.pathName);
        }

        public String toString(){
            return (this.pathName + "\t" + removeSciNot(this.p_value) + "\t" + removeSciNot(this.q_value_Bonf) + "\t" + removeSciNot(this.q_value_BenHoc) + "\t" + this.mappedMetabolites.toString() + "\t" + this.nb_mapped + "\t" + this.coverage + "\n");
        }
    }

    public String removeSciNot(double value) {
        String tmp = (String.valueOf(value));
        if(value < 1e-3){
            String[] splitting = tmp.split("E");
            String power = splitting[1];
            if(power.charAt(1) == '0') {//remove the O after the E if exists
                power = power.substring(2,power.length());
            }
            else power = power.substring(1,power.length());
            String[] number = splitting[0].split("\\.");//obtain the integer and the decimal parts
            return "0." + (new String(new char[Integer.parseInt(power)-1]).replace("\0", "0")) + number[0] + number[1];
        }
        return tmp;
    }

    public String round(double value) {
        return String.valueOf((double) Math.round(value * 100) / 100);
    }

    public void runTime(long elapsedTime){
        long min = elapsedTime / 60000000000L;
        long sec = elapsedTime / 1000000000L - (min * 60L);
        System.out.println("Time to run the process : " + min + "min " + sec + "s");
    }

    public void exec(long startTime) throws IOException{
        extractData();
        mapping();
        pathwayEnrichment();
        runTime(System.nanoTime() - startTime);
    }
}