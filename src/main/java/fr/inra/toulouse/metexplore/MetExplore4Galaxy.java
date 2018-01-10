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

    public String inputFile;
    public String outputFile1;
    public String outputFile2;
    public String outputFile3;
    public int nameColumn;
    public int chebiColumn;
    public int inchiColumn;
    public int idSBMLColumn;
    public int filteredColumns;
    public String[] inchiLayers;
    public BioNetwork bn;
    private BufferedWriter f3;
    public HashMap <String, String[]> parsedFile= new HashMap<String, String[]>(); //input file after formating and filtering
    public Set <BioPhysicalEntity> mappingList = new HashSet<BioPhysicalEntity>(); //list of mapped metabolites used for analysis
    //Set is used to avoid metabolites duplicates (no need to used Set now, could be refactored)
    public List <MappingElement> mappingElementList = new ArrayList<MappingElement>(); //list of mapped metabolites used only for writing mapping output into a file
    public HashMap<String, String[]> remainingMetabolites = new HashMap<String, String[]>(); //list of non-mapped metabolites
    public List<HashMap<BioPathway, Double>> pathwayEnrList = new ArrayList<HashMap<BioPathway, Double>>(); //list of pathway containing mapped metabolites, p-value and corrections values


    public MetExplore4Galaxy(BioNetwork bn, String inputFile, String outputFile1, String outputFile2, String outputFile3, int nameColumn, int chebiColumn, int inchiColumn, int idSBMLColumn, int filteredColumns, String[] inchiLayers) throws IOException{
       this.bn=bn;
       this.inputFile=inputFile;
       this.outputFile1=outputFile1;
       this.outputFile2=outputFile2;
       this.outputFile3=outputFile3;
        this.nameColumn=nameColumn;
       this.chebiColumn=chebiColumn;
       this.inchiColumn=inchiColumn;
       this.idSBMLColumn=idSBMLColumn;
       this.filteredColumns=filteredColumns;
       this.inchiLayers=inchiLayers;
       if (outputFile3 != "") {
           File fo3 = new File(outputFile3);
           fo3.createNewFile();
           this.f3 = new BufferedWriter(new FileWriter(fo3));
       }
    }

    public void  extractData() throws IOException {

        Boolean filtered = (filteredColumns >= 0 ) ? true : false;
        BufferedReader f = new BufferedReader(new FileReader(new File(inputFile)));
        String line;
        int i = 0;

        while ((line = f.readLine()) != null) {//Loop on each lines from the input file
            if (i != 0) {//exception on header
                String[] values = line.replaceAll("\"","").split("\t");//splitting by tabulation
                try {
                    if (filtered ==  false || values[filteredColumns] != ""){ //optional filtering on a specified column
                        parsedFile.put(""+i, values);//add to hashmap
                    }
                } catch (ArrayIndexOutOfBoundsException e) {//avoid errors with filtering functionality containing empty values
                }
            }
            i++;
        }
        if (f != null) f.close();
        if (parsedFile.size() < 1 ){//no extraction = error generation
            System.err.println("File badly formatted");
            exit(1);
        }
    }

    public void mapping() throws IOException{

            Boolean mapping = false; Boolean mappingBpe = false;
            remainingMetabolites = (HashMap<String, String[]>) parsedFile.clone();//will contain non-mapped metabolites
            int occurences;
            int i = 1;
            BioPhysicalEntity mappedBpe;

            //Loop on each metabolite from the input file
            for (String[] entry : parsedFile.values()) {
                mapping = false;
                mappedBpe = new BioPhysicalEntity();
                occurences = 0;//identification of multiple mapping

                //Loop for each metabolite from the SBML
                for (BioPhysicalEntity bpe : bn.getPhysicalEntityList().values()) {
                    mappingBpe = false;

                    //Mapping on metabolite identifier in the SBML
                    // Discarding mapping on NA and blank values
                    try {
                        if ((idSBMLColumn >= 0) && !(entry[idSBMLColumn]).equals("NA") && !(entry[idSBMLColumn]).equals("")
                                && (bpe.getId().equals(entry[idSBMLColumn]))) {
                            addMappingElement2List(entry,bpe,idSBMLColumn,bpe.getId());
                            mappingBpe = true; mapping = true; occurences++; mappedBpe = bpe;
                        }
                    } catch (ArrayIndexOutOfBoundsException e) {//avoid errors with idSBML column containing empty values
                    }

                    //Testing CHEBI mapping
                    //Discarding mapping on NA and blank values
                    if(!mappingBpe){
                        if ((chebiColumn >= 0) && !(entry[chebiColumn]).equals("NA") && !((entry[chebiColumn]).equals(""))) {

                            //Loop on attribute of the metabolite from the SBML
                            for (Map.Entry<String, Set<BioRef>> ref : bpe.getRefs().entrySet()) {
                                if (ref.getKey().equals("chebi")) {//researching the one positioned on chebi
                                    //Sometimes different CHEBI can be associated
                                    for (BioRef val : ref.getValue()) {
                                        if (entry[chebiColumn].equals(val.id)) {
                                            //add a mappingElement to the List
                                            addMappingElement2List(entry,bpe,chebiColumn,val.id);
                                            mappingBpe = true; mapping = true; occurences++; mappedBpe = bpe;
                                            break;
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }

                    //Testing InChI mapping (if there is no CHEBI mapping yet)
                    if (!mappingBpe) {
                        //Discarding mapping on NA and blank values
                        //Call to the "equal" function of the InChI4Galaxy class (allow mapping on selected layers functionality)
                        if ((inchiColumn >= 0) && !(entry[inchiColumn]).equals("NA") && !(entry[inchiColumn]).equals("")
                                && (new InChI4Galaxy(bpe.getInchi(), inchiLayers)).equals(new InChI4Galaxy(entry[inchiColumn], inchiLayers))) {
                            addMappingElement2List(entry,bpe,inchiColumn,bpe.getInchi());
                            mapping = true; occurences++; mappedBpe = bpe;
                        }
                    }

                }

                if (mapping) {
                    //If there is no doublets, add the mapped metabolite to the mapped metabolite list (mappingList variable) used for pathway enrichment
                    //Warning: in any case, the multiple matches will be written in the mapping output file (mappingElementList variable) (but not used in the analysis)
                    if(occurences == 1) mappingList.add(mappedBpe);
                    else System.out.println("###Warning: There is " + occurences + " possible matches for " + entry[0] + ". Please, check the corresponding lines in the mapping output file. These duplicates will be discarded from the pathway analysis.");

                    //Remove the mapped metabolite from remaining list
                    remainingMetabolites.remove(""+i);
                }
                i++;
            }
        writeOutputMapping();
    }

    public void addMappingElement2List(String[] entry, BioPhysicalEntity bpe, int matchValueInFile, String matchValueSbml){

        String nameInInputFile = (nameColumn >= 0) ? entry[nameColumn] : "";
        mappingElementList.add(new MappingElement(true,nameInInputFile,bpe.getName(),bpe.getId(),entry[matchValueInFile],matchValueSbml));
    }

    public void writeOutputMapping() throws IOException{

        File fo1 = new File(outputFile1);
        fo1.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo1));

        f.write("Mapped\tName (Input File)\tName (SBML)\tSBML ID\tMatched value (input File)\tMatched value (SBML)\n");

        String printingMessage = (parsedFile.size() - remainingMetabolites.size()) + " metabolites have been mapped (on " + parsedFile.size() + ").";
        System.out.println(printingMessage);

        //if "writing console output in a file" functionality is activated
        if (outputFile3 != "") {
            this.f3.write(printingMessage);
        }

        //Add non-mapped metabolites to the mapping output file
        for (String[] entry : remainingMetabolites.values()) {
            mappingElementList.add(new MappingElement(false, entry[0],"", "", "", ""));
        }

        //Sorting by input file metabolites names (and by true/false mapping) and writing in output file
        Collections.sort(mappingElementList, new MappedComparator());
        for (int i = 0; i < mappingElementList.size(); i++) {
            MappingElement m = mappingElementList.get(i);
            f.write(m.mapped + "\t" + m.inFileName + "\t" + m.sbmlName + "\t" + m.ID + "\t" + m.inFileVal + "\t" + m.sbmlVal + "\n");
        }
        f.close();
    }

    public void pathwayEnrichment() throws IOException{

        System.out.println("Pathway enrichment in progress...");
        PathwayEnrichment enr = new PathwayEnrichment(bn, mappingList);
        HashMap<BioPathway, Double> res = enr.computeEnrichment(); //obtaining p-values for mapped pathway

        pathwayEnrList.add(res);
        pathwayEnrList.add(enr.bonferroniCorrection(res)); //obtaining Bonferroni q-values
        pathwayEnrList.add(enr.benjaminiHochbergCorrection(res));//same for Benjamini Hochberg

        String printingMessage = pathwayEnrList.get(0).size() + " pathways are concerned among the network (on " + bn.getPathwayList().size() + ").";
        System.out.println(printingMessage);

        //if "writing console output in a file" functionality is activated
        if (outputFile3 != "") {
            this.f3.write("\n" + printingMessage);
            this.f3.close();
        }
        writeOutputPathEnr();
    }

    public void writeOutputPathEnr() throws IOException{

        File fo2 = new File(outputFile2);
        fo2.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo2));
        List <PathwayEnrichment4Galaxy> listPathwayEnr = new ArrayList<PathwayEnrichment4Galaxy>(); //list of pathway enrichment instantiation for sorting
        List <String> listPathwayMetabolites = new ArrayList<String>();

        f.write("Pathway enrichment\tFisher's p-value\tBonferroni correction\tBenjamini-Hochberg correction\tMapped metabolites\tNb of mapped\tCoverage (%)\n");
        HashMap<BioPathway, Double> result = pathwayEnrList.get(0);//get pathway enrichment without corrections
        Iterator itBonCorr = pathwayEnrList.get(1).values().iterator(); //Create a loop on Bonferonni values
        Iterator itBenHocCorr = pathwayEnrList.get(2).values().iterator();//Same for Benjamini Hochberg

        for (Map.Entry<BioPathway, Double> entry : result.entrySet()) {//Loop on pathway enrichment without corrections
            listPathwayMetabolites = new ArrayList<String>();
            BioPathway path = entry.getKey();

            int j = 0; //number of mapped metabolites contained in a BioPathway

            //Extracting metabolites from the mapping list contained in BioPathway
            for (BioPhysicalEntity bpe : mappingList) {
                if (path.getListOfInvolvedMetabolite().containsValue(bpe)) {
                   listPathwayMetabolites.add(bpe.getName());
                   j++;
                }
            }
            Collections.sort(listPathwayMetabolites);
            String coverage = round((double) j / (double) path.getListOfInvolvedMetabolite().size() * (double) 100);
            listPathwayEnr.add(new PathwayEnrichment4Galaxy(entry.getKey().getName(),entry.getValue(),(double)itBonCorr.next(),(double)itBenHocCorr.next(),listPathwayMetabolites,j,coverage));
        }
        Collections.sort(listPathwayEnr);
        for (int i=0;i< listPathwayEnr.size();i++){
            f.write(listPathwayEnr.get(i).toString());
        }
        if (f != null) {
            f.close();
        }
    }

    public class MappingElement implements Comparable <MappingElement>{

        public Boolean mapped;
        public String inFileName;
        public String sbmlName;
        public String ID;
        public String inFileVal;
        public String sbmlVal;


        public MappingElement(Boolean mapped, String inFileName, String sbmlName, String ID, String inFileVal, String sbmlVal){
            this.mapped=mapped;
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

    class MappedComparator implements Comparator <MappingElement> {

        public int compare(MappingElement m1, MappingElement m2) {
            if (m1.mapped == m2.mapped) return m1.inFileName.compareToIgnoreCase(m2.inFileName);
            if (m1.mapped == true && m2.mapped == false) return -1;
            return 1;
        }
    }

    public class PathwayEnrichment4Galaxy implements Comparable <PathwayEnrichment4Galaxy>{

        public String pathName;
        public double p_value;
        public double q_value_Bonf;
        public double q_value_BenHoc;
        public String mappedMetabolites;
        public int nb_mapped;
        public String coverage;

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
            String[] splitted = tmp.split("E");
            String power = splitted[1];
            if(power.charAt(1) == '0') {
                power = power.substring(2,power.length());
            }
            else power = power.substring(1,power.length());
            String[] number = splitted[0].split("\\.");
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