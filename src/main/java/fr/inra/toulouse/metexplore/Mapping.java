package fr.inra.toulouse.metexplore;

import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPhysicalEntity;
import parsebionet.biodata.BioRef;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import static java.lang.System.exit;

public class Mapping {
    //Performs mapping on InChI, CHEBI, SBML_ID, PubChem_ID, SMILES, KEGG_ID and InChIKey
    
    public String outFileMapping, outFileInfo, text4outputFileInfo="";
    public String[] inchiLayers;
    public BioNetwork network;
    public HashMap<String, String[]> list_fingerprint= new HashMap<String, String[]>(); //input file after formating and filtering
    public HashMap<String, String[]> list_unmappedMetabolites = new HashMap<String, String[]>(); //list of non-mapped metabolites
    public Set<BioPhysicalEntity> list_mappedMetabolites = new HashSet<BioPhysicalEntity>(); //list of mapped metabolites used for analysis
    //Set type is used to avoid metabolites duplicates (no need to used Set now, could be refactored)
    public List<Mapping.MappingElement> list_mappingElement = new ArrayList<Mapping.MappingElement>(); //list of mapped metabolites used only for writing mapping output into a file
    //public Boolean ifGalaxy = false;
    
    public Mapping (BioNetwork network, HashMap<String, String[]> list_fingerprint, String[] inchiLayers, String outFileMapping,
                    String outFileInfo) throws IOException{
        this.network=network;
        this.list_fingerprint=list_fingerprint;
        this.inchiLayers=inchiLayers;
        this.outFileMapping=outFileMapping;
        this.outFileInfo=outFileInfo;
        
        //if(!this.ifGalaxy) this.performMapping();
        //else this.quickMapping();
    }

    public void performMapping() throws IOException {
        //Performs mapping on InChI, CHEBI, SBML_ID, PubChem_ID, SMILES, KEGG_ID and InChIKey
        //Remove doublets for analysis and prints warnings
        
        //TODO: add name mapping

        Boolean isMapped = false, isMappedCurrentBpe = false;
        list_unmappedMetabolites = (HashMap<String, String[]>) list_fingerprint.clone();//will contain non-mapped metabolites
        int mappingOccurrences, id = 1;
        BioPhysicalEntity mappedBpe;
        String warningsDoublets="";

        //Loop on each metabolite from the input file
        for (String[] lineInFile : list_fingerprint.values()) {
            isMapped = false;
            mappedBpe = new BioPhysicalEntity();
            mappingOccurrences = 0;//identification of multiple mapping

            //Loop for each metabolite from the SBML
            for (BioPhysicalEntity bpe : network.getPhysicalEntityList().values()) {
                isMappedCurrentBpe = false;

                //Mapping on metabolite identifier associated with a bionetwork, InChI, SMILES and PubCHEM_ID
                String[] associatedValueInSbml = {bpe.getId(), bpe.getInchi(), bpe.getSmiles(), bpe.getPubchemCID(), bpe.getMolecularWeight()};
                int[] mappingColumnInfile = {1, 2, 4, 5, 7};
                for (int i = 0; i < associatedValueInSbml.length; i++){
                    if (!isMappedCurrentBpe) {
                        isMappedCurrentBpe = mapping4AttributesCase(lineInFile,associatedValueInSbml[i],mappingColumnInfile[i],bpe);
                    }
                }

                //Mapping on CHEBI, InChIKey or KEGG
                String[] associatedValueInSbml2 = {"chebi", "inchikey", "kegg.compound", "hmdb", "chemspider"};
                //TODO?: regex to take account for SBML diversity
                int[] mappingColumnInfile2 = {3, 6, 8, 9, 10};
                for (int i = 0; i < associatedValueInSbml2.length; i++){
                    if(!isMappedCurrentBpe) {
                        isMappedCurrentBpe = mapping4BiorefCase(lineInFile, associatedValueInSbml2[i], mappingColumnInfile2[i], bpe);
                    }
                }

                //Updating variables for the end of the loop
                if (isMappedCurrentBpe){
                    isMapped = true; mappingOccurrences++; mappedBpe = bpe;
                }
            }

            if (isMapped) {
                //If there is no doublets, add the mapped metabolite to the mapped metabolite list (mappingList variable) used for pathway enrichment
                //Warning: in any case, the multiple matches will be written in the mapping output file (mappingElementList variable) (but not used in the analysis)
                if(mappingOccurrences == 1) list_mappedMetabolites.add(mappedBpe);
                else {
                    if (outFileMapping != "") {
                        warningsDoublets = "###Warning: There are " + mappingOccurrences + " possible matches for " + lineInFile[0] + ".\n";
                        writeLog(warningsDoublets);
                    }
                }

                //Remove the mapped metabolite from unmapped list
                list_unmappedMetabolites.remove(""+id);
            }
            id++;
        }
        if (outFileMapping != "") {
            if (warningsDoublets != "")
                writeLog("###Warning: Please, check the corresponding lines in the mapping output file. These duplicates will be discarded from the pathway analysis.\n");
            writeOutputMapping();
            writeOutputInfo();
        }
        if (list_mappedMetabolites.size() == 0 ){
            System.err.println("There is no match for this network. \nCommon mistakes: mapping parameters ");
            exit(1);
        }
    }

    public Boolean mapping4AttributesCase (String[] lineInFile, String associatedValueInSbml, int mappingColumnInfile, BioPhysicalEntity bpe) {
    //Mapping case for values which are accessible directly through the attributes of 
        // BioPhysicalEntity
        
        Boolean ifEquals;

        try {
            if (ifNotBlankValue(lineInFile, mappingColumnInfile)) {
                ifEquals = (mappingColumnInfile == 2) ?
                    (new InChI4Galaxy(bpe.getInchi(), inchiLayers)).equals(new InChI4Galaxy(lineInFile[2], inchiLayers)) 
                        : associatedValueInSbml.equals(lineInFile[mappingColumnInfile]);                        
                //Call to the "equal" function of the InChI4Galaxy class (allow mapping on selected layers functionality)
                if (ifEquals) {
                    addMappingElement2List(lineInFile, bpe, mappingColumnInfile, associatedValueInSbml);
                    return true;
                }
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            //avoid errors with an only column in the input file containing empty values
        }
        return false;
    }

    public Boolean mapping4BiorefCase (String[] lineInFile, String associatedValueInSbml, int mappingColumnInfile, BioPhysicalEntity bpe) {
    //Mapping case for values which needs to be found into the BioRefs of the BioPhysicalEntity
        
        try {
            if (ifNotBlankValue(lineInFile, mappingColumnInfile)) {

                //Loop on attribute of the metabolite from the SBML
                for (Map.Entry<String, Set<BioRef>> key : bpe.getRefs().entrySet()) {
                    if (key.getKey().equals(associatedValueInSbml)) {//researching the one positioned on chebi
                        //Sometimes different values can be associated to one key
                        for (BioRef val : key.getValue()) {
                            if (lineInFile[3].equals(val.id)) {
                                //add a mappingElement to the list
                                addMappingElement2List(lineInFile, bpe, mappingColumnInfile, val.id);
                                return true;
                            }
                        }
                        return false; //case: find the key but values did not match
                    }
                }
                return false; //case: did not find the key

            }
        } catch (ArrayIndexOutOfBoundsException e) {
            //avoid errors with a column in the input file containing empty values
        }
        return false; //case: outOfBounds or blank values
    }

    public Boolean ifNotBlankValue(String[] lineInFile, int mappingColumnInfile){
        //Test if mapping is allowed with this parameters and discard mapping on NA and blank values
        return ((mappingColumnInfile >= 0) && !(lineInFile[mappingColumnInfile]).equals("NA") && !(lineInFile[mappingColumnInfile]).equals(""));
    }

    public void addMappingElement2List(String[] lineInFile, BioPhysicalEntity bpe, int mappingColumnInfile, String associatedValueInSbml){
        //Create a mappingElement and add it to the mapped metabolites list
        
        list_mappingElement.add(new Mapping.MappingElement(true,lineInFile[0],bpe.getName(),bpe.getId(),lineInFile[mappingColumnInfile],associatedValueInSbml));
    }

    public void writeOutputMapping() throws IOException{
        //Create an output file for mapping containing for each metabolites, the success of the mapping
        // (true/false), their names in the fingerprint dataset and in the SBML, their ID in the SBML,
        // and the corresponding matched  values (in the dataset & in the SBML)
        
        File fo1 = new File(outFileMapping);
        fo1.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo1));
        int nbMappedMetabolites = list_fingerprint.size() - list_unmappedMetabolites.size();
        double coverageInFile = nbMappedMetabolites / list_fingerprint.size();
        double coverageSBML = nbMappedMetabolites / network.getPhysicalEntityList().size();

        //File header
        f.write("Mapped\tName_(Input_File)\tName_(SBML)\tSBML_ID\tMatched_value_(Input_File)\tMatched_value_(SBML)\n");

        //Print on screen and writing in log
        writeLog( nbMappedMetabolites + " metabolites have been mapped on " + list_fingerprint.size() + "in the fingerprint dataset (" + round(coverageInFile) + ") and " + "in the network ( "+ round(coverageSBML) + ").\n");

        //Add non-mapped metabolites to the mapping output file
        for (String[] unmappedMetabolites : list_unmappedMetabolites.values()) {
            list_mappingElement.add(new Mapping.MappingElement(false, unmappedMetabolites[0],"", "", "", ""));
        }

        //Sorting by input file metabolites names (and by true/false mapping) and writing in output file
        Collections.sort(list_mappingElement, new Mapping.MappedComparator());
        for (int i = 0; i < list_mappingElement.size(); i++) {
            Mapping.MappingElement m = list_mappingElement.get(i);
            f.write(m.isMapped + "\t" + m.inFileName + "\t" + m.sbmlName + "\t" + m.ID + "\t" + m.inFileVal + "\t" + m.sbmlVal + "\n");
        }
        f.close();
    }

    public void writeOutputInfo() throws IOException {
        //if (this.ifGalaxy) {//if "writing console output in a file" functionality is activated
            File f = new File(outFileInfo);
            f.createNewFile();
            BufferedWriter b = new BufferedWriter(new FileWriter(f));
            b.write(text4outputFileInfo);
            b.close();
        //}
    }

    public void quickMapping() {

        for (String[] metabolite : list_fingerprint.values()) {
            list_mappedMetabolites.add(network.getBioPhysicalEntityById(metabolite[1]));
        }

        if (list_mappedMetabolites.size() == 0 ){
            System.err.println("There is no match for this network. \nCommon mistakes: mapping parameters ");
            exit(1);
        }
    }

    public class MappingElement{

        public Boolean isMapped;
        public String inFileName, sbmlName, ID, inFileVal, sbmlVal;


        public MappingElement(Boolean isMapped, String inFileName, String sbmlName, String ID, 
                              String inFileVal, String sbmlVal){
            this.isMapped=isMapped;
            this.inFileName=inFileName;
            this.sbmlName=sbmlName;
            this.ID=ID;
            this.inFileVal=inFileVal;
            this.sbmlVal=sbmlVal;
        }
    }

    public class MappedComparator implements Comparator <Mapping.MappingElement> {
        public int compare(Mapping.MappingElement m1, Mapping.MappingElement m2) {
            if (m1.isMapped == m2.isMapped) return m1.inFileName.compareToIgnoreCase(m2.inFileName);
            else if (m1.isMapped == true && m2.isMapped == false) return -1;
            return 1;
        }
    }

    public void writeLog(String message){
        System.out.println(message.replaceAll("\n", ""));
        text4outputFileInfo+=message;
    }

    public String round(double value) {
        return String.valueOf((double) Math.round(value * 100) / 100);
    }
}