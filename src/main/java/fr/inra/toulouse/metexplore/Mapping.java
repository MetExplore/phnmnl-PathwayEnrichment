package fr.inra.toulouse.metexplore;

import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPhysicalEntity;
import parsebionet.biodata.BioRef;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Mapping {
    public String outFileMapping, outFileInfo, text4outputFileInfo="";
    public String[] inchiLayers;
    public BioNetwork network;
    public HashMap<String, String[]> list_fingerprint= new HashMap<String, String[]>(); //input file after formating and filtering
    public HashMap<String, String[]> list_unmappedMetabolites = new HashMap<String, String[]>(); //list of non-mapped metabolites
    public Set<BioPhysicalEntity> list_mappingBpe = new HashSet<BioPhysicalEntity>(); //list of mapped metabolites used for analysis
    //Set type is used to avoid metabolites duplicates (no need to used Set now, could be refactored)
    public List<Mapping.MappingElement> list_mappingElement = new ArrayList<Mapping.MappingElement>(); //list of mapped metabolites used only for writing mapping output into a file

    public Mapping (BioNetwork network, HashMap<String, String[]> list_fingerprint, String outFileInfo, String outFileMapping){
        this.network=network;
        this.list_fingerprint=list_fingerprint;
        this.outFileInfo=outFileInfo;
        this.outFileMapping=outFileMapping;
    }

    public void performMapping() throws IOException {

        Boolean isMapped = false, isMappedCurrentBpe = false;
        list_unmappedMetabolites = (HashMap<String, String[]>) list_fingerprint.clone();//will contain non-mapped metabolites
        int mappingOccurences, id = 1;
        BioPhysicalEntity mappedBpe;
        String warningsDoublets="";

        //Loop on each metabolite from the input file
        for (String[] lineInFile : list_fingerprint.values()) {
            isMapped = false;
            mappedBpe = new BioPhysicalEntity();
            mappingOccurences = 0;//identification of multiple mapping

            //Loop for each metabolite from the SBML
            for (BioPhysicalEntity bpe : network.getPhysicalEntityList().values()) {
                isMappedCurrentBpe = false;

                //Mapping on metabolite identifier associated with a bionetwork
                isMappedCurrentBpe = mapping4ID_INCHI(lineInFile,"ID",bpe);

                //InChI mapping
                if (!isMappedCurrentBpe) {
                    isMappedCurrentBpe = mapping4ID_INCHI(lineInFile,"INCHI",bpe);
                }

                //CHEBI mapping
                if(!isMappedCurrentBpe){
                    try{
                        if (ifNotBlankValue(lineInFile,3)) {

                            //Loop on attribute of the metabolite from the SBML
                            for (Map.Entry<String, Set<BioRef>> key : bpe.getRefs().entrySet()) {
                                if (key.getKey().equals("chebi")) {//researching the one positioned on chebi
                                    //Sometimes different CHEBI can be associated
                                    for (BioRef val : key.getValue()) {
                                        if (lineInFile[3].equals(val.id)) {
                                            //add a mappingElement to the List
                                            addMappingElement2List(lineInFile,bpe,3,val.id);
                                            isMappedCurrentBpe = true; break;
                                        }
                                    } break;
                                }
                            }
                        }
                    }catch (ArrayIndexOutOfBoundsException e) {
                        //avoid errors with an only column in the input file containing empty values
                    }
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
                    warningsDoublets="###Warning: There are " + mappingOccurences + " possible matches for " + lineInFile[0] + ".\n";
                    writeLog(warningsDoublets);
                }

                //Remove the mapped metabolite from unmapped list
                list_unmappedMetabolites.remove(""+id);
            }
            id++;
        }
        if (warningsDoublets!="") writeLog("###Warning: Please, check the corresponding lines in the mapping output file. These duplicates will be discarded from the pathway analysis.\n");
        writeOutputMapping();
    }

    public Boolean mapping4ID_INCHI (String[] lineInFile, String mappingType, BioPhysicalEntity bpe) {

        String matchValueSbml;
        int mappingColumnInfile;
        Boolean ifEquals;

        if (mappingType =="ID"){
            matchValueSbml = bpe.getId();
            mappingColumnInfile = 1;
        }else{
            matchValueSbml = bpe.getInchi();
            mappingColumnInfile = 2;
        }
        try {
            if (ifNotBlankValue(lineInFile, mappingColumnInfile)) {
                ifEquals = (mappingType == "ID") ? bpe.getId().equals(lineInFile[1]) :
                        (new InChI4Galaxy(bpe.getInchi(), inchiLayers)).equals(new InChI4Galaxy(lineInFile[2], inchiLayers));
                //Call to the "equal" function of the InChI4Galaxy class (allow mapping on selected layers functionality)
                if (ifEquals) {
                    addMappingElement2List(lineInFile, bpe, mappingColumnInfile, matchValueSbml);
                    return true;
                }
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            //avoid errors with an only column in the input file containing empty values
        }
        return false;
    }

    public Boolean ifNotBlankValue(String[] lineInFile, int matchingColumnInFile){
        //Test if mapping is allowed with this parameters and discard mapping on NA and blank values
        return ((matchingColumnInFile >= 0) && !(lineInFile[matchingColumnInFile]).equals("NA") && !(lineInFile[matchingColumnInFile]).equals(""));
    }

    public void addMappingElement2List(String[] lineInFile, BioPhysicalEntity bpe, int mappingColumnInfile, String matchValueSbml){
        list_mappingElement.add(new Mapping.MappingElement(true,lineInFile[0],bpe.getName(),bpe.getId(),lineInFile[mappingColumnInfile],matchValueSbml));
    }

    public void writeOutputMapping() throws IOException{

        File fo1 = new File(outFileMapping);
        fo1.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo1));

        f.write("Mapped\tName (Input File)\tName (SBML)\tSBML ID\tMatched value (input File)\tMatched value (SBML)\n");
        writeLog((list_fingerprint.size() - list_unmappedMetabolites.size()) + " metabolites have been mapped (on " + list_fingerprint.size() + ").\n");

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

    public class MappingElement{

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
}