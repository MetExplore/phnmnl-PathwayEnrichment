package fr.inra.toulouse.metexplore;

import parsebionet.biodata.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import static java.lang.System.exit;

public class Mapping extends Omics {
    //Performs mapping on InChI, CHEBI, SBML_ID, PubChem_ID, SMILES, KEGG_ID and InChIKey

    protected String outFileMapping;
    protected String[] inchiLayers;
    protected HashMap<String, String[]> list_unmappedEntities; //list of non-mapped metabolites
    protected List<MappingElement> list_mappingElement = new ArrayList<MappingElement>(); //list of mapped metabolites used only for writing mapping output into a file
    protected BioEntity mappedBpe = new BioPhysicalEntity();

    public Mapping(BioNetwork network, HashMap<String, String[]> list_fingerprint,
                   String[] inchiLayers, String outFileMapping, Boolean ifGalaxy,
                   int bioEntityType) throws IOException {
        super(ifGalaxy, list_fingerprint, network, bioEntityType);
        this.inchiLayers = inchiLayers;
        this.outFileMapping = outFileMapping;
        this.list_unmappedEntities = (HashMap<String, String[]>) list_fingerprint.clone();
        //will contain non-mapped metabolites

        if (this.outFileMapping != "") this.performMapping();
        else this.quickMapping();
    }

    public Boolean mapEntity(String[] lineInFile, int typeMapping, Collection entitySet){
        Iterator it = entitySet.iterator();
        int mappingColumnInFile;
        String valueInSbml;

        while(it.hasNext()) {
            BioEntity e = (BioEntity)it.next();

            if(typeMapping == 1){
                mappingColumnInFile = 1;
                valueInSbml = e.getId();
            }else{
                mappingColumnInFile = 0;
                valueInSbml = e.getName();
            }

            if(lineInFile[mappingColumnInFile].equals(valueInSbml)){
                this.mappedBpe = e;
                return mapping4AttributesCase(lineInFile,valueInSbml,mappingColumnInFile,e);
            }//TODO: mapping on partial string of the name in SBML
        }
        return false;
    }

    public void performMapping() throws IOException {
        //Performs mapping on InChI, CHEBI, SBML_ID, PubChem_ID, SMILES, KEGG_ID and InChIKey
        //Remove doublets for analysis and prints warnings
        Boolean isMapped = false, isMappedCurrentBpe = false;
        int mappingOccurrences = 0, id = 1;
        String warningsDoublets = "";

        //Loop on each metabolite from the input file
        for (String[] lineInFile : this.list_fingerprint.values()) {
            isMapped = false;
            mappingOccurrences = 0;//identification of multiple mapping
            //System.out.println(Arrays.toString(lineInFile));

            //Mapping on other bioEntity than mapping
            Collection entitySet = this.methods.getEntitySetInNetwork().values();
            if (lineInFile[1] != "") {
                isMapped = mapEntity(lineInFile, 1, entitySet);
            }
            if (!isMapped && lineInFile[0]!="") {
                isMapped = mapEntity(lineInFile, 2, entitySet);
            }

            //Mapping on metabolites
            // Loop for each metabolite from the SBML
            if (this.bioEntityType == 1 && !isMapped) {
                for (BioPhysicalEntity bpe : this.network.getPhysicalEntityList().values()) {
                    isMappedCurrentBpe = false;

                    //Mapping on metabolite identifier associated with a bionetwork, InChI, SMILES and PubCHEM_ID
                    String[] associatedValueInSbml = {bpe.getInchi(), bpe.getSmiles(), bpe.getPubchemCID(), bpe.getMolecularWeight()};
                    int[] mappingColumnInfile = {2, 4, 5, 10};
                    for (int i = 0; i < associatedValueInSbml.length; i++) {
                        if (!isMappedCurrentBpe) {
                            isMappedCurrentBpe = mapping4AttributesCase(lineInFile, associatedValueInSbml[i], mappingColumnInfile[i], bpe);
                        }
                    }

                    //Mapping on CHEBI, InChIKey or KEGG
                    String[] associatedValueInSbml2 = {"chebi", "inchikey", "kegg.compound", "hmdb", "chemspider"};
                    //TODO?: regex to take account for SBML diversity
                    int[] mappingColumnInfile2 = {3, 6, 7, 8, 9};
                    for (int i = 0; i < associatedValueInSbml2.length; i++) {
                        if (!isMappedCurrentBpe) {
                            isMappedCurrentBpe = mapping4BiorefCase(lineInFile, associatedValueInSbml2[i], mappingColumnInfile2[i], bpe);
                        }
                    }

                    //Updating variables for the end of the loop
                    if (isMappedCurrentBpe) {
                        isMapped = true;
                        mappingOccurrences++;
                        this.mappedBpe = bpe;
                    }
                }
            }

            if(isMapped){
                //If there is no doublets, add the mapped metabolite to the mapped metabolite list (mappingList variable) used for pathway enrichment
                //Warning: in any case, the multiple matches will be written in the mapping output file (mappingElementList variable) (but not used in the analysis)
                //System.out.println(mappedBpe.getName());
                if (mappingOccurrences <= 1) {
                    this.list_mappedEntities.add(this.mappedBpe);
                } else {
                    if (this.outFileMapping != "") {
                        warningsDoublets = "###Warning: There are " + mappingOccurrences + " possible matches for " + lineInFile[0] + ".\n";
                        writeLog(warningsDoublets);
                    }
                }
                //Remove the mapped metabolite from unmapped list
                this.list_unmappedEntities.remove("" + id);
            }
            id++;
        }
        if (this.outFileMapping != "") {
            if (warningsDoublets != "")
                writeLog("###Warning: Please, check the corresponding lines in the mapping output file. " +
                        "These duplicates will be discarded from the pathway analysis.\n");
            writeOutputMapping();
            writeOutputInfo();
        }
        if (this.list_mappedEntities.size() == 0) {
            if (warningsDoublets != "")
                System.err.println("There is multiple possible match for your whole list of metabolites. " +
                        "Please choose the ID of the desired metabolites among those proposed in the output file. " +
                        "Then you can re-run the analysis by adding them into a new column of your input dataset and " +
                        "enter the number of this added column into your program settings.");
            else System.err.println("There is no match for this network. \nCommon mistakes: wrong type of mapping " +
                    "(by default on InChI only), wrong number of column from the dataset or wrong type of bioEntity." +
                    " Please check your settings" +
                    " and rerun the analysis.");
            exit(1);
        }
    }

    public Boolean mapping4AttributesCase(String[] lineInFile, String associatedValueInSbml, int mappingColumnInfile, BioEntity bpe) {
        //Mapping case for values which are accessible directly through the attributes of BioPhysicalEntity

        Boolean ifEquals;

        try {
            if (ifNotBlankValue(lineInFile, mappingColumnInfile)) {
                ifEquals = (mappingColumnInfile == 2) ?
                        (new InChI4Galaxy(((BioPhysicalEntity)bpe).getInchi(), this.inchiLayers)).equals(new InChI4Galaxy(lineInFile[2], this.inchiLayers))
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

    public Boolean mapping4BiorefCase(String[] lineInFile, String associatedValueInSbml,
                                      int mappingColumnInfile, BioPhysicalEntity bpe) {
        //Mapping case for values which needs to be found into the BioRefs of the BioPhysicalEntity

        try {
            if (ifNotBlankValue(lineInFile, mappingColumnInfile)) {

                //Loop on attribute of the metabolite from the SBML
                for (Map.Entry<String, Set<BioRef>> key : bpe.getRefs().entrySet()) {
                    if (key.getKey().equals(associatedValueInSbml)) {//researching the one positioned on chebi
                        //Sometimes different values can be associated to one key
                        for (BioRef val : key.getValue()) {
                            if (lineInFile[mappingColumnInfile].equals(val.id)) {
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

    public Boolean ifNotBlankValue(String[] lineInFile, int mappingColumnInfile) {
        //Test if mapping is allowed with this parameters and discard mapping on NA and blank values
        return (!(lineInFile[mappingColumnInfile]).equals("NA") && !(lineInFile[mappingColumnInfile]).equals(""));
    }

    public void addMappingElement2List(String[] lineInFile, BioEntity bpe, int mappingColumnInfile, String associatedValueInSbml) {
        //Create a mappingElement and add it to the mapped metabolites list
        this.list_mappingElement.add(new MappingElement(true, lineInFile[0], bpe.getName(), bpe.getId(), lineInFile[mappingColumnInfile], associatedValueInSbml));
        //TODO: write associated SBML value only with InChI mapping
    }

    public void writeOutputMapping() throws IOException {
        //Create an output file for mapping containing for each metabolites, the success of the mapping
        // (true/false), their names in the fingerprint dataset and in the SBML, their ID in the SBML,
        // and the corresponding matched  values (in the dataset & in the SBML)

        File fo1 = new File(this.outFileMapping);
        fo1.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo1));
        int nbEntityInNetwork = this.methods.getEntitySetInNetwork().size();
        int nbMappedMetabolites = this.list_fingerprint.size() - this.list_unmappedEntities.size();
        String coverageInFile = calculPercent(nbMappedMetabolites, this.list_fingerprint.size());
        String coverageSBML = calculPercent(nbMappedMetabolites, nbEntityInNetwork);

        //File header
        f.write("Mapped\tName_(Input_File)\tName_(SBML)\tSBML_ID\tMatched_value_(Input_File)\tMatched_value_(SBML)\n");

        //Print on screen and writing in log
        writeLog(nbMappedMetabolites + " metabolites have been mapped on " + this.list_fingerprint.size() + " in the fingerprint dataset ("
                + coverageInFile + "%) and on " + nbEntityInNetwork + " in the network (" +
                coverageSBML + "%).\n");

        //Add non-mapped metabolites to the mapping output file
        for (String[] unmappedMetabolites : this.list_unmappedEntities.values()) {
            this.list_mappingElement.add(new MappingElement(false, unmappedMetabolites[0], "", "", "", ""));
        }

        //Sorting by input file metabolites names (and by true/false mapping) and writing in output file
        Collections.sort(this.list_mappingElement);
        for (int i = 0; i < this.list_mappingElement.size(); i++) {
            MappingElement m = this.list_mappingElement.get(i);
            f.write(m.toString());
        }
        f.close();
    }

    public String calculPercent(int numerator, int denominator) {
        this.writingComportment.round(1.124);
        return this.writingComportment.round(((double) numerator / (double) denominator) * 100);
    }

    public void quickMapping() {
        for (String[] metabolite : this.list_fingerprint.values()) {
            BioEntity entity=(BioEntity) methods.getEntitySetInNetwork().get(metabolite[1]);

            if (entity != null) this.list_mappedEntities.add(entity);
            else System.out.println("##Warning: " + metabolite[1] + " has not been mapped.");
        }

        if (this.list_mappedEntities.size() == 0) {
            System.err.println("No metabolite have been extracted from the network. \nCheck the list the format of the list of ID provided by the mapping module.");
            exit(1);
        }
    }
}
