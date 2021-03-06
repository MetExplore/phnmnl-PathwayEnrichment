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

package fr.inra.toulouse.metexplore.omics;

import fr.inra.toulouse.metexplore.biodata.InChI;
import fr.inra.toulouse.metexplore.omicsComponents.MappedEntity;
import parsebionet.biodata.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static java.lang.System.exit;

public class Mapping extends Omics {
    //Maps bio-entities from a fingerprint to a SBML files based on database information among:
    //SBML_ID, InChI, InChIKey, ChEBI, PubChem_ID, SMILES and KEGG_ID

    protected String outFileMapping;
    protected String[] inchiLayers;
    protected ArrayList<String[]> list_unmappedEntities; //list of non-mapped metabolites
    protected List<MappedEntity> list_MappedEntity = new ArrayList<MappedEntity>(); //list of mapped metabolites used only for writing mapping output into a file
    protected int nameMapping;
    protected int weightPrecision;

    //for performMapping function
    protected ArrayList <String> matchedValues;
    protected ArrayList <String> matchedValuesSBML;
    protected Boolean isMappedBpe;

    public Mapping(String logContent, BioNetwork network, ArrayList<String[]> list_fingerprint,
                   String[] inchiLayers, int nameMapping, int weightPrecision, String outFileMapping, String galaxy,
                   int entityType2Map) throws IOException {
        super(logContent, galaxy, list_fingerprint, network, entityType2Map);
        this.weightPrecision = weightPrecision;
        this.inchiLayers = inchiLayers;
        this.outFileMapping = outFileMapping;
        this.nameMapping = nameMapping;
        this.list_unmappedEntities = (ArrayList<String[]>) list_fingerprint.clone();
        if (!this.outFileMapping.equals("")) this.performMapping();
        else this.quickMapping();
    }


    public void performMapping() throws IOException {
        //Performs mapping on InChI, CHEBI, SBML_ID, PubChem_ID, SMILES, KEGG_ID and InChIKey
        //Remove doublets for analysis and prints warnings
        String warningsDoublets = "";

        //Loop on each metabolite from the input file
        for (String[] lineInFile : this.list_fingerprint) {

            Boolean isMapped = false;
            int nbOccurencesBpe = 0;//identification of multiple mapping
            BioEntity mappedBpe = new BioPhysicalEntity();
            //System.out.println(Arrays.toString(lineInFile));

            //Mapping on metabolites
            // Loop for each metabolite from the SBML
            for (BioEntity e : (Collection<BioEntity>) getEntitySetInNetwork().values()) {

                this.matchedValues = new ArrayList<>();
                this.matchedValuesSBML = new ArrayList<>();
                this.isMappedBpe = false;

                //Mapping on metabolite identifier associated with a bionetwork, InChI, SMILES and PubCHEM_ID
                ArrayList<String> associatedValueInSbml = new ArrayList<>(Collections.singletonList(e.getId()));
                ArrayList<Integer> mappingColumnInfile = new ArrayList<>(Collections.singletonList(1));
                if (this.nameMapping > 0) {
                    associatedValueInSbml.add(e.getName());
                    mappingColumnInfile.add(0);
                }
                if (this.entityType2Map == 1) {
                    BioPhysicalEntity bpe = (BioPhysicalEntity) e;
                    bpe.setMolecularWeight();
                    associatedValueInSbml.addAll(Arrays.asList(bpe.getInchi(), bpe.getSmiles(), bpe.getPubchemCID(), bpe.getMolecularWeight()));
                    mappingColumnInfile.addAll(Arrays.asList(2, 4, 5, 10));
                }
                //System.out.println(Arrays.toString(mappingColumnInfile.toArray()));
                //System.out.println(Arrays.toString(associatedValueInSbml.toArray()));
                Iterator colNum = mappingColumnInfile.iterator();
                for (String sbmlVal : associatedValueInSbml) {
                    mapping4AttributesCase(lineInFile, sbmlVal, (int) colNum.next(), e);
                }
                //TODO: mapping on partial string of the name in SBML

                if (this.entityType2Map == 1) {
                    //Mapping on CHEBI, InChIKey or KEGG
                    String[] associatedValueInSbml2 = {"inchikey", "chebi", "kegg.compound", "hmdb", "chemspider", "pubchem.compound"};
                    //TODO?: regex to take account for SBML diversity
                    //TODO: regex for multiple chebi (CHEBI:[1-9]*)
                    int[] mappingColumnInfile2 = {6, 3, 7, 8, 9, 5};
                    for (int i = 0; i < associatedValueInSbml2.length; i++) {
                        mapping4BiorefCase(lineInFile, associatedValueInSbml2[i], mappingColumnInfile2[i], (BioPhysicalEntity) e);
                    }
                }

                if (this.isMappedBpe) {
                    this.list_MappedEntity.add(new MappedEntity(true, lineInFile[0], e.getName(), e.getId(), String.join(";", matchedValues), String.join(";", matchedValuesSBML)));
                    isMapped = true;
                    mappedBpe = e;
                    nbOccurencesBpe++;
                }
            }

            //If there is no doublets, add the mapped metabolite to the mapped metabolite list (mappingList variable) used for pathway enrichment
            //Warning: in any case, the multiple matches will be written in the mapping output file (MappedEntityList variable) (but not used in the analysis)
            //System.out.println(mappedBpe.getName());
            if (isMapped) {
                if (nbOccurencesBpe > 1) {
                    if (!this.outFileMapping.equals("")) {
                        warningsDoublets = "[WARNING] There are " + nbOccurencesBpe + " possible matches for " + lineInFile[0] + ".\n";
                        writeLog(warningsDoublets);
                    }
                } else {
                    this.list_mappedEntities.put(mappedBpe, lineInFile[0]);
                }

                //Remove the mapped metabolite from unmapped list
                this.list_unmappedEntities.remove(lineInFile);
            }
        }

        if (this.list_mappedEntities.size() == 0) {
            String message="";
            if (!warningsDoublets.equals("")) {
                message = "[FATAL] There is multiple possible match for your whole list of metabolites !\n " +
                        "Please, choose the ID of the desired metabolites among those proposed in the output file.\n " +
                        "Then you can re-run the analysis by adding them into a new column of your input dataset and " +
                        "enter the number of this added column into your program settings.";
            }else {
                message = "[FATAL] There is no match for this network !\nCommon mistakes: wrong type of mapping " +
                        "(by default on the SBML ID and the name of the metabolites), wrong number of column from the dataset" +
                        " or wrong type of bioEntity (or bad SBML).\nPlease check your settings and rerun the analysis.";
            }
            sysExit(this.logContent,message,galaxy,20);
        }

        if (!this.outFileMapping.equals("")) {
            if (!warningsDoublets.equals(""))
                writeLog("[WARNING] Please, check the corresponding lines in the mapping output file.\n" +
                        "[WARNING] These duplicates will be discarded from the pathway analysis.\n");
            writeOutputMapping();
        }
    }

    public void mapping4AttributesCase(String[] lineInFile, String associatedValueInSbml, int mappingColumnInfile, BioEntity bpe) {
        //Mapping case for values which are accessible directly through the attributes of BioPhysicalEntity

        Boolean ifEquals = false;
        try {
                if (ifNotBlankValue(lineInFile, mappingColumnInfile)) {
                    //Splitting database ID values from the fingerprint dataset, if there is more than one per case
                    String[] list_id;
                    if (mappingColumnInfile != 2 ) {
                        list_id = lineInFile[mappingColumnInfile].split(";");
                   }else{
                        list_id = new String[]{lineInFile[mappingColumnInfile]};
                    }

                    for (String id : list_id) {

                        if (mappingColumnInfile == 2) {
                            try {
                                ifEquals = (new InChI(((BioPhysicalEntity) bpe).getInchi(), this.inchiLayers)).equals(new InChI(id, this.inchiLayers));
                                //need to print SBML value only if it is a InChI mapping
                                // (layers selection can compare two different string)
                            } catch (NullPointerException e) {
                                ifEquals = false;
                                //System.out.println("#Warning: " + lineInFile[0] + " have encounter an error with an InChI format. Please, check it validity.");
                                //TODO: check that in Fingerprint class
                            }
                        }else if (mappingColumnInfile == 10 && !associatedValueInSbml.equals("NA")) {
                            try {
                                ifEquals = compareMass(Double.parseDouble(id),Double.parseDouble(associatedValueInSbml),weightPrecision);
                            }catch (NumberFormatException e){
                                //catch conversion in double errors if mass column is wrongly set (e.g., by default on 2nd column)
                            }
                        }else if(mappingColumnInfile == 1 && (this.entityType2Map==5 || this.entityType2Map == 4)){
                            ifEquals = compareProtEnz(id,associatedValueInSbml);
                        }else if(mappingColumnInfile == 1 && this.entityType2Map==6){
                            ifEquals = compareGenes(id,associatedValueInSbml);
                        }else{
                            ifEquals = associatedValueInSbml.equals(id);
                        }

                        //Call to the "equal" function of the InChI4Galaxy class (allow mapping on selected layers functionality)
                        if (ifEquals) {
                            this.matchedValuesSBML.add(associatedValueInSbml);
                            this.matchedValues.add(id);
                            this.isMappedBpe = true;
                            return;
                        }
                    }
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            //avoid errors with an only column in the input file containing empty values
        }
    }

    public void mapping4BiorefCase(String[] lineInFile, String associatedValueInSbml,
                                      int mappingColumnInfile, BioPhysicalEntity bpe) {
        //Mapping case for values which needs to be found into the BioRefs of the BioPhysicalEntity

        try {
            if (ifNotBlankValue(lineInFile, mappingColumnInfile)) {
                //Splitting database ID values from the fingerprint dataset, if there is more than one per case
                String[] list_id = lineInFile[mappingColumnInfile].split(";");
                //Loop on attribute of the metabolite from the SBML
                for (Map.Entry<String, Set<BioRef>> key : bpe.getRefs().entrySet()) {
                    if (key.getKey().equals(associatedValueInSbml)) {//researching the one positioned on chebi
                        //Sometimes different values can be associated to one key

                        for (String id : list_id) {
                            for (BioRef val : key.getValue()) {
                                if (id.equals(val.id)) {
                                    this.matchedValuesSBML.add(val.id);
                                    this.matchedValues.add(id);
                                    this.isMappedBpe = true;
                                    return;
                                }
                            }
                        }
                        return; //case: find the key but values did not match
                    }
                }
                return; //case: did not find the key
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            //avoid errors with a column in the input file containing empty values
        }
        return; //case: outOfBounds or blank values
    }

    public Boolean ifNotBlankValue(String[] lineInFile, int mappingColumnInfile) {
        //Test if mapping is allowed with this parameters and discard mapping on NA and blank values
        return (!(lineInFile[mappingColumnInfile]).equals("NA") && !(lineInFile[mappingColumnInfile]).equals(""));
    }

    public static Boolean compareMass(Double query, Double hit, int error) {
        Double res = error * 10e-6 * hit;
        if (query <= hit + res && query >= hit - res) {
            return true;
        }else {
            return false;
        }
    }

    public Boolean compareGenes(String query, String hit){
        if (hit.startsWith("hsa:")) {
            hit = hit.substring(4, hit.length());
            if (Pattern.compile("(.+)\\.\\d$").matcher(query).matches()) {
                query = query.substring(0, query.length() - 2);
            }
        }
        return query.equals(hit);
    }

    public Boolean compareName(String query, String hit){
        if (Pattern.compile("").matcher(query).matches()) {
            hit = hit.substring(4, hit.length());
            if (Pattern.compile("(.+)\\.\\d$").matcher(query).matches()) {
                query = query.substring(0, query.length() - 2);
            }
        }
        return query.equals(hit);
    }

    public Boolean compareProtEnz(String query, String hit) {

        //two format could be encounter in Recon2:
        //(i) first case: where the id from SBML begin by  "_HSA" (e.g. "_HSA:AL038231")
        // the motif (AL038231) would not include the splicing version (e.g., 1, 2 or 3)
        // from the Metexplore website ID (e.g. "_AL038231_1_c")
        //(ii) second one: where it begin by "_" only (e.g. "_514_2")
        //the motif (514_2) would take account to splicing version (e.g. "_514_2_c")

        Matcher m = Pattern.compile("_HSA:(.+)").matcher(hit);
        if (m.matches()) {
            hit = m.group(1);
            m = Pattern.compile("_(.+)(_\\w){2}").matcher(query);
            if (m.matches()) {
                query = m.group(1);
                return query.equals(hit);
            }

        } else {

                if (hit.startsWith("_")) {
                    hit = hit.substring(1, hit.length());
                    m = Pattern.compile("_(.+)_\\w$").matcher(query);
                    if (m.matches()) {
                        query = m.group(1);
                        return query.equals(hit);
                    }
                }

        }
        return false;
    }

    public void writeOutputMapping() throws IOException {
        //Create an output file for mapping containing for each metabolites, the success of the mapping
        // (true/false), their names in the fingerprint dataset and in the SBML, their ID in the SBML,
        // and the corresponding matched  values (in the dataset & in the SBML)

        File fo1 = new File(this.outFileMapping);
        fo1.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo1));
        int nbEntityInNetwork = getEntitySetInNetwork().size();
        int nbMappedMetabolites = this.list_fingerprint.size() - this.list_unmappedEntities.size();
        String coverageInFile = calculPercent(nbMappedMetabolites, this.list_fingerprint.size());
        String coverageSBML = calculPercent(nbMappedMetabolites, nbEntityInNetwork);

        //File header
        //TODO: if(!this.inchiMapping) f.write("Mapped\tName_(Input_File)\tName_(SBML)\tSBML_ID\tMatched_value\n");
        //else
        f.write("Mapped\tName (Fingerprint)\tName (SBML)\tSBML ID\tMatched value (Fingerprint)\tMatched value (SBML)\n");

        //Print on screen and writing in log
        String plural = (nbMappedMetabolites > 1) ? "s have": " has";

        writeLog(nbMappedMetabolites + " " + typeOfMappedEntity + plural + " been mapped on " + this.list_fingerprint.size() + " in the fingerprint dataset ("
                + coverageInFile + "%) and on " + nbEntityInNetwork + " in the network (" +
                coverageSBML + "%).\n");

        //Add non-mapped metabolites to the mapping output file
        for (String[] unmappedMetabolites : this.list_unmappedEntities) {
            this.list_MappedEntity.add(new MappedEntity(false, unmappedMetabolites[0], "", "", "", ""));
        }

        //Sorting by input file metabolites names (and by true/false mapping) and writing in output file
        Collections.sort(this.list_MappedEntity);
        for (MappedEntity m : this.list_MappedEntity) {
            f.write(m.toString());
        }
        f.close();
    }

    public void quickMapping() {
        for (String[] metabolite : this.list_fingerprint) {

            BioEntity entity = (BioEntity) getEntitySetInNetwork().get(metabolite[1]);

            if (entity != null) this.list_mappedEntities.put(entity, metabolite[0]);
            else System.out.println("[WARNING] " + metabolite[0] + " has not been mapped. Please, check the ID: " + metabolite[1] + ".");

        }

        if (this.list_mappedEntities.size() == 0) {
            String message ="[FATAL] No metabolite have been extracted from the network. \nPlease, check the format of database identifier validity.";
            sysExit(this.logContent,message,galaxy,21);
        }
    }
}
