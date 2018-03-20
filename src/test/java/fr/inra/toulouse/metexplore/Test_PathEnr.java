package fr.inra.toulouse.metexplore;

import junit.framework.TestCase;
import org.apache.commons.io.FileUtils;
import parsebionet.biodata.BioEntity;
import parsebionet.biodata.BioNetwork;

import java.io.*;

import java.security.Permission;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

public class Test_PathEnr extends TestCase {
    protected String separator, outputFile, galaxy, logFile="temp/information.txt", dummyFile="temp/dummy.tsv";
    protected int filteredColumn, bioEntityType, entityType2Enrich;
    protected Boolean ifNoHeader, nameMapping, noFormatCheck;
    protected String[] inchiLayers;
    protected int[] mappingColumn;
    protected List<BioEntity> expectedMappedMetabolite;
    // ID biosource:3223
    protected static BioNetwork network = (new JSBML2Bionetwork4Galaxy("data/recon2.02_without_compartment.xml")).getBioNetwork();
    protected Fingerprint fingerprint;
    protected File file;
    protected BufferedReader buffer;
    protected Mapping mapping;
    protected fr.inra.toulouse.metexplore.PathwayEnrichment pathEnr;

    public void setUp() throws Exception {
    //Initialization of the parameters before each tests
        super.setUp();
        //this.setSecurityManager();
        //WritingComportment write = new WritingComportment(this.galaxy);
        (new File("temp")).mkdir();
        this.setDefaultInChILayers();
        this.setDefaultMappingColumn();
        this.noFormatCheck = false;
        this.nameMapping = false;
        this.fingerprint = null;
        this.mapping = null;
        this.pathEnr = null;
        this.separator="\t";
        this.ifNoHeader = false;
        this.galaxy = "";
        this.outputFile = "temp/output.tsv";
        this.filteredColumn = -1;
        this.bioEntityType = 1;
        this.entityType2Enrich = 3;
    }

    protected void tearDown() throws Exception {
        super.tearDown();
        File folder = new File("temp");
        FileUtils.forceDelete(folder);
        folder.delete();
    }

    /*******************************
     *             SETTINGS
     *****************************/

    public static void setSecurityManager(){
    //Allow to catch exit(1) when object is called in a main function (or by other object)
        //but not when running this class as testunit
        System.setSecurityManager(new SecurityManager(){

            @Override
            public void checkPermission(Permission perm) {}

            @Override
            public void checkExit(int status) {
                throw new ThreadDeath();
            }

        });
    }

    public void setBufferReader(String outFile) {
        this.file = new File(outFile);
        try{
            this.buffer = new BufferedReader(new FileReader(this.file));
        }catch (IOException e){
            e.printStackTrace();
        }
    }

    /**********setInChILayers**************/

    public void setDefaultInChILayers() {
        this.setInChILayers("c", "h", "");
    }

    public void setInChILayers(String c, String h, String p) {
        this.inchiLayers = new String[]{c,h,p};
    }

    /************setMappingColumn************/

    public void setDefaultMappingColumn() {
        this.mappingColumn = new int[]{-1, 4, -1, -1, -1, -1, -1, -1, -1, -1};
    }

    public void setMappingColumn(int mappingType, int mappingColumn) {
        this.mappingColumn[1] = -1; //Desactivate inchi mapping by default
        this.mappingColumn[mappingType] = mappingColumn;
        /*mappingType:
        0: id
        1: inchi
        2: chebi
        3: smiles
        4: pubchem
        5: inchikey
        6: kegg
        7: hmdb
        8: chemspider
        9: weight
         */
    }

    public void setMappingAllColumn() {
        for (int i = 0; i < 10; i++) {
            this.mappingColumn[i] = i + 1;
        }
    }

    /************setMapping************/

    public void setMapping(String bpe) {
    //Extract the expected metabolite and process a mapping for comparison

        this.expectedMappedMetabolite = new ArrayList<>();

        try {
            this.mapping = new Mapping(this.network, this.fingerprint.list_entities, this.inchiLayers, this.nameMapping, this.outputFile,
                    this.galaxy, this.bioEntityType);
            this.file = new File(this.outputFile);
            OmicsMethods methods = new OmicsMethods(this.mapping.list_mappedEntities,network,this.bioEntityType);
            this.expectedMappedMetabolite.add((BioEntity)methods.getEntitySetInNetwork().get(bpe));
            assertEquals(this.expectedMappedMetabolite.iterator().next().getName(),
                    this.mapping.list_mappedEntities.keySet().iterator().next().getName());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void setMapping4OneColumnFileByID(String inputLine, String bpe){
        this.setMapping4OneColumnFile(0,0,inputLine,bpe);
    }

    public void setMapping4OneColumnFile(int mappingType, int mappingColumn, String inputLine, String bpe){
        this.setMappingColumn(mappingType,mappingColumn);
        this.createDummyFileWithOnlyColumn(inputLine);
        this.setMapping(bpe);
    }

    public void setMapping4MultipleColumnFile(String inputLine, String bpe) {
        this.createDummyFileWithMultipleColumns(inputLine);
        this.setMapping(bpe);
    }

    public void setMapping4Testo(){
        this.setMapping4MultipleColumnFile("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t" +
                        "[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O" +
                        "\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1" +
                        "\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969" +
                        "\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\n",
                "M_tststeroneglc");
    }

    /************setWriteOutput************/

    public void setBufferTest(String fileName, String header, String line){
        setBufferReader(fileName);
        try {
            if(fileName.equals(this.logFile)) this.pathEnr.write.writeOutputInfo();
            assertEquals(buffer.readLine(), header);
            assertEquals(buffer.readLine(), line);
            assertEquals(buffer.readLine(), null);
        }catch (IOException e){
            e.printStackTrace();
        }
    }

    /***********createDummyFile*************/

    public void createDummyFile(String inputLine, String header){
    //Create a dummy fingerprint dataset for tests

        try {
            this.file = new File(this.dummyFile);
            this.file.createNewFile();
            BufferedWriter dummyFile = new BufferedWriter(new FileWriter(this.file));
            dummyFile.write(header + "\n");
            dummyFile.write(inputLine);
            dummyFile.close();
            this.fingerprint = new Fingerprint(false,this.noFormatCheck,this.dummyFile, this.ifNoHeader, this.separator,";",0,
                    this.mappingColumn,this.inchiLayers,this.filteredColumn);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void createDummyFileWithMultipleColumns(String inputLine){
        this.createDummyFile(inputLine, "\"variableMetadata\\tdatabase_identifier\\tchemical_formula\\tsmiles\\tinchi" +
                "\\tmetabolite_identification\\tmass_to_charge\\tmass_of_proton\\tmass\\tfragmentation\\tmodifications\\tcharge" +
                "\\tretention_time\\treliability\\tsample_mean\\tsample_sd\\tsample_CV\\tpool_mean\\tpool_sd\\tpool_CV\\" +
                "tpoolCV_over_sampleCV\\tchebi.id\\tchemspider.id\\tbiodb.compound.name");
    }

    public void createDummyFileWithOnlyColumn(String inputLine){
        this.createDummyFile(inputLine, "SBML_ID");
    }

    /*******************************
     *              TESTS
     *****************************/

    /*************ExtractData***********/

    public void testExtractData() {
    //Test that each possible mapping values are correctly extracted
        this.setMappingAllColumn();
        this.createDummyFileWithMultipleColumns("nameMetabolite\tidSBML\tinchi\tchebi\tsmiles\tpubchem\tinchikeys\tkegg\thmd\tchemspider\tweight");
        String[] expectedLine = {"nameMetabolite","idSBML","inchi","chebi","smiles","pubchem","inchikeys","kegg","hmd","chemspider","weight"};
        assertEquals(Arrays.toString(expectedLine), Arrays.toString((this.fingerprint.list_entities).iterator().next()));
    }

    public void testSeparator() {
    //Test that each possible mapping values are correctly extracted
        this.setMappingAllColumn();
        this.noFormatCheck=true;
        this.separator=";";
        this.createDummyFileWithMultipleColumns("nameMetabolite;idSBML;inchi;chebi;smiles;pubchem;inchikeys;kegg;hmd;chemspider;weight");
        String[] expectedLine = {"nameMetabolite","idSBML","inchi","chebi","smiles","pubchem","inchikeys","kegg","hmd","chemspider","weight"};
        assertEquals(Arrays.toString(expectedLine), Arrays.toString((this.fingerprint.list_entities).iterator().next()));
    }

    public void testHeader() {
    //Test that each possible mapping values are correctly extracted
        this.ifNoHeader=true;
        this.createDummyFileWithOnlyColumn("M_taur");
        assertTrue(this.fingerprint.list_entities.size() == 2);
    }

    public void testFiltered () {
    //Test that line(s) with empty values in the filtered column are discarded from the parsing
        this.filteredColumn=24;
        this.createDummyFileWithMultipleColumns("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\tplsda|randomforest|svm\n" +
                "Pantothenic acid\tCHEBI:7916\tC9H17NO5\tCC(C)(CO)C(O)C(=O)NCCC(O)=O\tInChI=1S/C9H17NO5/c1-9(2,5-11)7(14)8(15)10-4-3-6(12)13/h7,11,14H,3-5H2,1-2H3,(H,10,15)(H,12,13)\tPantothenic acid\t218,102478\t1,00727647\t219,10975447\tNA\t[(M-H)]-\t1\t4,77\t5\t3,5599610222\t0,2982536819\t0,0837800414\t3012837,77207209\t131428,160471926\t0,043622714\t0,5206814567\t7916\tNA\tpantothenic acid");
        assertTrue(this.fingerprint.list_entities.size()==1);

        //Without the filtering, test that all the lines have been extracted
        this.filteredColumn=-1;
        this.createDummyFileWithMultipleColumns("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\tplsda|randomforest|svm\n" +
                "Pantothenic acid\tCHEBI:7916\tC9H17NO5\tCC(C)(CO)C(O)C(=O)NCCC(O)=O\tInChI=1S/C9H17NO5/c1-9(2,5-11)7(14)8(15)10-4-3-6(12)13/h7,11,14H,3-5H2,1-2H3,(H,10,15)(H,12,13)\tPantothenic acid\t218,102478\t1,00727647\t219,10975447\tNA\t[(M-H)]-\t1\t4,77\t5\t3,5599610222\t0,2982536819\t0,0837800414\t3012837,77207209\t131428,160471926\t0,043622714\t0,5206814567\t7916\tNA\tpantothenic acid");
        assertTrue(this.fingerprint.list_entities.size()==2);
    }

    /************Mapping************/

    /***InChI***/

    public void testMappingInChI () {
    //Test the success of a mapping with on first two layers (c and h) of InChI from a dataset containing multiple columns

        this.setMapping4MultipleColumnFile("Taurine\tNA\tC2H7NO3S\tC(CS(O)(=O)=O)N\tInChI=1S/C2H7NO3S/c3-1-2-7(4,5)6/h1-3H2,(H,4,5,6)" +
                        "\tTaurine\t124,006693\t1,00727647\t125,01396947\tNA\t[(M-H)]-\t1\t0,88\t5\t2,6122895216\t0,6358457794\t0,2434055545" +
                        "\t387859,346882448\t11652,3712684191\t0,0300427755\t0,1234268278\t15891\tNA\ttaurine",
                "M_taur");
    }

    public void testMappingLayerP () {
    //Test the success of a mapping with an extra p layer on InChI String

        //Positive test
        this.setInChILayers("c", "h", "p");
        this.setMapping4MultipleColumnFile("Glyceric acid\tCHEBI:33508\tC3H6O4\tOCC(O)C(O)=O\t" +
                        "InChI=1S/C3H6O4/c4-1-2(5)3(6)7/h2,4-5H,1H2,(H,6,7)/p-1\tGlyceric acid\t105,018901\tNA\t[(M-H)]-" +
                        "\t1\t0,92\t5\t2,9214987592\t0,2421744543\t0,0828939097\t646166,879585113\t24995,1780580671" +
                        "\t0,0386822334\t0,4666474721\t0,0574485787\t-0,0609303111\t-0,1375521553\t0,0010775982\t1\t-0,1318119477" +
                        "\t-0,0796214185\t1,3556765679\t0,6630650171\t-0,053608993",
                "M_glyc_R");

        //Negative test
//        this.setDefaultInChILayers();
//        try {
//            this.setMapping4MultipleColumnFile("Glyceric acid\tCHEBI:33508\tC3H6O4\tOCC(O)C(O)=O\t
// InChI=1S/C3H6O4/c4-1-2(5)3(6)7/h2,4-5H,1H2,(H,6,7)/p-1\tGlyceric acid\t105,018901\tNA\t[(M-H)]-\t1\t0,92\t5\t2,9214987592
// \t0,2421744543\t0,0828939097\t646166,879585113\t24995,1780580671\t0,0386822334\t0,4666474721\t0,0574485787\t-0,0609303111
// \t-0,1375521553\t0,0010775982\t1\t-0,1318119477\t-0,0796214185\t1,3556765679\t0,6630650171\t-0,053608993",
//                    "M_glyc_R");
//        }catch (ThreadDeath e){}//Multiple match are expected
    }

    public void testMappingLayerFormula () {
    //Test the success of a mapping with formula layer only on InChI String

        //Positive test
        this.setInChILayers("","","");
        this.setMapping4MultipleColumnFile("Cinnamoylglycine\tCHEBI:68616\tC11H11NO3\tOC(=O)CNC(=O)C=Cc1ccccc1\t" +
                        "InChI=1S/C11H11NO3/c13-10(12-8-11(14)15)7-6-9-4-2-1-3-5-9/h1-7H,8H2,(H,12,13)(H,14,15)/b7-6+\t" +
                        "Cinnamoylglycine\t204,065452\t1,00727647\t205,07272847\tNA\t[(M-H)]-\t1\t7,03\t5\t4,0160399219" +
                        "\t0,5133270871\t0,1278192192\t11742041,4996239\t2134365,54261586\t0,1817712484\t1,4220963763\t68616" +
                        "\tNA\tN-cinnamoylglycine",
                "M_5moxact");

        //Negative test
//        this.setDefaultInChILayers();
//        try{
//            this.setMapping4MultipleColumnFile("Cinnamoylglycine\tCHEBI:68616\tC11H11NO3\tOC(=O)CNC(=O)C=Cc1ccccc1\tInChI=1S/C11H11NO3/c13-10(12-8-11(14)15)7-6-9-4-2-1-3-5-9/h1-7H,8H2,(H,12,13)(H,14,15)/b7-6+\tCinnamoylglycine\t204,065452\t1,00727647\t205,07272847\tNA\t[(M-H)]-\t1\t7,03\t5\t4,0160399219\t0,5133270871\t0,1278192192\t11742041,4996239\t2134365,54261586\t0,1817712484\t1,4220963763\t68616\tNA\tN-cinnamoylglycine",
//                    "M_5moxact");
//        }catch (ThreadDeath  e){}//No match is expected
    }

    /***Others***/

    public void testMappingCHEBI () {
    //Test the success of a mapping with CHEBI
        this.setMappingColumn(2,1);
        this.setMapping4Testo();
    }

    public void testMappingHMDB () {
        //Test the success of a mapping with the HMDB ID
        this.setMapping4OneColumnFile(7,0,"HMDB00251", "M_taur");
    }

    public void testMappingInchiKey () {
        //Test the success of a mapping with the Inchikey
        this.setMapping4OneColumnFile(5,0,"XOAAWQZATWQOTB-UHFFFAOYSA-N", "M_taur");
    }

    public void testMappingKegg () {
        //Test the success of a mapping with the KEGG ID
        this.setMapping4OneColumnFile(6,0,"C00160", "M_glyclt");
    }

    public void testMappingID () {
        //Test the success of a mapping with the ID of the network from a dataset containing an only column
        this.setMapping4OneColumnFileByID("M_taur", "M_taur");
    }

    public void testMappingName () {
        //Test the success of a mapping with the name of the network from a dataset containing an only column
        this.nameMapping = true;
        this.setMapping4OneColumnFileByID("Taurine", "M_taur");
    }

    public void testMappingIDReaction () {
        //Test the success of a mapping with the ID of a reaction
        this.bioEntityType = 2;
        this.setMapping4OneColumnFileByID("R_FUM", "R_FUM");
    }

    public void testMappingNameReaction () {
        //Test the success of a mapping with the ID of the network from a dataset containing an only column
        this.nameMapping = true;
        this.bioEntityType = 2;
        this.setMapping4OneColumnFileByID("fumarase", "R_FUM");
    }

    public void testMappingIDPathway() {
        //Test the success of a mapping with the ID of a pathway
        this.bioEntityType = 3;
        this.setMapping4OneColumnFileByID("Fatty acid oxidation", "Fatty acid oxidation");
    }

    public void testMappingNamePathway () {
        //Test the success of a mapping with the name of a pathway
        this.nameMapping = true;
        this.bioEntityType = 3;
        this.setMapping4OneColumnFileByID("Fatty acid oxidation", "Fatty acid oxidation");
    }

    /*************Writing output***********/
    public void testWriteOutputMapping(){
        //Test the expected format of the output file obtained by mapping
        this.setMapping4Testo();
        this.setBufferTest(this.outputFile,
                "Mapped\tName_(Input_File)\tName_(SBML)\tSBML_ID\tMatched_value_(Input_File)\tMatched_value_(SBML)",
                "true\tTestosterone glucuronide\ttestosterone 3-glucosiduronic acid\tM_tststeroneglc" +
                        "\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1" +
                        "\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1");
    }

    public void setWriteOutputPathEnr(String pathFile, String typeOfMappedEntity,String typeOfEnrichedEntity, String galColumn, String line) {
    //Test the expected format of the output file obtained by pathway enrichment
        try {
                this.pathEnr = new fr.inra.toulouse.metexplore.PathwayEnrichment(network,this.fingerprint.list_entities,
                        this.mapping.list_mappedEntities, pathFile,this.galaxy,this.bioEntityType, this.entityType2Enrich);
            }catch (IOException e ){
                e.printStackTrace();
            }
            this.setBufferTest(
                    pathFile,
                    typeOfEnrichedEntity + "_name\tp-value\tBonferroni_corrected_p_value\tBH_corrected_p_value\t" +
                            "Mapped_" + typeOfMappedEntity + "_(SBML)\tMapped_" + typeOfMappedEntity + "_(fingerprint)\tMapped_" + typeOfMappedEntity + "_ID\t" +
                            "Nb. of mapped\tCoverage (%)"+galColumn,
                    line);
    }

    public void setGalaxy(){
        this.galaxy=this.logFile;
        this.file = new File(this.galaxy);
        try{
            this.file.createNewFile();
        }catch (IOException e){
            e.printStackTrace();
        }
        WritingComportment.text4outputFileInfo="";
    }

    public void testWriteOutputPathEnr() {
        this.setMapping4Testo();
        this.setWriteOutputPathEnr(
                this.outputFile,
                "metabolites",
                "Pathway",
                "",
                "Steroid metabolism\t0.02314814814814815\t0.02314814814814815\t0.02314814814814815" +
                        "\ttestosterone 3-glucosiduronic acid\tTestosterone glucuronide\tM_tststeroneglc\t1\t1.67");
    }

    public void testWriteOutputReacEnr() {
        this.setMapping4Testo();
        this.entityType2Enrich=2;
        this.setWriteOutputPathEnr(
                this.outputFile,
                "metabolites",
                "Reaction",
                "",
                "UDP-glucuronosyltransferase 1-10 precursor, microsomal\t0.00154320987654321\t0.00154320987654321\t0.00154320987654321\ttestosterone 3-glucosiduronic acid\tTestosterone glucuronide\tM_tststeroneglc\t1\t25.0");
    }

        public void testWriteOutputPathEnr4Galaxy() {
        this.setGalaxy();
        this.setMapping4Testo();
        this.setWriteOutputPathEnr(
                this.outputFile,
                "metabolites",
                "Pathway",
                "\tNb. of unmapped in pathway\tNb. of unmapped in fingerprint\tNb. of remaining in network",
                "Steroid metabolism\t0.02314814814814815\t0.02314814814814815\t0.02314814814814815" +
                        "\ttestosterone 3-glucosiduronic acid\tTestosterone glucuronide\tM_tststeroneglc\t1\t1.67\t59\t0\t2532");
        this.setBufferTest(this.galaxy,
                "1 metabolite has been mapped on 1 in the fingerprint dataset (100.0%) and on 2592 in the network (0.04%).",
                "1 pathway is concerned among the network (on 97 in the network; 1.03%).");
    }

    public void testWriteOutputPathEnrWithReaction() {
        //Test the expected format of the output file obtained by pathway enrichment and with a reaction
        this.bioEntityType = 2;
        this.nameMapping = true;
        this.setMapping4OneColumnFileByID("fumarase", "R_FUM");
        this.setWriteOutputPathEnr(
                this.outputFile,
                "reactions",
                "Pathway",
                "",
                "Citric acid cycle\t0.004750593824228029\t0.004750593824228029\t0.004750593824228029\tfumarase\tfumarase\tR_FUM\t1\t5.0");
    }

    public void testWriteOutput4GalaxyWithReaction() {
        //Test the expected format of the output file obtained by pathway enrichment and with a reaction
        this.setGalaxy();
        this.outputFile="";
        this.bioEntityType = 2;
        this.setMapping4OneColumnFileByID("R_FUM", "R_FUM");
        this.setWriteOutputPathEnr(
                "temp/pathEnr.tsv",
                "reactions",
                "Pathway",
                "\tNb. of unmapped in pathway\tNb. of unmapped in fingerprint\tNb. of remaining in network",
                "Citric acid cycle\t0.004750593824228029\t0.004750593824228029\t0.004750593824228029\tfumarase\tR_FUM\tR_FUM\t1\t5.0\t19\t0\t4190");
        this.setBufferTest(this.galaxy,
                "1 pathway is concerned among the network (on 97 in the network; 1.03%).",
                null);
    }

    /*************Tests with Recon2.02 containing enzymes, protein & genes***********/

    /***Mapping***/

    public void testMappingIDGene() {
        network = new JSBML2Bionetwork4Galaxy("data/recon2.02.xml").getBioNetwork();
        this.bioEntityType = 6;
        this.setMapping4OneColumnFileByID("10026.1", "10026.1");
    }

    public void testMappingNameGene () {
        this.nameMapping = true;
        this.bioEntityType = 6;
        this.setMapping4OneColumnFileByID("10026.1", "10026.1");
    }

    public void testMappingIDEnzyme() {
        this.bioEntityType = 4;
        this.setMapping4OneColumnFileByID("_9415_1_c", "_HSA:9415");
    }

    public void testMappingIDProtein() {
        this.bioEntityType = 5;
        this.setMapping4OneColumnFileByID("_9415_1_c", "_HSA:9415");
    }

   /* public void testMappingNameProtein () {
        this.bioEntityType = 5;
        this.setMapping4OneColumnFileByID("FADS2", "_9415_1");
    }

    public void testMappingNameEnzyme () {
        this.bioEntityType = 4;
        this.setMapping4OneColumnFileByID("FADS2", "_9415_1");
    }*/
    //BUG: there is no name in SBML

    /**Pathway Enrichment**/

    public void testWriteOutputPathEnrWithGene() {
        //Test the expected format of the output file obtained by pathway enrichment and with a reaction
       testMappingNameGene();
       this.setWriteOutputPathEnr(
                this.outputFile,
               "genes",
                "Pathway",
                "",
                "Phosphatidylinositol phosphate metabolism\t0.02823018458197611\t0.02823018458197611\t0.02823018458197611\t10026.1\t10026.1\t10026.1\t1\t1.92");
    }

    public void testWriteOutput4GalaxyWithGene() {
        //Test the expected format of the output file obtained by pathway enrichment and with a reaction
        this.setGalaxy();
        this.outputFile="";
        testMappingNameGene();
        this.setWriteOutputPathEnr(
                "temp/pathEnr.tsv",
                "genes",
                "Pathway",
                "\tNb. of unmapped in pathway\tNb. of unmapped in fingerprint\tNb. of remaining in network",
                "Phosphatidylinositol phosphate metabolism\t0.02823018458197611\t0.02823018458197611\t0.02823018458197611\t10026.1\t10026.1\t10026.1\t1\t1.92\t51\t0\t1790");
        this.setBufferTest(this.galaxy,
                "1 pathway is concerned among the network (on 100 in the network; 1.0%).",
                null);
    }

    public void testWriteOutputPathEnrWithProtein() {
        //Test the expected format of the output file obtained by pathway enrichment and with a reaction
       testMappingIDProtein();
        this.setWriteOutputPathEnr(
                this.outputFile,
                "proteins",
                "Pathway",
                "",
                "Fatty acid synthesis\t0.03203040173724213\t0.03203040173724213\t0.03203040173724213\thsa:9415 (TH)\t_9415_1_c\t_HSA:9415\t1\t1.69");
    }

    public void itestWriteOutput4GalaxyWithProtein() {
        //Test the expected format of the output file obtained by pathway enrichment and with a reaction
        this.setGalaxy();
        this.outputFile="";
        testMappingIDProtein();
        this.setWriteOutputPathEnr(
                "temp/pathEnr.tsv",
                "proteins",
                "Pathway",
                "\tNb. of unmapped in pathway\tNb. of unmapped in fingerprint\tNb. of remaining in network",
                "Fatty acid synthesis\t0.03203040173724213\t0.03203040173724213\t0.03203040173724213\thsa:9415 (TH)\t_9415_1_c\t_HSA:9415\t1\t1.69\t19\t0\t4190");
        this.setBufferTest(this.galaxy,
                "1 pathways are concerned among the network (on 100 in the network).",
                null);
    }

    public void testWriteOutputPathEnrWithEnzyme() {
        //Test the expected format of the output file obtained by pathway enrichment and with a reaction
       testMappingIDEnzyme();
        this.setWriteOutputPathEnr(
                this.outputFile,
                "enzymes",
                "Pathway",
                "",
                "Fatty acid synthesis\t0.02199850857568978\t0.02199850857568978\t0.02199850857568978\thsa:9415 (TH)\t_9415_1_c\t_HSA:9415\t1\t1.69");
    }

    public void itestWriteOutput4GalaxyWithEnzyme() {
        //Test the expected format of the output file obtained by pathway enrichment and with a reaction
        this.setGalaxy();
        this.outputFile="";
        testMappingIDEnzyme();
        this.setWriteOutputPathEnr(
                "temp/pathEnr.tsv",
                "enzymes",
                "Pathway",
                "\tNb. of unmapped in pathway\tNb. of unmapped in fingerprint\tNb. of remaining in network",
                "Fatty acid synthesis\t0.02199850857568978\t0.02199850857568978\t0.02199850857568978\thsa:9415 (TH)\t_9415_1_c\t_HSA:9415\t1\t1.69\t4190");
        this.setBufferTest(this.galaxy,
                "1 pathways are concerned among the network (on 100 in the network).",
                null);
    }
}
