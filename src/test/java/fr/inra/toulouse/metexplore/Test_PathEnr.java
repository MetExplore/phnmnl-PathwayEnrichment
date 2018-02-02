package fr.inra.toulouse.metexplore;

import junit.framework.TestCase;
import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPhysicalEntity;

import java.io.*;

import java.security.Permission;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

public class Test_PathEnr extends TestCase {
    protected static String outputMappingFile, separator;
    protected static int filteredColumn;
    protected static Boolean ifHeader, ifGalaxy;
    protected static String[] inchiLayers = {"c", "h", ""};
    protected static int[] mappingColumn = {-1, 4, -1, -1, -1, -1, -1, -1, -1, -1};
    protected static List<BioPhysicalEntity> expectedMap = new ArrayList<BioPhysicalEntity>();
    //curl  http://metexplore.toulouse.inra.fr:8080/metExploreWebService/biosources/3223
    protected static BioNetwork network = (new JSBML2Bionetwork4Galaxy("data/recon2.v03_ext_noCompartment_noTransport_v2.xml")).getBioNetwork();
    protected static Fingerprint fingerprint;
    protected File file;
    protected static Mapping mapping;
    protected static fr.inra.toulouse.metexplore.PathwayEnrichment pathEnr;
    //protected BioNetwork network = null;

    public void setUp() throws Exception {
        //Initialization of the parameters before each tests
        super.setUp();
        this.setSecurityManager();
        if (this.file != null) this.file.delete();
        this.setDefaultInChILayers();
        this.setDefaultMappingColumn();
        this.fingerprint = null;
        this.mapping = null;
        this.pathEnr = null;
        this.separator="\t";
        this.ifHeader=true;
        this.ifGalaxy = false;
        this.outputMappingFile = "output.tsv";
        this.filteredColumn = -1;
    }

    /*******************************
     *             SETTINGS
     *****************************/

    public static void setSecurityManager(){
        System.setSecurityManager(new SecurityManager() {

            @Override
            public void checkPermission(Permission perm)
            {
            }

            @Override
            public void checkExit(int status)
            {
                throw new SecurityException();
            }

        });
    }

    /**********setInChILayers**************/

    public void setDefaultInChILayers() {
        this.setInChILayers("c", "h", "");
    }

    public void setInChILayers(String c, String h, String p) {
        this.inchiLayers[0] = c;
        this.inchiLayers[1] = h;
        this.inchiLayers[2] = p;
    }

    /************setMappingColumn************/

    public void setDefaultMappingColumn() {
        this.setMappingColumn(-1, 4, -1);
    }

    public void setMappingColumn(int idColumn, int inchiColumn, int chebiColumn) {
        this.mappingColumn[0] = idColumn;
        this.mappingColumn[1] = inchiColumn;
        this.mappingColumn[2] = chebiColumn;
    }

    public void setMappingColumn2(int mappingType, int mappingColumn) {
        this.mappingColumn[1] = -1;
        this.mappingColumn[mappingType] = mappingColumn;
    }

    public void setMappingAllColumn() {
        for (int i = 0; i < 10; i++) {
            this.mappingColumn[i] = i + 1;
        }
    }

    /************setMapping************/

    public void setMapping(String bpe) {
        this.expectedMap = new ArrayList<BioPhysicalEntity>();
        this.expectedMap.add(this.network.getBioPhysicalEntityById(bpe));
        try {
            this.mapping = new Mapping(this.network, this.fingerprint.list_metabolites, this.inchiLayers, this.outputMappingFile, this.ifGalaxy);
            this.file = new File(this.outputMappingFile);
            assertEquals(this.expectedMap.iterator().next().getName(), this.mapping.list_mappedMetabolites.iterator().next().getName());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void setMapping4OneColumnFile(String inputLine, String bpe){
        this.createDummyFileWithOnlyColumn(inputLine);
        this.setMapping(bpe);
    }

    public void setMapping4MultipleColumnFile(String inputLine, String bpe) {
        this.createDummyFileWithMultipleColumns(inputLine);
        this.setMapping(bpe);
    }

    /***********createDummyFile*************/

    public void createDummyFile(String inputLine, String header){
        //Create a dummy fingerprint dataset for tests

        try {
            this.file = new File("dummyFile.tsv");
            this.file.createNewFile();
            BufferedWriter dummyFile = new BufferedWriter(new FileWriter(this.file));
            dummyFile.write(header + "\n");
            dummyFile.write(inputLine);
            dummyFile.close();
            this.fingerprint = new Fingerprint("dummyFile.tsv", this.ifHeader, this.separator,0,this.mappingColumn,this.filteredColumn);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void createDummyFileWithMultipleColumns(String inputLine){
        this.createDummyFile(inputLine, "\"variableMetadata\\tdatabase_identifier\\tchemical_formula\\tsmiles\\tinchi\\tmetabolite_identification\\tmass_to_charge\\tmass_of_proton\\tmass\\tfragmentation\\tmodifications\\tcharge\\tretention_time\\treliability\\tsample_mean\\tsample_sd\\tsample_CV\\tpool_mean\\tpool_sd\\tpool_CV\\tpoolCV_over_sampleCV\\tchebi.id\\tchemspider.id\\tbiodb.compound.name");
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
        this.createDummyFileWithMultipleColumns("nameMetabolite\tidSBML\tinchi\tchebi\tsmiles\tpubchemColum\tinchikeys\tkegg\thmd\tchemspider\tweight");
        String[] expectedLine = {"nameMetabolite","idSBML","inchi","chebi","smiles","pubchemColum","inchikeys","kegg","hmd","chemspider","weight"};
        assertEquals(Arrays.toString(expectedLine), Arrays.toString((this.fingerprint.list_metabolites).values().iterator().next()));
    }

    public void testSeparator() {
        //Test that each possible mapping values are correctly extracted
        this.setMappingAllColumn();
        this.separator=";";
        this.createDummyFileWithMultipleColumns("nameMetabolite;idSBML;inchi;chebi;smiles;pubchemColum;inchikeys;kegg;hmd;chemspider;weight");
        String[] expectedLine = {"nameMetabolite","idSBML","inchi","chebi","smiles","pubchemColum","inchikeys","kegg","hmd","chemspider","weight"};
        assertEquals(Arrays.toString(expectedLine), Arrays.toString((this.fingerprint.list_metabolites).values().iterator().next()));
    }

    public void testHeader() {
        //Test that each possible mapping values are correctly extracted
        this.ifHeader=false;
        this.createDummyFileWithOnlyColumn("M_taur");
        assertTrue(this.fingerprint.list_metabolites.size() == 2);
    }



    public void testFiltered () {
        //Test that line(s) with empty values in the filtered column are discarded from the parsing
        this.filteredColumn=24;
        this.createDummyFileWithMultipleColumns("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\tplsda|randomforest|svm\n" +
                "Pantothenic acid\tCHEBI:7916\tC9H17NO5\tCC(C)(CO)C(O)C(=O)NCCC(O)=O\tInChI=1S/C9H17NO5/c1-9(2,5-11)7(14)8(15)10-4-3-6(12)13/h7,11,14H,3-5H2,1-2H3,(H,10,15)(H,12,13)\tPantothenic acid\t218,102478\t1,00727647\t219,10975447\tNA\t[(M-H)]-\t1\t4,77\t5\t3,5599610222\t0,2982536819\t0,0837800414\t3012837,77207209\t131428,160471926\t0,043622714\t0,5206814567\t7916\tNA\tpantothenic acid");
        assertTrue(this.fingerprint.list_metabolites.size()==1);

        //Without the filtering, test that all the lines have been extracted
        this.filteredColumn=-1;
        this.createDummyFileWithMultipleColumns("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\tplsda|randomforest|svm\n" +
                "Pantothenic acid\tCHEBI:7916\tC9H17NO5\tCC(C)(CO)C(O)C(=O)NCCC(O)=O\tInChI=1S/C9H17NO5/c1-9(2,5-11)7(14)8(15)10-4-3-6(12)13/h7,11,14H,3-5H2,1-2H3,(H,10,15)(H,12,13)\tPantothenic acid\t218,102478\t1,00727647\t219,10975447\tNA\t[(M-H)]-\t1\t4,77\t5\t3,5599610222\t0,2982536819\t0,0837800414\t3012837,77207209\t131428,160471926\t0,043622714\t0,5206814567\t7916\tNA\tpantothenic acid");
        assertTrue(this.fingerprint.list_metabolites.size()==2);
    }

    /************Mapping************/

    public void testMappingInChI () {
        //Test the success of a mapping with on first two layers (c and h) of InChI from a dataset containing multiple columns

        this.setMapping4MultipleColumnFile("Taurine\tNA\tC2H7NO3S\tC(CS(O)(=O)=O)N\tInChI=1S/C2H7NO3S/c3-1-2-7(4,5)6/h1-3H2,(H,4,5,6)\tTaurine\t124,006693\t1,00727647\t125,01396947\tNA\t[(M-H)]-\t1\t0,88\t5\t2,6122895216\t0,6358457794\t0,2434055545\t387859,346882448\t11652,3712684191\t0,0300427755\t0,1234268278\t15891\tNA\ttaurine",
                "M_taur");
    }

    public void testMappingLayerP () {
        //Test the success of a mapping with an extra p layer on InChI String

        //Positive test
        this.setInChILayers("c", "h", "p");
        this.setMapping4MultipleColumnFile("Glyceric acid\tCHEBI:33508\tC3H6O4\tOCC(O)C(O)=O\tInChI=1S/C3H6O4/c4-1-2(5)3(6)7/h2,4-5H,1H2,(H,6,7)/p-1\tGlyceric acid\t105,018901\tNA\t[(M-H)]-\t1\t0,92\t5\t2,9214987592\t0,2421744543\t0,0828939097\t646166,879585113\t24995,1780580671\t0,0386822334\t0,4666474721\t0,0574485787\t-0,0609303111\t-0,1375521553\t0,0010775982\t1\t-0,1318119477\t-0,0796214185\t1,3556765679\t0,6630650171\t-0,053608993",
                "M_glyc_R");

        //Negative test
        this.setDefaultInChILayers();
        try {
            this.setMapping4MultipleColumnFile("Glyceric acid\tCHEBI:33508\tC3H6O4\tOCC(O)C(O)=O\tInChI=1S/C3H6O4/c4-1-2(5)3(6)7/h2,4-5H,1H2,(H,6,7)/p-1\tGlyceric acid\t105,018901\tNA\t[(M-H)]-\t1\t0,92\t5\t2,9214987592\t0,2421744543\t0,0828939097\t646166,879585113\t24995,1780580671\t0,0386822334\t0,4666474721\t0,0574485787\t-0,0609303111\t-0,1375521553\t0,0010775982\t1\t-0,1318119477\t-0,0796214185\t1,3556765679\t0,6630650171\t-0,053608993",
                    "M_glyc_R");
        }catch (SecurityException e){}//Multiple match are expected
    }

    public void testMappingLayerFormula () {
        //Test the success of a mapping with formula layer only on InChI String

        //Positive test
        this.setInChILayers("","","");
        this.setMapping4MultipleColumnFile("Cinnamoylglycine\tCHEBI:68616\tC11H11NO3\tOC(=O)CNC(=O)C=Cc1ccccc1\tInChI=1S/C11H11NO3/c13-10(12-8-11(14)15)7-6-9-4-2-1-3-5-9/h1-7H,8H2,(H,12,13)(H,14,15)/b7-6+\tCinnamoylglycine\t204,065452\t1,00727647\t205,07272847\tNA\t[(M-H)]-\t1\t7,03\t5\t4,0160399219\t0,5133270871\t0,1278192192\t11742041,4996239\t2134365,54261586\t0,1817712484\t1,4220963763\t68616\tNA\tN-cinnamoylglycine",
                "M_5moxact");

        //Negative test
        this.setDefaultInChILayers();
        try{
            this.setMapping4MultipleColumnFile("Cinnamoylglycine\tCHEBI:68616\tC11H11NO3\tOC(=O)CNC(=O)C=Cc1ccccc1\tInChI=1S/C11H11NO3/c13-10(12-8-11(14)15)7-6-9-4-2-1-3-5-9/h1-7H,8H2,(H,12,13)(H,14,15)/b7-6+\tCinnamoylglycine\t204,065452\t1,00727647\t205,07272847\tNA\t[(M-H)]-\t1\t7,03\t5\t4,0160399219\t0,5133270871\t0,1278192192\t11742041,4996239\t2134365,54261586\t0,1817712484\t1,4220963763\t68616\tNA\tN-cinnamoylglycine",
                    "M_5moxact");
        }catch (SecurityException e){}//No match is expected
    }

    public void testMappingCHEBI () {
        //Test the success of a mapping with CHEBI
        this.setMappingColumn(-1,-1,1);
        this.setMapping4MultipleColumnFile("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tNA\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\n",
                "M_tststeroneglc");
    }

    public void testMappingID () {
        //Test the success of a mapping with the ID of the network from a dataset containing an only column
        this.setMappingColumn(0,-1,-1);
        this.setMapping4OneColumnFile("M_taur", "M_taur");
    }

    /*************Writing output***********/

    public void testWriteOutputMapping() {
        //Test the expected format of the output file obtained by mapping
        this.setMapping4MultipleColumnFile("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\n",
                "M_tststeroneglc");
        try {
            BufferedReader bo = new BufferedReader(new FileReader(this.outputMappingFile));
            bo.readLine();
            assertEquals(bo.readLine(), "true\tTestosterone glucuronide\ttestosterone 3-glucosiduronic acid\tM_tststeroneglc\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void testWriteOutputPathEnr() {
        //Test the expected format of the output file obtained by pathway enrichment
        this.setMapping4MultipleColumnFile("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\n",
                "M_tststeroneglc");
        try {
            this.pathEnr = new fr.inra.toulouse.metexplore.PathwayEnrichment(this.network,this.mapping.list_mappedMetabolites, "pathwayEnr.tsv",this.ifGalaxy);
            this.file = new File("pathwayEnr.tsv");
            BufferedReader bo = new BufferedReader(new FileReader(this.file));
            bo.readLine();
            assertEquals(bo.readLine(), "Steroid metabolism\t0.01805225653206651\t0.01805225653206651\t0.01805225653206651\ttestosterone 3-glucosiduronic acid\tM_tststeroneglc\t1\t1.67");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
