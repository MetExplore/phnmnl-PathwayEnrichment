package fr.inra.toulouse.metexplore;

import junit.framework.TestCase;
import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPhysicalEntity;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.File;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;

import java.util.List;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.ArrayList;

public class Test_MetExplore4Galaxy extends TestCase {

    static String outputFile =  "output1.tsv";
    static MetExplore4Galaxy met = null;
    static BioNetwork bionet = null;
    static File fd;
    static BufferedWriter dummyFile;
    static Set <BioPhysicalEntity> map = new HashSet<BioPhysicalEntity>();
    static List<BioPhysicalEntity> expectedMap = new ArrayList<BioPhysicalEntity>();
    static String[] inchiLayers = new String[8];

    protected void setUp() throws Exception {
        super.setUp();
        if (met == null) {
            inchiLayers[0] = "c";
            inchiLayers[1] = "h";
            bionet = (new JSBML2Bionetwork4Galaxy("data/recon2.v03_ext_noCompartment_noTransport_v2.xml")).getBioNetwork();
            met = new MetExplore4Galaxy(bionet,"dummyFile.tsv", outputFile,outputFile,"",0,1,4,-1, -1,inchiLayers);
        }
    }

    protected void createdDummyFile(String inputLine){
        try {
            fd = new File("dummyFile.tsv");
            fd.createNewFile();
            dummyFile = new BufferedWriter(new FileWriter(fd));
            dummyFile.write("variableMetadata\tdatabase_identifier\tchemical_formula\tsmiles\tinchi\tmetabolite_identification\tmass_to_charge\tmass_of_proton\tmass\tfragmentation\tmodifications\tcharge\tretention_time\treliability\tsample_mean\tsample_sd\tsample_CV\tpool_mean\tpool_sd\tpool_CV\tpoolCV_over_sampleCV\tchebi.id\tchemspider.id\tbiodb.compound.name\n");
            dummyFile.write(inputLine);
            dummyFile.close();
            met.parsedFile = new HashMap<String, String[]>();
            met.extractData();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void testExtractData() {
        createdDummyFile("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tNA\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\n");

        HashMap<String, String[]> obtainedParsing = met.parsedFile;
        String[] expectedArray = {"Testosterone glucuronide","CHEBI:28835","C25H36O8","[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O","NA","Testosterone glucuronide","463,2329","1,00727647","464,24017647","NA","[(M-H)]-","1","7,9","4","2,1475578771","0,5701078279","0,265467969","178149,617939526","12351,5841321731","0,0693326445","0,2611714128","28835","NA","testosterone 17-glucosiduronic acid"};
        HashMap<String, String[]> expectedParsing = new HashMap<String, String[]>();
        expectedParsing.put("Testosterone glucuronide", expectedArray);

        String[] expectedLine = expectedParsing.values().iterator().next();
        String[] obtainedLine = obtainedParsing.values().iterator().next();

        for(int i = 0 ; i < obtainedLine.length; i++){
                assertEquals(expectedLine[i], obtainedLine[i]);
        }
    }

    public void setMappingComparaison (String bpe) {
        expectedMap = new ArrayList<BioPhysicalEntity>();
        expectedMap.add(bionet.getBioPhysicalEntityById(bpe));
        met.mappingList=new HashSet<BioPhysicalEntity>();
        try{
            met.mapping();
            map = met.mappingList;
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void testMappingCHEBIMap () throws SecurityException {
        createdDummyFile("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tNA\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\n");
        setMappingComparaison("M_tststeroneglc");
        assertEquals(expectedMap.iterator().next().getName(), map.iterator().next().getName());
    }

    public void testMappingLayerPSelect () throws SecurityException {
        createdDummyFile("Malic acid\tCHEBI:6650\tC4H6O5\tOC(CC(O)=O)C(O)=O\tInChI=1S/C4H6O5/c5-2(4(8)9)1-3(6)7/h2,5H,1H2,(H,6,7)(H,8,9)\tMalic acid\t133,013511\t1,00727647\t134,02078747\tNA\t[(M-H)]-\t1\t1,08\t5\t2,705885808\t0,2514229882\t0,0929170726\t455610,344205852\t48548,3017483153\t0,1065566275\t1,1467927746\t6650\tNA\tmalic acid");
        inchiLayers[3]="p";
        met.inchiLayers=inchiLayers;
        setMappingComparaison("M_mal_L");
        assertTrue(map.size() == 0);
        inchiLayers[3]="";
        met.inchiLayers=inchiLayers;
        setMappingComparaison("M_mal_L");
        assertFalse(expectedMap.size() == 0);
        assertEquals(expectedMap.iterator().next().getName(), map.iterator().next().getName());
    }

    public void testMappingLayerFormulaSelect () throws SecurityException {
        createdDummyFile("Cinnamoylglycine\tCHEBI:68616\tC11H11NO3\tOC(=O)CNC(=O)C=Cc1ccccc1\tInChI=1S/C11H11NO3/c13-10(12-8-11(14)15)7-6-9-4-2-1-3-5-9/h1-7H,8H2,(H,12,13)(H,14,15)/b7-6+\tCinnamoylglycine\t204,065452\t1,00727647\t205,07272847\tNA\t[(M-H)]-\t1\t7,03\t5\t4,0160399219\t0,5133270871\t0,1278192192\t11742041,4996239\t2134365,54261586\t0,1817712484\t1,4220963763\t68616\tNA\tN-cinnamoylglycine");
        inchiLayers[0]="";
        inchiLayers[1]="";
        met.inchiLayers=inchiLayers;
        setMappingComparaison("M_5moxact");
        assertFalse(expectedMap.size() == 0);
        assertEquals(expectedMap.iterator().next().getName(), map.iterator().next().getName());
        inchiLayers[0] = "c";
        inchiLayers[1] = "h";
        met.inchiLayers=inchiLayers;
        setMappingComparaison("M_5moxact");
        assertTrue(map.size() == 0);
    }

    public void testMappingInChIMap () throws SecurityException {
        createdDummyFile("Taurine\tNA\tC2H7NO3S\tC(CS(O)(=O)=O)N\tInChI=1S/C2H7NO3S/c3-1-2-7(4,5)6/h1-3H2,(H,4,5,6)\tTaurine\t124,006693\t1,00727647\t125,01396947\tNA\t[(M-H)]-\t1\t0,88\t5\t2,6122895216\t0,6358457794\t0,2434055545\t387859,346882448\t11652,3712684191\t0,0300427755\t0,1234268278\t15891\tNA\ttaurine");
        setMappingComparaison("M_taur");
        assertEquals(expectedMap.iterator().next().getName(), map.iterator().next().getName());
    }

    public void testMappingFail () throws SecurityException {
        createdDummyFile("Pantothenic acid\tNA\tC9H17NO5\tCC(C)(CO)C(O)C(=O)NCCC(O)=O\tNA\tPantothenic acid\t218,102478\t1,00727647\t219,10975447\tNA\t[(M-H)]-\t1\t4,77\t5\t3,5599610222\t0,2982536819\t0,0837800414\t3012837,77207209\t131428,160471926\t0,043622714\t0,5206814567\t7916\tNA\tpantothenic acid");
        setMappingComparaison("M_pnto_R");

        BioPhysicalEntity bpe = bionet.getBioPhysicalEntityById("M_pnto_R");
        String[] valuesList = met.parsedFile.values().iterator().next();

        assertTrue(map.size() == 0);
        System.out.println("InChI's parsed file: " + valuesList[4] + "\n InChI's MetExplore: " + bpe.getInchi());
    }

    public void testMappingFiltered () throws SecurityException {
        try{
            createdDummyFile("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\tplsda|randomforest|svm\n" +
                    "Pantothenic acid\tCHEBI:7916\tC9H17NO5\tCC(C)(CO)C(O)C(=O)NCCC(O)=O\tInChI=1S/C9H17NO5/c1-9(2,5-11)7(14)8(15)10-4-3-6(12)13/h7,11,14H,3-5H2,1-2H3,(H,10,15)(H,12,13)\tPantothenic acid\t218,102478\t1,00727647\t219,10975447\tNA\t[(M-H)]-\t1\t4,77\t5\t3,5599610222\t0,2982536819\t0,0837800414\t3012837,77207209\t131428,160471926\t0,043622714\t0,5206814567\t7916\tNA\tpantothenic acid");

            met.filteredColumns=24;
            met.extractData();
            setMappingComparaison("M_tststeroneglc");

            BioPhysicalEntity bpe = bionet.getBioPhysicalEntityById("M_tststeroneglc");
            String[] valuesList = met.parsedFile.values().iterator().next();
            String msg = "InChI's parsed file: " + valuesList[4] + "\n InChI's MetExplore: " + bpe.getInchi();

            assertTrue(msg, map.size() > 0);
            System.out.println(msg);
            met.filteredColumns=-1;
        }catch (IOException e){
            e.printStackTrace();
        }
    }

    public void testWriteOutput() {
        createdDummyFile("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\n");
        setMappingComparaison("M_tststeroneglc");
        try {
            met.pathwayEnrichment();
            File fo = new File(outputFile);
            BufferedReader bo = new BufferedReader(new FileReader(fo));
            bo.readLine();
            assertEquals(bo.readLine(), "Steroid metabolism\t0.01805225653206651\t0.01805225653206651\t0.01805225653206651\ttestosterone 3-glucosiduronic acid\t1\t1.67");
            fo.delete();
        } catch (IOException e) {
            e.printStackTrace();
        }
        fd.delete();
    }
}