package fr.inra.toulouse.metexplore;

import junit.framework.TestCase;
import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPathway;
import parsebionet.biodata.BioPhysicalEntity;
import parsebionet.io.JSBMLToBionetwork;
import parsebionet.statistics.PathwayEnrichment;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

public class Test_MetExplore4Galaxy extends TestCase {

    static String dir = "/home/bmerlet/Documents/PathwayEnrichment/";
    static String sbml = dir + "recon2.v03_ext_noCompartment_noTransport_v2.xml";
    static String outputFile = dir + "output4.tsv";
    static MetExplore4Galaxy met = null;
    static BioNetwork bionet = null;
    static HashMap<String, String[]> parsedFile = new HashMap <String, String[]> ();
    static BufferedWriter dummyFile;
    static Set<BioPhysicalEntity> map = new HashSet();
    static Set<BioPhysicalEntity> expectedMap = new HashSet();

    protected void setUp() throws Exception {

        super.setUp();
        if (met == null) {
            bionet = (new JSBMLToBionetwork(sbml)).getBioNetwork();
            met = new MetExplore4Galaxy();
            createdDummyFile("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\n");
        }
    }

    protected void createdDummyFile(String inputLine){
        try {
            File f = new File(dir + "dummyFile.tsv");
            f.createNewFile();
            dummyFile = new BufferedWriter(new FileWriter(f));
            dummyFile.write("variableMetadata\tdatabase_identifier\tchemical_formula\tsmiles\tinchi\tmetabolite_identification\tmass_to_charge\tmass_of_proton\tmass\tfragmentation\tmodifications\tcharge\tretention_time\treliability\tsample_mean\tsample_sd\tsample_CV\tpool_mean\tpool_sd\tpool_CV\tpoolCV_over_sampleCV\tchebi.id\tchemspider.id\tbiodb.compound.name\n");
            dummyFile.write(inputLine);
            dummyFile.close();
            parsedFile = met.extractData(dir + "dummyFile.tsv", false, -1);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void testExtractData() {
            HashMap<String, String[]> obtainedParsing = parsedFile;

            String[] expectedArray = {"Testosterone glucuronide","CHEBI:28835","C25H36O8","[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O","InChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1","Testosterone glucuronide","463,2329","1,00727647","464,24017647","NA","[(M-H)]-","1","7,9","4","2,1475578771","0,5701078279","0,265467969","178149,617939526","12351,5841321731","0,0693326445","0,2611714128","28835","NA","testosterone 17-glucosiduronic acid"};
            HashMap<String, String[]> expectedParsing = new HashMap<String, String[]>();
            expectedParsing.put("Testosterone glucuronide", expectedArray);

            String[] expectedLine = expectedParsing.values().iterator().next();
            String[] obtainedLine = obtainedParsing.values().iterator().next();

            for(int i = 0 ; i < obtainedLine.length; i++){
                    assertEquals(expectedLine[i], obtainedLine[i]);
            }
    }

    public void setMappingComparaison (String bpe){
        expectedMap.add(bionet.getBioPhysicalEntityById(bpe));
        map = met.mapping(bionet,parsedFile,4);
    }

    public void testMappingWin () {
        setMappingComparaison("M_tststeroneglc");
        assertEquals(expectedMap.iterator().next().getName(), map.iterator().next().getName());
    }

    public void testMappingFail () {
        createdDummyFile("Pantothenic acid\tCHEBI:7916\tC9H17NO5\tCC(C)(CO)C(O)C(=O)NCCC(O)=O\tInChI=1S/C9H17NO5/c1-9(2,5-11)7(14)8(15)10-4-3-6(12)13/h7,11,14H,3-5H2,1-2H3,(H,10,15)(H,12,13)\tPantothenic acid\t218,102478\t1,00727647\t219,10975447\tNA\t[(M-H)]-\t1\t4,77\t5\t3,5599610222\t0,2982536819\t0,0837800414\t3012837,77207209\t131428,160471926\t0,043622714\t0,5206814567\t7916\tNA\tpantothenic acid");
        setMappingComparaison("M_pnto_R");

        BioPhysicalEntity bpe = bionet.getBioPhysicalEntityById("M_pnto_R");
        String[] valuesList = parsedFile.values().iterator().next();
        String msgFail = "InChI's parsed file: " + valuesList[4] + "\n InChI's MetExplore: " + bpe.getInchi();

        assertFalse(msgFail, map.size() > 0);
        System.out.println(msgFail);
    }

    public void obtainedPathwayEnrichment () {
        PathwayEnrichment enr = new PathwayEnrichment(bionet, map);
        ArrayList<HashMap<BioPathway, Double>> resultList = new ArrayList<HashMap<BioPathway, Double>>();

        resultList.add(enr.computeEnrichment());
        for (int i = 0; i < 3; i++) {
            resultList.add(enr.computeEnrichment(i));
        }
        try{
            met.writeOutput(resultList, map, outputFile);
        }catch (IOException e){
            e.printStackTrace();
        }
    }

    public void testWriteOutput() {
        createdDummyFile("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\n");
        setMappingComparaison("M_tststeroneglc");
        obtainedPathwayEnrichment();

        try {
            BufferedReader f = new BufferedReader(new FileReader(new File(outputFile)));
            f.readLine();
            assertEquals(f.readLine(), "Steroid metabolism\t0.02\t0.02\t0.02\t0.02\ttestosterone 3-glucosiduronic acid\t1\t1.67");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}