package fr.inra.toulouse.metexplore;

import junit.framework.TestCase;
import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPhysicalEntity;
import parsebionet.io.JSBMLToBionetwork;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

public class Test_MetExplore4Galaxy extends TestCase {

    static MetExplore4Galaxy met = null;
    static BioNetwork bionet = null;
    static HashMap<String, String[]> parsedFile = new HashMap <String, String[]> ();
    static String inputFile = "sacurineVariableMetadataEnhanced.tsv";
    //static String inputFile = dir + "Galaxy15-[Biosigner_Multivariate_Univariate_Multivariate_variableMetadata.tsv].tabular";
    static String dir = "/home/bmerlet/Documents/PathwayEnrichment/";
    static String sbml = dir + "recon2.v03_ext_noCompartment_noTransport_v2.xml";
    static String outputFile = dir + "output.tsv";
    static BufferedWriter dummyFile;


    protected void setUp() throws Exception {

        super.setUp();
        if (met == null) {
            try {
                bionet = (new JSBMLToBionetwork(sbml)).getBioNetwork();
                met = new MetExplore4Galaxy();
                createdDummyFile("Taurine\tCHEBI:15891\tC2H7NO3S\tC(CS(O)(=O)=O)N\tInChI=1S/C2H7NO3S/c3-1-2-7(4,5)6/h1-3H2,(H,4,5,6)\tTaurine\t124,006693\t1,00727647\t125,01396947\tNA\t[(M-H)]-\t1\t0,88\t5\t2,6122895216\t0,6358457794\t0,2434055545\t387859,346882448\t11652,3712684191\t0,0300427755\t0,1234268278\t15891\tNA\ttaurine\n");
                parsedFile = met.extractData(dir + "dummyFile.tsv", false, -1);
            } catch (IOException e) {
                e.printStackTrace();
            }
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
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void testExtractData() {
            HashMap<String, String[]> obtainedParsing = parsedFile;

            String[] expectedArray = {"Taurine","CHEBI:15891","C2H7NO3S","C(CS(O)(=O)=O)N","InChI=1S/C2H7NO3S/c3-1-2-7(4,5)6/h1-3H2,(H,4,5,6)","Taurine","124,006693","1,00727647","125,01396947","NA","[(M-H)]-","1","0,88","5","2,6122895216","0,6358457794","0,2434055545","387859,346882448","11652,3712684191","0,0300427755","0,1234268278","15891","NA","taurine"};
            HashMap<String, String[]> expectedParsing = new HashMap<String, String[]>();
            expectedParsing.put("Taurine", expectedArray);

            String[] expectedLine = expectedParsing.values().iterator().next();
            String[] obtainedLine = obtainedParsing.values().iterator().next();

            for(int i = 0 ; i < obtainedLine.length; i++){
                    assertEquals(expectedLine[i], obtainedLine[i]);
            }
    }

    public void testMappingWin () {
        Set<BioPhysicalEntity> expectedMap = new HashSet();
        expectedMap.add(bionet.getBioPhysicalEntityById("M_taur"));
        Set<BioPhysicalEntity> obtainedMap = met.mapping(bionet,parsedFile,4);
        assertEquals(expectedMap.iterator().next().getName(), obtainedMap.iterator().next().getName());
    }

    public void testMappingFail () {
        createdDummyFile("Pantothenic acid\tCHEBI:7916\tC9H17NO5\tCC(C)(CO)C(O)C(=O)NCCC(O)=O\tInChI=1S/C9H17NO5/c1-9(2,5-11)7(14)8(15)10-4-3-6(12)13/h7,11,14H,3-5H2,1-2H3,(H,10,15)(H,12,13)\tPantothenic acid\t218,102478\t1,00727647\t219,10975447\tNA\t[(M-H)]-\t1\t4,77\t5\t3,5599610222\t0,2982536819\t0,0837800414\t3012837,77207209\t131428,160471926\t0,043622714\t0,5206814567\t7916\tNA\tpantothenic acid");
        try {
            parsedFile = met.extractData(dir + "dummyFile.tsv", false, -1);

            Set<BioPhysicalEntity> expectedMap = new HashSet();
            BioPhysicalEntity bpe = bionet.getBioPhysicalEntityById("M_pnto_R");
            expectedMap.add(bpe);

            Set<BioPhysicalEntity> obtainedMap = met.mapping(bionet, parsedFile, 4);
            String[] valuesList = parsedFile.values().iterator().next();

            assertFalse("InChI's parsed file: " + valuesList[4] + "\n InChI's MetExplore: " + bpe.getInchi(), obtainedMap.size() > 0);
            System.out.println("InChI's parsed file: " + valuesList[4] + "\n InChI's MetExplore: " + bpe.getInchi());
        }
        catch (IOException e){
            e.printStackTrace();
        }
    }
}
