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
    static String sbml = dir + "recon2.v03_ext_noCompartment_noTransport.xml";
    static String outputFile = dir + "output.tsv";
    static BufferedWriter dummyFile;


    protected void setUp() throws Exception {

        super.setUp();
        if (met == null) {
            try {
                bionet = (new JSBMLToBionetwork(sbml)).getBioNetwork();
                met = new MetExplore4Galaxy();
                createdDummyFile();
                parsedFile = met.extractData(dir + "dummyFile.tsv", false, -1);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    protected void createdDummyFile(){
        try {
            File f = new File(dir + "dummyFile.tsv");
            f.createNewFile();
            dummyFile = new BufferedWriter(new FileWriter(f));
            dummyFile.write("variableMetadata\tdatabase_identifier\tchemical_formula\tsmiles\tinchi\tmetabolite_identification\tmass_to_charge\tmass_of_proton\tmass\tfragmentation\tmodifications\tcharge\tretention_time\treliability\tsample_mean\tsample_sd\tsample_CV\tpool_mean\tpool_sd\tpool_CV\tpoolCV_over_sampleCV\tchebi.id\tchemspider.id\tbiodb.compound.name\n");
            dummyFile.write("Taurine\tCHEBI:15891\tC2H7NO3S\tC(CS(O)(=O)=O)N\tInChI=1S/C2H7NO3S/c3-1-2-7(4,5)6/h1-3H2,(H,4,5,6)\tTaurine\t124,006693\t1,00727647\t125,01396947\tNA\t[(M-H)]-\t1\t0,88\t5\t2,6122895216\t0,6358457794\t0,2434055545\t387859,346882448\t11652,3712684191\t0,0300427755\t0,1234268278\t15891\tNA\ttaurine\n");
            //dummyFile.write("(2-methoxyethoxy)propanoic acid isomer\tCHEBI:67255\tC6H12O4\tCOCCOC(C)C(O)=O\tInChI=1S/C6H12O4/c1-5(6(7)8)10-4-3-9-2/h5H,3-4H2,1-2H3,(H,7,8)\t(2-methoxyethoxy)propanoic acid isomer\t147,065655\t1,00727647\t148,07293147\tNA\t[(M-H)]-\t1\t4,75\t4\t2,1374324398\t0,7609978352\t0,3560336322\t453243,359568937\t34664,1934521402\t0,0764803118\t0,2148120426\t67255\tNA\t2-(2-methoxyethoxy)propanoic acid\n");
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

    public void testMapping () {

        Set<BioPhysicalEntity> expectedMap = new HashSet();
        expectedMap.add(bionet.getBioPhysicalEntityById("M_taur"));
        Set<BioPhysicalEntity> obtainedMap = met.mapping(bionet,parsedFile,4);
        assertEquals(expectedMap.iterator().next().getName(), obtainedMap.iterator().next().getName());
    }
}
