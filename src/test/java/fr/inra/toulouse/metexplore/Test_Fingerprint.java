package fr.inra.toulouse.metexplore;

import fr.inra.toulouse.metexplore.io.Fingerprint;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

public class Test_Fingerprint extends Test{


    protected Fingerprint fingerprint;
    protected String checkingFile, separator, dummyFile = "temp/dummy.tsv";
    protected int filteredColumn;
    protected Boolean ifNoHeader, noFormatCheck, layerWarning;
    protected String[] inchiLayers;
    protected int[] mappingColumn;

    public void setUp() throws Exception {
        super.setUp();
        this.fingerprint = null;
        this.noFormatCheck = false;
        this.separator = "\t";
        this.ifNoHeader = false;
        this.layerWarning = false;
        this.checkingFile = "";
        this.filteredColumn = -1;
        this.setDefaultInChILayers();
        this.setDefaultMappingColumn();
    }


    /**********setInChILayers**************/

    public void setDefaultInChILayers() {
        this.setInChILayers("c", "h", "");
    }

    public void setInChILayers(String c, String h, String p) {
        this.inchiLayers = new String[]{c, h, p};
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

    /***********createDummyFile*************/

    public void createDummyFile(String inputLine, String header) {
        //Create a dummy fingerprint dataset for tests

        try {
            this.file = new File(this.dummyFile);
            this.file.createNewFile();
            BufferedWriter dummyFile = new BufferedWriter(new FileWriter(this.file));
            dummyFile.write(header + "\n");
            dummyFile.write(inputLine);
            dummyFile.close();
            this.fingerprint = new Fingerprint(this.logContent,this.layerWarning, this.noFormatCheck, this.checkingFile, this.dummyFile, this.ifNoHeader, this.separator, ";", 0,
                    this.mappingColumn, this.inchiLayers, this.filteredColumn);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void createDummyFileWithMultipleColumns(String inputLine) {
        this.createDummyFile(inputLine, "\"variableMetadata\\tdatabase_identifier\\tchemical_formula\\tsmiles\\tinchi" +
                "\\tmetabolite_identification\\tmass_to_charge\\tmass_of_proton\\tmass\\tfragmentation\\tmodifications\\tcharge" +
                "\\tretention_time\\treliability\\tsample_mean\\tsample_sd\\tsample_CV\\tpool_mean\\tpool_sd\\tpool_CV\\" +
                "tpoolCV_over_sampleCV\\tchebi.id\\tchemspider.id\\tbiodb.compound.name");
    }

    public void createDummyFileWithOnlyColumn(String inputLine) {
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
        String[] expectedLine = {"nameMetabolite", "idSBML", "inchi", "chebi", "smiles", "pubchem", "inchikeys", "kegg", "hmd", "chemspider", "weight"};
        assertEquals(Arrays.toString(expectedLine), Arrays.toString((this.fingerprint.getEntityList()).iterator().next()));
    }

    public void testSeparator() {
        //Test that each possible mapping values are correctly extracted
        this.setMappingAllColumn();
        this.noFormatCheck = true;
        this.separator = ";";
        this.createDummyFileWithMultipleColumns("nameMetabolite;idSBML;inchi;chebi;smiles;pubchem;inchikeys;kegg;hmd;chemspider;weight");
        String[] expectedLine = {"nameMetabolite", "idSBML", "inchi", "chebi", "smiles", "pubchem", "inchikeys", "kegg", "hmd", "chemspider", "weight"};
        assertEquals(Arrays.toString(expectedLine), Arrays.toString((this.fingerprint.getEntityList()).iterator().next()));
    }

    public void testHeader() {
        //Test that each possible mapping values are correctly extracted
        this.ifNoHeader = true;
        this.createDummyFileWithOnlyColumn("M_taur");
        assertTrue(this.fingerprint.getEntityList().size() == 2);
    }

    public void testFiltered() {
        //Test that line(s) with empty values in the filtered column are discarded from the parsing
        this.filteredColumn = 24;
        this.createDummyFileWithMultipleColumns("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\tplsda|randomforest|svm\n" +
                "Pantothenic acid\tCHEBI:7916\tC9H17NO5\tCC(C)(CO)C(O)C(=O)NCCC(O)=O\tInChI=1S/C9H17NO5/c1-9(2,5-11)7(14)8(15)10-4-3-6(12)13/h7,11,14H,3-5H2,1-2H3,(H,10,15)(H,12,13)\tPantothenic acid\t218,102478\t1,00727647\t219,10975447\tNA\t[(M-H)]-\t1\t4,77\t5\t3,5599610222\t0,2982536819\t0,0837800414\t3012837,77207209\t131428,160471926\t0,043622714\t0,5206814567\t7916\tNA\tpantothenic acid");
        assertTrue(this.fingerprint.getEntityList().size() == 1);

        //Without the filtering, test that all the lines have been extracted
        this.filteredColumn = -1;
        this.createDummyFileWithMultipleColumns("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\tplsda|randomforest|svm\n" +
                "Pantothenic acid\tCHEBI:7916\tC9H17NO5\tCC(C)(CO)C(O)C(=O)NCCC(O)=O\tInChI=1S/C9H17NO5/c1-9(2,5-11)7(14)8(15)10-4-3-6(12)13/h7,11,14H,3-5H2,1-2H3,(H,10,15)(H,12,13)\tPantothenic acid\t218,102478\t1,00727647\t219,10975447\tNA\t[(M-H)]-\t1\t4,77\t5\t3,5599610222\t0,2982536819\t0,0837800414\t3012837,77207209\t131428,160471926\t0,043622714\t0,5206814567\t7916\tNA\tpantothenic acid");
        assertTrue(this.fingerprint.getEntityList().size() == 2);
    }

    public void testCheckingFormat() {
        this.checkingFile="temp/format_checking.tsv";
        this.layerWarning = true;
        this.setInChILayers("c", "h", "t");
        this.createDummyFileWithMultipleColumns("randomEntity\t\t\t\tInChI=1S/C25H36O8/c1-24-9-7-13(26)1/h11,14-21,23");
        this.setBufferTest(this.checkingFile,
                "Line number\tBioentity name\tDatabase\tInChI layer\tWrong value",
                "1\trandomEntity\tInChI\ttetraStereo (t)\tInChI=1S/C25H36O8/c1-24-9-7-13(26)1/h11,14-21,23");
    }

}
