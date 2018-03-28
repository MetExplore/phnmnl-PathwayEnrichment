package fr.inra.toulouse.metexplore;

import fr.inra.toulouse.metexplore.io.Fingerprint;
import fr.inra.toulouse.metexplore.io.JSBML2Bionetwork;
import fr.inra.toulouse.metexplore.io.WritingComportment;
import fr.inra.toulouse.metexplore.omics.Mapping;
import fr.inra.toulouse.metexplore.omics.PathwayEnrichment;
import junit.framework.TestCase;
import org.apache.commons.io.FileUtils;
import parsebionet.biodata.BioEntity;
import parsebionet.biodata.BioNetwork;

import java.io.*;

import java.security.Permission;
import java.util.ArrayList;
import java.util.List;

public class Test extends TestCase implements WritingComportment {
    protected String separator, outputFile, galaxy, checkingFile, logContent, logFile = "temp/information.txt", dummyFile = "temp/dummy.tsv";
    protected int nameMapping, filteredColumn, bioEntityType, entityType2Enrich, weightPrecision;
    protected Boolean ifNoHeader, noFormatCheck, layerWarning;
    protected String[] inchiLayers;
    protected int[] mappingColumn;
    protected List<BioEntity> expectedMappedMetabolite;
    // ID biosource:3223
    protected static BioNetwork network;
    protected Fingerprint fingerprint;
    protected File file;
    protected BufferedReader buffer;
    protected Mapping mapping;
    protected PathwayEnrichment pathEnr;
    protected static int nbTestsRecFlat = 0, nbTestsRecFull = 0;


    //Instantiation mode to deals with appropriate SBML and both Maven and hand by hand tests running modes
    public Test(){
        String classe = this.getClass().getSimpleName();
        if(classe.equals("Test_Mapping") || classe.equals("Test_PathEnr")){
            nbTestsRecFlat++;
            nbTestsRecFull = 0;
        }else if(classe.equals("Test_GPR")){
            nbTestsRecFlat = 0;
            nbTestsRecFull++;
        }
        if (nbTestsRecFull == 1){
            network = new JSBML2Bionetwork("data/recon2.02.xml").getBioNetwork();
        }else if(nbTestsRecFlat == 1){
            network = (new JSBML2Bionetwork("data/recon2.02_without_compartment.xml")).getBioNetwork();
        }
    }


    public void setUp() throws Exception {
        //Initialization of the parameters before each tests
        super.setUp();
        //this.setSecurityManager();
        //WritingComportment write = new WritingComportment(this.galaxy);
        (new File("temp")).mkdir();
        this.setDefaultInChILayers();
        this.setDefaultMappingColumn();
        this.noFormatCheck = false;
        this.nameMapping = -1;
        this.fingerprint = null;
        this.mapping = null;
        this.pathEnr = null;
        this.separator = "\t";
        this.ifNoHeader = false;
        this.layerWarning = false;
        this.galaxy = "";
        this.outputFile = "temp/output.tsv";
        this.checkingFile = "";
        this.logContent = "";
        this.filteredColumn = -1;
        this.bioEntityType = 1;
        this.entityType2Enrich = 3;
        this.weightPrecision = 2;
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

    public static void setSecurityManager() {
        //Allow to catch exit(1) when object is called in a main function (or by other object)
        //but not when running this class as testunit
        System.setSecurityManager(new SecurityManager() {

            @Override
            public void checkPermission(Permission perm) {
            }

            @Override
            public void checkExit(int status) {
                throw new ThreadDeath();
            }

        });
    }

    /************setWriteOutput************/
    public void setBufferReader(String outFile) {
        this.file = new File(outFile);
        try {
            this.buffer = new BufferedReader(new FileReader(this.file));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void setBufferTest(String fileName, String header, String line) {
        setBufferReader(fileName);
        try {
            if (fileName.equals(this.logFile)) {
                this.logContent = this.pathEnr.getLogContent();
                writeOutput(this.logContent,this.file);
            }
            assertEquals(buffer.readLine(), header);
            assertEquals(buffer.readLine(), line);
            assertEquals(buffer.readLine(), null);
        } catch (IOException e) {
            e.printStackTrace();
        }
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
            this.fingerprint = new Fingerprint(this.logContent,this.layerWarning, this.noFormatCheck, this.dummyFile, this.ifNoHeader, this.separator, ";", 0,
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

    /************setMapping************/

    public void setMapping(String bpe) {
        //Extract the expected metabolite and process a mapping for comparison

        this.expectedMappedMetabolite = new ArrayList<>();

        try {
            this.mapping = new Mapping(this.logContent,this.network, this.fingerprint.getEntityList(),
                    this.inchiLayers, this.nameMapping, this.weightPrecision, this.outputFile,
                    this.galaxy, this.bioEntityType);
            this.logContent = mapping.getLogContent();
            this.file = new File(this.outputFile);
            this.expectedMappedMetabolite.add((BioEntity) mapping.getEntitySetInNetwork().get(bpe));
            assertEquals(this.expectedMappedMetabolite.iterator().next().getName(),
                    this.mapping.getList_mappedEntities().keySet().iterator().next().getName());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void setMapping4OneColumnFileByID(String inputLine, String bpe) {
        this.setMapping4OneColumnFile(0, 0, inputLine, bpe);
    }

    public void setMapping4OneColumnFile(int mappingType, int mappingColumn, String inputLine, String bpe) {
        this.setMappingColumn(mappingType, mappingColumn);
        this.createDummyFileWithOnlyColumn(inputLine);
        this.setMapping(bpe);
    }

    public void setMapping4MultipleColumnFile(String inputLine, String bpe) {
        this.createDummyFileWithMultipleColumns(inputLine);
        this.setMapping(bpe);
    }

    public void setMapping4Testo() {
        this.setMapping4MultipleColumnFile("Testosterone glucuronide\tCHEBI:28835\tC25H36O8\t" +
                        "[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O" +
                        "\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1" +
                        "\tTestosterone glucuronide\t463,2329\t1,00727647\t464,24017647\tNA\t[(M-H)]-\t1\t7,9\t4\t2,1475578771\t0,5701078279\t0,265467969" +
                        "\t178149,617939526\t12351,5841321731\t0,0693326445\t0,2611714128\t28835\tNA\ttestosterone 17-glucosiduronic acid\n",
                "M_tststeroneglc");
    }

    /*************Writing output***********/
    public void setWriteOutputPathEnr(String pathFile, String typeOfMappedEntity,String typeOfEnrichedEntity, String galColumn, String line) {
        //Test the expected format of the output file obtained by pathway enrichment
        try {
            this.pathEnr = new PathwayEnrichment(this.logContent,network,this.fingerprint.getEntityList(),
                    this.mapping.getList_mappedEntities(), pathFile,this.galaxy,this.bioEntityType, this.entityType2Enrich);
        }catch (IOException e ){
            e.printStackTrace();
        }
        this.setBufferTest(
                pathFile,
                typeOfEnrichedEntity + " name\tCoverage (%)\tNb. of mapped\tP-value\tBonferroni corrected p-value\tBH corrected p-value\t" +
                        "Mapped " + typeOfMappedEntity + " (SBML)\tMapped " + typeOfMappedEntity + " (fingerprint)\tMapped "
                        + typeOfMappedEntity + " ID" +galColumn,
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
    }

    /*******************************
     *              TEST
     *****************************/

    public void testToPreventMavenErrorDuringBuild(){
        assertTrue(true);
    }
}
