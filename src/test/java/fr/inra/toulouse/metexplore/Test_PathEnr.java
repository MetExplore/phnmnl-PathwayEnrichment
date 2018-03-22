package fr.inra.toulouse.metexplore;


import fr.inra.toulouse.metexplore.io.JSBML2Bionetwork;
import fr.inra.toulouse.metexplore.omics.PathwayEnrichment;

import java.io.File;
import java.io.IOException;

public class Test_PathEnr extends Test_Mapping{

    protected int entityType2Enrich;

    public void setUp() throws Exception {
        super.setUp();
        this.pathEnr = null;
        this.entityType2Enrich = 3;
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
        network = new JSBML2Bionetwork("data/recon2.02.xml").getBioNetwork();
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

    public void itestWriteOutputProtEnrWithGene() {
        //Test the expected format of the output file obtained by pathway enrichment and with a reaction
        testMappingNameGene();
        this.entityType2Enrich=5;
        this.setWriteOutputPathEnr(
                this.outputFile,
                "genes",
                "Protein",
                "",
                "10026.1 (TH)\t0.0005428881650380022\t0.002714440825190011\t0.002714440825190011\t10026.1\t10026.1\t10026.1\t1\t100.0");
    }

    public void testWriteOutputGeneEnrWithProt() {
        //Test the expected format of the output file obtained by pathway enrichment and with a reaction
        testMappingIDProtein();
        this.entityType2Enrich=6;
        this.setWriteOutputPathEnr(
                this.outputFile,
                "proteins",
                "Gene",
                "",
                "hsa:9415\t0.0005428881650380022\t0.0005428881650380022\t0.0005428881650380022\thsa:9415 (TH)\t_9415_1_c\t_HSA:9415\t1\t100.0");
    }
}
