package fr.inra.toulouse.metexplore;

public class Test_GPR extends Test {

    /***Mapping***/

    public void testMappingIDGene() {
        //network = new JSBML2Bionetwork("data/recon2.02.xml").getBioNetwork();

        this.bioEntityType = 6;
        this.setMapping4OneColumnFileByID("10026.1", "10026.1");
    }

    public void testMappingNameGene() {
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

    /**
     * Pathway Enrichment
     **/

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
        this.outputFile = "";
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
        this.outputFile = "";
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
        this.outputFile = "";
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
        this.entityType2Enrich = 5;
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
        this.entityType2Enrich = 6;
        this.setWriteOutputPathEnr(
                this.outputFile,
                "proteins",
                "Gene",
                "",
                "hsa:9415\t0.0005428881650380022\t0.0005428881650380022\t0.0005428881650380022\thsa:9415 (TH)\t_9415_1_c\t_HSA:9415\t1\t100.0");

    }

}
