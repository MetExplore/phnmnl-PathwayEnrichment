package fr.inra.toulouse.metexplore;

public class Test_PathEnr extends Test{


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

}
