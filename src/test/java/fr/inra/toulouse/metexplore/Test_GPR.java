/*******************************************************************************
 * Copyright INRA
 *
 *  Contact: ludovic.cottret@toulouse.inra.fr
 *
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *  In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *  The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *******************************************************************************/

package fr.inra.toulouse.metexplore;

public class Test_GPR extends Test {

    /***Mapping***/

    public void testMappingIDGene() {
        //network = new JSBML2Bionetwork("data/recon2.02.xml").getBioNetwork();

        this.bioEntityType = 6;
        this.setMapping4OneColumnFileByID("10026.1", "10026.1");
    }

    public void testMappingNameGene() {
        this.nameMapping = -1;
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


    /*** Pathway Enrichment ***/

    public void testWriteOutputPathEnrWithGene() {
        //Test the expected format of the output file obtained by pathway enrichment and with a reaction
        testMappingNameGene();
        this.setWriteOutputPathEnr(
                this.outputFile,
                "genes",
                "Pathway",
                "",
                "Phosphatidylinositol phosphate metabolism\t1.92\t1\t0.02823018458197611\t0.02823018458197611\t0.02823018458197611\t10026.1\t10026.1\t10026.1");
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
                "\tNb. of unmapped (pathway)\tNb. of unmapped (fingerprint)\tNb. of remaining (network)",
                "Phosphatidylinositol phosphate metabolism\t1.92\t1\t0.02823018458197611\t0.02823018458197611\t0.02823018458197611\t10026.1\t10026.1\t10026.1\t51\t0\t1790");
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
                "Fatty acid synthesis\t1.69\t1\t0.03203040173724213\t0.03203040173724213\t0.03203040173724213\thsa:9415 (TH)\t_9415_1_c\t_HSA:9415");
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
                "Fatty acid synthesis\t1.69\t1\t0.02199850857568978\t0.02199850857568978\t0.02199850857568978\thsa:9415 (TH)\t_9415_1_c\t_HSA:9415");
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
                "hsa:9415\t100.0\t1\t0.0005428881650380022\t0.0005428881650380022\t0.0005428881650380022\thsa:9415 (TH)\t_9415_1_c\t_HSA:9415");

    }

}
