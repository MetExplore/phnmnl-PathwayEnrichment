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

public class Test_PathEnr extends Test{


    public void testWriteOutputPathEnr() {
        this.setMapping4Testo();
        this.setWriteOutputPathEnr(
                this.outputFile,
                "metabolites",
                "Pathway",
                "",
                "Steroid metabolism\t1.67\t1\t0.02314814814814815\t0.02314814814814815\t0.02314814814814815\t" +
                        "testosterone 3-glucosiduronic acid\tTestosterone glucuronide\tM_tststeroneglc");
    }

    public void testWriteOutputReacEnr() {
        this.setMapping4Testo();
        this.entityType2Enrich=2;
        this.setWriteOutputPathEnr(
                this.outputFile,
                "metabolites",
                "Reaction",
                "",
                "UDP-glucuronosyltransferase 1-10 precursor, microsomal\t25.0\t1\t0.00154320987654321\t0.00154320987654321\t0.00154320987654321\ttestosterone 3-glucosiduronic acid\tTestosterone glucuronide\tM_tststeroneglc");
    }

    public void testWriteOutputPathEnr4Galaxy() {
        this.setGalaxy();
        this.setMapping4Testo();
        this.setWriteOutputPathEnr(
                this.outputFile,
                "metabolites",
                "Pathway",
                "\tNb. of unmapped (pathway)\tNb. of unmapped (fingerprint)\tNb. of remaining (network)",
                "Steroid metabolism\t1.67\t1\t0.02314814814814815\t0.02314814814814815\t0.02314814814814815\ttestosterone 3-glucosiduronic acid\tTestosterone glucuronide\tM_tststeroneglc\t59\t0\t2532");
        this.setBufferTest(this.galaxy,
                "1 metabolite has been mapped on 1 in the fingerprint dataset (100.0%) and on 2592 in the network (0.04%).",
                "1 pathway is concerned among the network (on 97 in the network; 1.03%).");
    }

    public void testWriteOutputPathEnrWithReaction() {
        //Test the expected format of the output file obtained by pathway enrichment and with a reaction
        this.bioEntityType = 2;
        this.nameMapping = 1;
        this.setMapping4OneColumnFileByID("fumarase", "R_FUM");
        this.setWriteOutputPathEnr(
                this.outputFile,
                "reactions",
                "Pathway",
                "",
                "Citric acid cycle\t5.0\t1\t0.004750593824228029\t0.004750593824228029\t0.004750593824228029\tfumarase\tfumarase\tR_FUM");
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
                "\tNb. of unmapped (pathway)\tNb. of unmapped (fingerprint)\tNb. of remaining (network)",
                "Citric acid cycle\t5.0\t1\t0.004750593824228029\t0.004750593824228029\t0.004750593824228029\tfumarase\tR_FUM\tR_FUM\t19\t0\t4190");
        this.setBufferTest(this.galaxy,
                "1 pathway is concerned among the network (on 97 in the network; 1.03%).",
                null);
    }

}
