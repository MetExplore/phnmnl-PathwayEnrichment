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

import java.util.Arrays;

public class Test_Fingerprint extends Test{

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
        this.checkingFile="temp/checking_format.tsv";
        this.noFormatCheck=false;
        this.layerWarning = true;
        this.setInChILayers("c", "h", "t");
        this.createDummyFileWithMultipleColumns("randomEntity\t\t\t\tInChI=1S/C25H36O8/c1-24-9-7-13(26)1/h11,14-21,23/t...");
    }

}
