package fr.inra.toulouse.metexplore;

public class Test_Mapping extends Test{

    /***InChI***/

    public void testMappingInChI () {
        //Test the success of a mapping with on first two layers (c and h) of InChI from a dataset containing multiple columns

        this.setMapping4MultipleColumnFile("Taurine\tNA\tC2H7NO3S\tC(CS(O)(=O)=O)N\tInChI=1S/C2H7NO3S/c3-1-2-7(4,5)6/h1-3H2,(H,4,5,6)" +
                        "\tTaurine\t124,006693\t1,00727647\t125,01396947\tNA\t[(M-H)]-\t1\t0,88\t5\t2,6122895216\t0,6358457794\t0,2434055545" +
                        "\t387859,346882448\t11652,3712684191\t0,0300427755\t0,1234268278\t15891\tNA\ttaurine",
                "M_taur");
    }

    public void testMappingLayerP () {
        //Test the success of a mapping with an extra p layer on InChI String

        //Positive test
        this.setInChILayers("c", "h", "p");
        this.setMapping4MultipleColumnFile("Glyceric acid\tCHEBI:33508\tC3H6O4\tOCC(O)C(O)=O\t" +
                        "InChI=1S/C3H6O4/c4-1-2(5)3(6)7/h2,4-5H,1H2,(H,6,7)/p-1\tGlyceric acid\t105,018901\tNA\t[(M-H)]-" +
                        "\t1\t0,92\t5\t2,9214987592\t0,2421744543\t0,0828939097\t646166,879585113\t24995,1780580671" +
                        "\t0,0386822334\t0,4666474721\t0,0574485787\t-0,0609303111\t-0,1375521553\t0,0010775982\t1\t-0,1318119477" +
                        "\t-0,0796214185\t1,3556765679\t0,6630650171\t-0,053608993",
                "M_glyc_R");

        //Negative test
//        this.setDefaultInChILayers();
//        try {
//            this.setMapping4MultipleColumnFile("Glyceric acid\tCHEBI:33508\tC3H6O4\tOCC(O)C(O)=O\t
// InChI=1S/C3H6O4/c4-1-2(5)3(6)7/h2,4-5H,1H2,(H,6,7)/p-1\tGlyceric acid\t105,018901\tNA\t[(M-H)]-\t1\t0,92\t5\t2,9214987592
// \t0,2421744543\t0,0828939097\t646166,879585113\t24995,1780580671\t0,0386822334\t0,4666474721\t0,0574485787\t-0,0609303111
// \t-0,1375521553\t0,0010775982\t1\t-0,1318119477\t-0,0796214185\t1,3556765679\t0,6630650171\t-0,053608993",
//                    "M_glyc_R");
//        }catch (ThreadDeath e){}//Multiple match are expected
    }

    public void testMappingLayerFormula () {
        //Test the success of a mapping with formula layer only on InChI String

        //Positive test
        this.setInChILayers("","","");
        this.setMapping4MultipleColumnFile("Cinnamoylglycine\tCHEBI:68616\tC11H11NO3\tOC(=O)CNC(=O)C=Cc1ccccc1\t" +
                        "InChI=1S/C11H11NO3/c13-10(12-8-11(14)15)7-6-9-4-2-1-3-5-9/h1-7H,8H2,(H,12,13)(H,14,15)/b7-6+\t" +
                        "Cinnamoylglycine\t204,065452\t1,00727647\t205,07272847\tNA\t[(M-H)]-\t1\t7,03\t5\t4,0160399219" +
                        "\t0,5133270871\t0,1278192192\t11742041,4996239\t2134365,54261586\t0,1817712484\t1,4220963763\t68616" +
                        "\tNA\tN-cinnamoylglycine",
                "M_5moxact");

        //Negative test
//        this.setDefaultInChILayers();
//        try{
//            this.setMapping4MultipleColumnFile("Cinnamoylglycine\tCHEBI:68616\tC11H11NO3\tOC(=O)CNC(=O)C=Cc1ccccc1\tInChI=1S/C11H11NO3/c13-10(12-8-11(14)15)7-6-9-4-2-1-3-5-9/h1-7H,8H2,(H,12,13)(H,14,15)/b7-6+\tCinnamoylglycine\t204,065452\t1,00727647\t205,07272847\tNA\t[(M-H)]-\t1\t7,03\t5\t4,0160399219\t0,5133270871\t0,1278192192\t11742041,4996239\t2134365,54261586\t0,1817712484\t1,4220963763\t68616\tNA\tN-cinnamoylglycine",
//                    "M_5moxact");
//        }catch (ThreadDeath  e){}//No match is expected
    }

    /***Others***/

    public void testMappingCHEBI () {
        //Test the success of a mapping with CHEBI
        this.setMappingColumn(2,1);
        this.setMapping4Testo();
    }

    public void testMappingHMDB () {
        //Test the success of a mapping with the HMDB ID
        this.setMapping4OneColumnFile(7,0,"HMDB00251", "M_taur");
    }

    public void testMappingInchiKey () {
        //Test the success of a mapping with the Inchikey
        this.setMapping4OneColumnFile(5,0,"XOAAWQZATWQOTB-UHFFFAOYSA-N", "M_taur");
    }

    public void testMappingMass () {
        //Test the success of a mapping with the Inchikey
        this.weightPrecision = 2;
        this.setMapping4OneColumnFile(9,0,"125.01453", "M_taur");
    }

    public void testMappingPubChem () {
        //Test the success of a mapping with the Inchikey
        this.setMapping4OneColumnFile(4,0,"4068592", "M_taur");
    }

    public void testMappingKegg () {
        //Test the success of a mapping with the KEGG ID
        this.setMapping4OneColumnFile(6,0,"C00160", "M_glyclt");
    }

    public void testMappingID () {
        //Test the success of a mapping with the ID of the network from a dataset containing an only column
        this.setMapping4OneColumnFileByID("M_taur", "M_taur");
    }

    public void testMappingName () {
        //Test the success of a mapping with the name of the network from a dataset containing an only column
        this.nameMapping = 1;
        this.setMapping4OneColumnFileByID("Taurine", "M_taur");
    }

    public void testMappingIDReaction () {
        //Test the success of a mapping with the ID of a reaction
        this.bioEntityType = 2;
        this.setMapping4OneColumnFileByID("R_FUM", "R_FUM");
    }

    public void testMappingNameReaction () {
        //Test the success of a mapping with the ID of the network from a dataset containing an only column
        this.nameMapping = 1;
        this.bioEntityType = 2;
        this.setMapping4OneColumnFileByID("fumarase", "R_FUM");
    }

    public void testMappingIDPathway() {
        //Test the success of a mapping with the ID of a pathway
        this.bioEntityType = 3;
        this.setMapping4OneColumnFileByID("Fatty acid oxidation", "Fatty acid oxidation");
    }

    public void testMappingNamePathway () {
        //Test the success of a mapping with the name of a pathway
        this.nameMapping = 1;
        this.bioEntityType = 3;
        this.setMapping4OneColumnFileByID("Fatty acid oxidation", "Fatty acid oxidation");
    }

    /*************Writing output***********/
    public void testWriteOutputMapping(){
        //Test the expected format of the output file obtained by mapping
        this.setMapping4Testo();
        this.setBufferTest(this.outputFile,
                "Mapped\tName_(Input_File)\tName_(SBML)\tSBML_ID\tMatched_value_(Input_File)\tMatched_value_(SBML)",
                "true\tTestosterone glucuronide\ttestosterone 3-glucosiduronic acid\tM_tststeroneglc" +
                        "\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1" +
                        "\tInChI=1S/C25H36O8/c1-24-9-7-13(26)11-12(24)3-4-14-15-5-6-17(25(15,2)10-8-16(14)24)32-23-20(29)18(27)19(28)21(33-23)22(30)31/h11,14-21,23,27-29H,3-10H2,1-2H3,(H,30,31)/t14-,15-,16-,17-,18-,19-,20+,21-,23+,24-,25-/m0/s1");
    }

}
