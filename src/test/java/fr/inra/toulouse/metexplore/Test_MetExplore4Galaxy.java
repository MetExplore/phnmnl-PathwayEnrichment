package fr.inra.toulouse.metexplore;

import junit.framework.TestCase;
import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPhysicalEntity;
import parsebionet.io.JSBMLToBionetwork;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.regex.Pattern;

public class Test_MetExplore4Galaxy extends TestCase {

    static MetExplore4Galaxy met = null;
    static BioNetwork bionet = null;
    static HashMap<String, String[]> parsedFile = new HashMap <String, String[]> ();
    static String inputFile;


    protected void setUp() throws Exception {
        super.setUp();
        String dir = "/home/bmerlet/Documents/PathwayEnrichment/";
        //inputFile = dir + "Galaxy15-[Biosigner_Multivariate_Univariate_Multivariate_variableMetadata.tsv].tabular";
        inputFile = "sacurineVariableMetadataEnhanced.tsv";
        String sbml = dir + "recon2.v03_ext_noCompartment_noTransport.xml";
        String outputFile = dir + "output.tsv";

        if (met == null) {
            try {

                bionet = (new JSBMLToBionetwork(sbml)).getBioNetwork();
                met = new MetExplore4Galaxy();
                parsedFile = met.extractData(inputFile, false, -1);

            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }




}
