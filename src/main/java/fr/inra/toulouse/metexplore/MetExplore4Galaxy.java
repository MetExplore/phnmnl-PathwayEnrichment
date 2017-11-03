package fr.inra.toulouse.metexplore;

import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPathway;
import parsebionet.biodata.BioPhysicalEntity;
import parsebionet.biodata.BioRef;
import parsebionet.io.JSBMLToBionetwork;
import parsebionet.statistics.PathwayEnrichment;

import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.File;

import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.Map;
import java.util.ArrayList;
import java.util.Iterator;

import static java.lang.System.exit;

public class MetExplore4Galaxy {

    private BufferedWriter f3;

    public MetExplore4Galaxy(String outputFile3){
       try {
           if (outputFile3 != "") {
               File fo3 = new File(outputFile3);
               fo3.createNewFile();
               this.f3 = new BufferedWriter(new FileWriter(fo3));
           }
       }catch (IOException e){
        e.printStackTrace();
        }
    }

    public HashMap <String, String[]>  extractData(String inputFile, Integer filteredColumns) throws IOException {

        Boolean filtered = true;

        if(filteredColumns < 1 ) filtered = false;

        HashMap <String, String[]> listMetabolites = new HashMap<String, String[]>();
        BufferedReader f = new BufferedReader(new FileReader(new File(inputFile)));
        String line;
        int i = 0; int j = 0;

        while ((line = f.readLine()) != null) {
            if (i != 0) {
                String[] values = line.replaceAll("\"","").split("\t");
                try {
                    if (filtered ==  false || values[filteredColumns] != ""){
                        listMetabolites.put(String.valueOf(j), values);
                        ++j;
                    }
                } catch (ArrayIndexOutOfBoundsException e) {
                }
            }
            ++i;
        }
        if (f != null) f.close();
        if (listMetabolites.size() < 1 ){
            System.err.println("File badly formatted");
            exit(1);
        }
       return listMetabolites;
    }

    public Set<BioPhysicalEntity> mapping(BioNetwork bn, String outputFile1, String outputFile3, HashMap <String, String[]> parsedFile, int chebiColumn, int inchiColumn, String[] inchiLayers){

        Set<BioPhysicalEntity> listMetabolites = new HashSet();
        int i = 0;

        try {
            File fo1 = new File(outputFile1);
            fo1.createNewFile();
            BufferedWriter f = new BufferedWriter(new FileWriter(fo1));

            HashMap<String, String[]> remainingMetabolites = (HashMap<String, String[]>) parsedFile.clone();
            Boolean mapping = false; Boolean mappingBpe = false;

            f.write("Mapped\tInputFile's name\tSBML's name\tInfile's val\tSBML's val\n");

            for (String[] entry : parsedFile.values()) {
                mapping = false;
                for (BioPhysicalEntity bpe : bn.getPhysicalEntityList().values()) {
                    mappingBpe = false;
                    if ((inchiColumn > 0) && !(entry[inchiColumn]).equals("NA") && !(entry[inchiColumn]).equals("")
                            && (new InChI4Galaxy(bpe.getInchi(), inchiLayers)).equals(new InChI4Galaxy(entry[inchiColumn], inchiLayers))) {
                        listMetabolites.add(bpe);
                        f.write("true\t" + entry[0] + "\t" + bpe.getName() + "\t" + entry[inchiColumn] + "\t" + bpe.getInchi() + "\n");
                        mappingBpe = true; mapping= true;
                    }

                    if (!mappingBpe) {

                        if ((chebiColumn > 0) && !(entry[chebiColumn]).equals("NA") && !((entry[chebiColumn]).equals(""))) {

                            for (Map.Entry<String, Set<BioRef>> ref : bpe.getRefs().entrySet()) {
                                if (ref.getKey().equals("chebi")) {
                                    for (BioRef val : ref.getValue()) {
                                        if (entry[chebiColumn].equals(val.id)) {
                                            listMetabolites.add(bpe);
                                            f.write( "true\t" + entry[0] + "\t" + bpe.getName() + "\t" + entry[chebiColumn] + "\t" + val.id + "\n");
                                            mappingBpe = true; mapping=true;
                                            break;
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
                if (mapping) remainingMetabolites.remove(String.valueOf(i));
                ++i;
            }
            String printingMessage = (parsedFile.size() - remainingMetabolites.size()) + " metabolites have been mapped (on " + parsedFile.size() + ").";
            System.out.println(printingMessage);
            if (outputFile3 != "") {
                this.f3.write(printingMessage);
            }
            for (String[] entry : remainingMetabolites.values()){
                f.write("false\t" + entry[0] +"\n");
            }
            f.close();
        }catch (IOException e){
            e.printStackTrace();
        }
        return listMetabolites;
    }

    public ArrayList<HashMap<BioPathway, Double>> pathwayEnrichment (BioNetwork bn, String outputFile3, Set<BioPhysicalEntity> map) {

        ArrayList<HashMap<BioPathway, Double>> resultList = new ArrayList<HashMap<BioPathway, Double>>();
        System.out.println("Pathway enrichment in progress...");
        PathwayEnrichment enr = new PathwayEnrichment(bn, map);
        HashMap<BioPathway, Double> res = enr.computeEnrichment();

        resultList.add(res);
        resultList.add(enr.bonferroniCorrection(res));
        resultList.add(enr.benjaminiHochbergCorrection(res));

        String printingMessage = resultList.get(0).size() + " pathways are concerned among the network (on " + bn.getPathwayList().size() + ").";
        System.out.println(printingMessage);

        if (outputFile3 != "") {
            try{
                this.f3.write("\n" + printingMessage);
                this.f3.close();
            }catch (IOException e) {
                e.printStackTrace();
            }
        }
        return resultList;
    }

    public void writeOutput(ArrayList<HashMap<BioPathway, Double>> resultList, Set<BioPhysicalEntity> map, String outputFile2) throws IOException{

        File fo2 = new File(outputFile2);
        fo2.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo2));

        f.write("Pathway enrichment\tFisher's p-value\tBonferroni correction\tBenjamini-Hochberg correction\tMapped metabolites\tNb of mapped\tCoverage (%)\n");
        HashMap<BioPathway, Double> result = resultList.get(0);
        Iterator itBonCorr = resultList.get(1).values().iterator();
        Iterator itBenHocCorr = resultList.get(2).values().iterator();

        for (Map.Entry<BioPathway, Double> entry : result.entrySet()) {
            BioPathway path = entry.getKey();
            String printed = entry.getKey().getName() + "\t" + removeSciNot((double)entry.getValue()) + "\t";
            printed += removeSciNot((double) itBonCorr.next()) + "\t" + removeSciNot((double) itBenHocCorr.next()) + "\t";

            int j = 0;

            for (BioPhysicalEntity bpe : map) {
                if (path.getListOfInvolvedMetabolite().containsValue(bpe)) {
                   printed += (j == 0) ? bpe.getName() : "; " + bpe.getName();
                   j++;
                }
            }
            String coverage = round((double) j / (double) path.getListOfInvolvedMetabolite().size() * (double) 100);
            f.write(printed + "\t" + j + "\t" + coverage + "\n");
        }
        if (f != null) {
            f.close();
        }
    }

    public String removeSciNot(double value) {
        String tmp = (String.valueOf(value));
        if(value < 1e-10){
            String[] splitted = tmp.split("E");
            String power = splitted[1];
            if(power.charAt(1) == '0') {
                power = power.substring(2,power.length());
            }
            else power = power.substring(1,power.length());
            String[] number = splitted[0].split("\\.");
            return "0." + (new String(new char[Integer.parseInt(power)-1]).replace("\0", "0")) + number[0] + number[1];
        }
        return tmp;
    }

    public String round(double value) {
        return String.valueOf((double) Math.round(value * 100) / 100);
    }

    public void runTime(long elapsedTime){
        if (elapsedTime / 60000000000L == 0L) {
            System.out.println("Time to run the process in seconde :" + elapsedTime / 1000000000L);
        } else {
            long min = elapsedTime / 60000000000L;
            long sec = elapsedTime / 1000000000L - (min * 60);
            System.out.println("Time to run the process : " + min + "min " + sec + "s");
        }
    }

    public static void main( String[] args ){

        //Parameters
        long startTime = System.nanoTime();
        String inputFile = "data/sacurineVariableMetadataEnhanced.tsv";
        String sbml = "data/recon2.v03_ext_noCompartment_noTransport_v2.xml";
        String outputFile1 = "mapping.tsv";
        String outputFile2 = "pathwayEnrichment.tsv";
        String outputFile3 = "info.txt";
        String[] inchiLayers = {"c","h"};

        //SBML parsing
        BioNetwork bionet = (new JSBMLToBionetwork(sbml)).getBioNetwork();
        bionet.printBioNetworkSizeToOut();

        MetExplore4Galaxy app = new MetExplore4Galaxy(outputFile3);
        try{

            //Enrichment
            HashMap <String, String[]> parsedFile = app.extractData(inputFile, -1);
            Set<BioPhysicalEntity> map = app.mapping(bionet, outputFile1, outputFile3, parsedFile, 1, 4, inchiLayers);
            app.writeOutput(app.pathwayEnrichment (bionet, outputFile3, map), map, outputFile2);
       }
        catch (IOException e){
            e.printStackTrace();
        }
        app.runTime(System.nanoTime() - startTime);
    }
}