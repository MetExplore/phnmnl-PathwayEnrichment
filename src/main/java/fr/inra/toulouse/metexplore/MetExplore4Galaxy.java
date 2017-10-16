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

    public HashMap <String, String[]>  extractData(String inputFile, Integer filteredColumns) throws IOException {

        Boolean filtered = true;

        if(filteredColumns < 1 ){
            filtered = false;
        }

        HashMap <String, String[]> listMetabolites = new HashMap<String, String[]>();
        BufferedReader f = new BufferedReader(new FileReader(new File(inputFile)));
        String line;
        int i = 0;

        while ((line = f.readLine()) != null) {
            if (i != 0) {
                String[] values = line.replaceAll("\"","").split("\t");
                try {
                    if (filtered ==  false || values[filteredColumns] != ""){
                        listMetabolites.put(String.valueOf(i), values);
                    }
                } catch (ArrayIndexOutOfBoundsException e) {
                }
            }
            ++i;
        }
        if (f != null) {
            f.close();
        }
        if (listMetabolites.size() < 1 ){
            System.err.println("File badly formatted");
            exit(1);
        }
       return listMetabolites;
    }

    public Set<BioPhysicalEntity> mapping(BioNetwork bn, String outputFile1, HashMap <String, String[]> parsedFile, int chebiColumn, int inchiColumn, String[] inchiLayers){

        Set<BioPhysicalEntity> listMetabolites = new HashSet();
        int i = 1;

        try {
            File fo1 = new File(outputFile1);
            fo1.createNewFile();
            BufferedWriter f = new BufferedWriter(new FileWriter(fo1));


            HashMap<String, String[]> remainingMetabolites = (HashMap<String, String[]>) parsedFile.clone();
            Boolean inchiMapping;

            f.write("Mapped\tInputFile's name\tMetExplore's name\n");

            for (String[] entry : parsedFile.values()) {
                for (BioPhysicalEntity bpe : bn.getPhysicalEntityList().values()) {
                    inchiMapping = false;
                    if (!(entry[inchiColumn]).equals("NA") && !(entry[inchiColumn]).equals("")
                            && (new InChI4Galaxy(bpe.getInchi(), inchiLayers)).equals(new InChI4Galaxy(entry[inchiColumn], inchiLayers))
                            && (inchiColumn > 0) ) {
                        listMetabolites.add(bpe);
                        remainingMetabolites.remove(String.valueOf(i));
                        f.write("true\t" + entry[0] + "\t" + bpe.getName() + "\n");
                        inchiMapping = true;
                    }

                    if (!inchiMapping) {

                        if (!(entry[chebiColumn]).equals("NA") && !((entry[chebiColumn]).equals("")) && (chebiColumn > 0) ) {

                            for (Map.Entry<String, Set<BioRef>> ref : bpe.getRefs().entrySet()) {
                                if (ref.getKey().equals("chebi")) {
                                    for (BioRef val : ref.getValue()) {
                                        if (entry[chebiColumn].equals(val.id)) {
                                            listMetabolites.add(bpe);
                                            remainingMetabolites.remove(String.valueOf(i));
                                            f.write( "true\t" + entry[0] + "\t" + bpe.getName() + "\n");
                                            break;
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
                ++i;
            }
            System.err.println((parsedFile.size() - remainingMetabolites.size()) + " metabolites have been mapped (on " + parsedFile.size() + ").\n");
            for (String[] entry : remainingMetabolites.values()){
                f.write("false\t" + entry[0] +"\n");
            }
            f.close();
        }catch (IOException e){
            e.printStackTrace();
        }
            return listMetabolites;
    }

    public ArrayList<HashMap<BioPathway, Double>> pathwayEnrichment (BioNetwork bn, Set<BioPhysicalEntity> map, int correction) {

        ArrayList<HashMap<BioPathway, Double>> resultList = new ArrayList<HashMap<BioPathway, Double>>();
        System.err.println("Pathway enrichment in progress...");
        PathwayEnrichment enr = new PathwayEnrichment(bn, map);

        resultList.add(enr.computeEnrichment());
        resultList.add(enr.computeEnrichment(correction));

        System.err.println(resultList.get(0).size() + " pathways are concerned among the network (on " + bn.getPathwayList().size() + ").");

        return resultList;
    }

    public void writeOutput(ArrayList<HashMap<BioPathway, Double>> resultList, Set<BioPhysicalEntity> map, String outputFile2) throws IOException{

        File fo2 = new File(outputFile2);
        fo2.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo2));

        f.write("Pathway enrichment\tFisher's p-value\tTest correction\tMapped metabolites\tNb of mapped\tCoverage (%)\n");
        HashMap<BioPathway, Double> result = resultList.get(0);
        Iterator itCorr = resultList.get(1).values().iterator();

        for (Map.Entry<BioPathway, Double> entry : result.entrySet()) {
            BioPathway path = entry.getKey();
            String printed = entry.getKey().getName() + "\t" + roundPval((double)entry.getValue()) + "\t";
            printed += roundPval((double) itCorr.next()) + "\t";

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

    public String roundPval(double value) {
        if(value < 0.01){
            String tmp = (String.valueOf(value));
            String[] splitted = tmp.split("E");
            double head = Double.parseDouble(splitted[0]);

            try{
                return  (String.valueOf(round(head)) + " E" + splitted[1] +" ");
            }catch (ArrayIndexOutOfBoundsException e){
                splitted = tmp.split("\\.");
                int j = 0;
                for (int i = 0 ; i < splitted[1].length(); i++){
                    if(splitted[1].charAt(i) == '0'){
                        j++;
                    }else{
                        break;
                    }
                }
                double tail = Double.parseDouble(splitted[1].substring(j, j+1) + "." + splitted[1].substring(j + 1));
                return (round(tail) + " E-" + String.valueOf(j + 1));
            }
        }
        return round(value);
    }

    public String round(double value) {
        return String.valueOf((double) Math.round(value * 100) / 100);
    }

    public void runTime(long elapsedTime){
        if (elapsedTime / 60000000000L == 0L) {
            System.err.println("Time to run the process in seconde :" + elapsedTime / 1000000000L);
        } else {
            long min = elapsedTime / 60000000000L;
            long sec = elapsedTime / 1000000000L - (min * 60);
            System.err.println("Time to run the process : " + min + "min " + sec + "s");
        }
    }

    public static void main( String[] args ){

        //Parameters
        long startTime = System.nanoTime();
        String inputFile = "data/sacurineVariableMetadataEnhanced.tsv";
        String sbml = "data/recon2.v03_ext_noCompartment_noTransport_v2.xml";
        String outputFile1 = "mapping.tsv";
        String outputFile2 = "pathwayEnrichment.tsv";
        String[] inchiLayers = {"c","h"};
        int correction = 1;

        //SBML parsing
        BioNetwork bionet = (new JSBMLToBionetwork(sbml)).getBioNetwork();
        bionet.printBioNetworkSizeToOut();

        MetExplore4Galaxy app = new MetExplore4Galaxy();
        try{

            //Enrichment
            HashMap <String, String[]> parsedFile = app.extractData(inputFile, -1);
            Set<BioPhysicalEntity> map = app.mapping(bionet, outputFile1, parsedFile, 1, 4, inchiLayers);
            app.writeOutput(app.pathwayEnrichment (bionet, map, correction), map, outputFile2);
       }
        catch (IOException e){
            e.printStackTrace();
        }
        app.runTime(System.nanoTime() - startTime);
    }
}