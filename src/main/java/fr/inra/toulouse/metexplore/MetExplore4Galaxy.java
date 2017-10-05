package fr.inra.toulouse.metexplore;

import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPathway;
import parsebionet.biodata.BioPhysicalEntity;
import parsebionet.io.JSBMLToBionetwork;
import parsebionet.statistics.PathwayEnrichment;
import parsebionet.utils.chemicalStructures.InChI;

import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.File;

import java.util.*;

public class MetExplore4Galaxy {

    public HashMap <String, String[]>  extractData(String inputFile, Boolean filterbyColumns, Integer filteredColumns) throws IOException {

        HashMap <String, String[]> listMetabolites = new HashMap<String, String[]>();
        BufferedReader f = new BufferedReader(new FileReader(new File(inputFile)));
        String line;
        int i = 0;

        while ((line = f.readLine()) != null) {
            if (i != 0) {
                String[] values = line.replaceAll("\"","").split("\t");
                try {
                    if ((filterbyColumns == true && values[filteredColumns] != "") || filterbyColumns == false) {
                        listMetabolites.put(values[0], values);
                    }
                } catch (ArrayIndexOutOfBoundsException e) {
                }
            }
            ++i;
        }
        if (f != null) {
            f.close();
        }
        return listMetabolites;
    }

    public Set<BioPhysicalEntity> mapping(BioNetwork bn, HashMap <String, String[]> parsedFile, int inchiColumn, String[] inchiLayers){//TODO: mapping on more InChI layout

        Set<BioPhysicalEntity> listMetabolites = new HashSet();
        HashMap <String, String[]> remainingMetabolites = parsedFile;


        System.out.println("MetExplore's name\tInputFile's name");//TODO: write a second file with these mapping informations or put them into dictionary or an attribute of a BioPhysicalEntity redefinition for the main file
        for (BioPhysicalEntity bpe : bn.getPhysicalEntityList().values()) {//TODO-BUG: does a metabolite/entry match with only one bioentity? If no, otherwhise remove() is inadapted
            for (String[] entry : parsedFile.values()) {//TODO: convert into an arrayList and sort the list alphabetic
                //System.out.println(bpe.getInchi() + " = " + entry[4]);
                if (!(entry[inchiColumn]).equals("NA") && !(bpe.getInchi()).equals("NA") && !(bpe.getInchi()).equals("") && !(entry[inchiColumn]).equals("") ){
                    //InChI bpeInchi = new InChI(bpe.getInchi());
                    //InChI entryInchi = new InChI(entry[inchiColumn]);
                    InChI4Galaxy bpeInchi = new InChI4Galaxy(bpe.getInchi(), inchiLayers);
                    InChI4Galaxy entryInchi = new InChI4Galaxy(entry[inchiColumn], inchiLayers);
                    if (bpeInchi.equals(entryInchi)) {
                        listMetabolites.add(bpe);
                        //remainingMetabolites.remove(entry[0]);
                        System.out.println(bpe.getName() + " = " + entry[0]);
                        //break;
                    }
                }
            }
        }
        System.err.println(listMetabolites.size() + " elements have been mapped (against " + remainingMetabolites.size() + " non-mapped).\n");
        return listMetabolites;
    }

    public void writeOutput(ArrayList<HashMap<BioPathway, Double>> resultList, Set<BioPhysicalEntity> map, String outputFile) throws IOException{
        BufferedWriter f = new BufferedWriter(new FileWriter(new File(outputFile)));

        f.write("Pathway enrichment\tFisher's p-value\tBonferroni correction\tBenjamini-Hochberg correction\tHolm-Bonferroni correction\tMapped metabolites\tNb of mapped\tCoverage\n");
        HashMap<BioPathway, Double> result = resultList.get(0);
        ArrayList <Iterator> iteratorList = new ArrayList<Iterator>();

        for (int i = 1; i < 4; i++) {
            iteratorList.add(resultList.get(i).values().iterator());
        }

        for (Map.Entry<BioPathway, Double> entry : result.entrySet()) {
            BioPathway path = entry.getKey();
            String printed = entry.getKey().getName() + "\t" + round((double)entry.getValue()) + "\t";

            for (int i = 0; i < 3; i++) {
                printed += round((double) iteratorList.get(i).next()) + "\t";
            }

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

    public static String round(double value) {
        if(value < 0.01){
            String tmp = (String.valueOf(value));
            String[] splitted = tmp.split("E");
            double head = Double.parseDouble(splitted[0]);

            head =  (double) Math.round(head * 100)/ 100;
            try{
                return  (String.valueOf(head) + "E" + splitted[1] +" ");
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
                tmp = String.valueOf((double) Math.round(tail * 100)/ 100);
                return (tmp + " E-" + String.valueOf(j + 1));
            }
        }
        return String.valueOf((double) Math.round(value * 100) / 100);
    }

    public static void main( String[] args ){

        //Parameters
        long startTime = System.nanoTime();
        String dir = "/home/bmerlet/Documents/PathwayEnrichment/";
        String inputFile = dir + "Galaxy15-[Biosigner_Multivariate_Univariate_Multivariate_variableMetadata.tsv].tabular";
        String sbml = dir + "recon2.v03_ext_noCompartment_noTransport.xml";
        String outputFile = dir + "output.tsv";
        String[] inchiLayers = {"c","h"};

        //SBML parsing
        JSBMLToBionetwork jsbml = new JSBMLToBionetwork(sbml);
        BioNetwork bionet = jsbml.getBioNetwork();
        bionet.printBioNetworkSizeToOut();


        MetExplore4Galaxy app = new MetExplore4Galaxy();
        try{

           //Mapping
            HashMap <String, String[]> parsedFile = app.extractData(inputFile, false,-1);
            Set<BioPhysicalEntity> map = app.mapping(bionet, parsedFile, 4, inchiLayers);

            //PathwayEnrichment
            PathwayEnrichment enr = new PathwayEnrichment(bionet, map);
            ArrayList <HashMap<BioPathway, Double>> resultList = new ArrayList<HashMap<BioPathway, Double> >();

            resultList.add(enr.computeEnrichment());
            for (int i = 0; i < 3; i++) {
                resultList.add(enr.computeEnrichment(i));
            }

            System.err.println(resultList.get(0).size() + " pathways are concerned among the network (on a total of " + bionet.getPathwayList().size() + ").");
            app.writeOutput(resultList, map, outputFile);

        }
        catch (IOException e){
            e.printStackTrace();
        }

        long elapsedTime = System.nanoTime() - startTime;//TODO: convert it in a class method
        if (elapsedTime / 1000000000L == 0L) {
            System.err.println("Time to run the process in miliseconde :" + elapsedTime / 1000000L);
        } else {
            System.err.println("Time to run the process in seconde :" + elapsedTime / 1000000000L);
        }

    }
}