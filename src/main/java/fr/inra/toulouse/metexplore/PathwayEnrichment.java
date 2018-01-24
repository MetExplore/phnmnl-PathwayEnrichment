package fr.inra.toulouse.metexplore;

import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPathway;
import parsebionet.biodata.BioPhysicalEntity;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class PathwayEnrichment {
    public List<HashMap<BioPathway, Double>> list_pathwayEnr = new ArrayList<HashMap<BioPathway, Double>>(); //list of pathway containing mapped metabolites, p-value and corrections values
    public BioNetwork network;
    public Set <BioPhysicalEntity> list_mappedMetabolites;
    public String outFilePathEnr, outFileInfo, text4outputFileInfo="";
    
    public PathwayEnrichment (BioNetwork network, Set <BioPhysicalEntity> list_mappedMetabolites, String outFileInfo, String outFilePathEnr){
        this.network=network;
        this.list_mappedMetabolites=list_mappedMetabolites;
        this.outFilePathEnr=outFilePathEnr;
        this.outFileInfo=outFileInfo;
    }

    public void computeEnrichment() throws IOException {

        System.out.println("Pathway enrichment in progress...");
        parsebionet.statistics.PathwayEnrichment pathEnr = new parsebionet.statistics.PathwayEnrichment(network, list_mappedMetabolites);
        HashMap<BioPathway, Double> pathEnrWhithPval = pathEnr.computeEnrichment(); //obtaining p-values for mapped pathway
        HashMap<BioPathway, Double> pathEnrBenHoc = pathEnr.benjaminiHochbergCorrection(pathEnrWhithPval);

        list_pathwayEnr.add(sortPathEnrByBenHocPath(pathEnrWhithPval, pathEnrBenHoc));//benjaminiHochberg function sorts biopath by pval, need to do the same here to join with it
        list_pathwayEnr.add(sortPathEnrByBenHocPath(pathEnr.bonferroniCorrection(pathEnrWhithPval), pathEnrBenHoc));
        list_pathwayEnr.add(pathEnrBenHoc);//same for Benjamini Hochberg

        writeLog(list_pathwayEnr.get(0).size() + " pathways are concerned among the network (on " + network.getPathwayList().size() + ").");
        writeOutputPathEnr();
    }

    public HashMap<BioPathway, Double> sortPathEnrByBenHocPath(HashMap<BioPathway, Double> disorderedPathEnr, HashMap<BioPathway, Double> pathEnrBenHoc) {
        ArrayList<BioPathway> pathBenHoc = new ArrayList(pathEnrBenHoc.keySet());
        HashMap<BioPathway, Double> orderedPathEnr = new HashMap();

        for (int i = 0; i < pathBenHoc.size(); ++i) {
            BioPathway path = (BioPathway) pathBenHoc.get(i);
            double pval = (Double) disorderedPathEnr.get(path);
            orderedPathEnr.put(path, pval);
        }
        return orderedPathEnr;
    }

    public void writeOutputPathEnr() throws IOException{

        File fo2 = new File(outFilePathEnr);
        fo2.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo2));
        List <PathwayEnrichmentElement> list_pathwayEnrElement = new ArrayList<PathwayEnrichmentElement>(); //list of pathway enrichment instantiation for sorting
        List <String> listPathwayMetabolites = new ArrayList<String>();

        f.write("Pathway enrichment\tFisher's p-value\tBonferroni correction\tBenjamini-Hochberg correction\tMapped metabolites\tNb of mapped\tCoverage (%)\n");
        HashMap<BioPathway, Double> result = list_pathwayEnr.get(0);//get pathway enrichment without corrections
        Iterator itBonCorr = list_pathwayEnr.get(1).values().iterator(); //Create a loop on Bonferonni values
        Iterator itBenHocCorr = list_pathwayEnr.get(2).values().iterator();//Same for Benjamini Hochberg

        for (Map.Entry<BioPathway, Double> pathEnrEntry : result.entrySet()) {//Loop on pathway enrichment without corrections
            listPathwayMetabolites = new ArrayList<String>();
            BioPathway path = pathEnrEntry.getKey();

            int j = 0; //number of mapped metabolites contained in a BioPathway

            //Extracting metabolites from the mapping list contained in BioPathway
            for (BioPhysicalEntity bpe : list_mappedMetabolites) {
                if (path.getListOfInvolvedMetabolite().containsValue(bpe)) {
                    listPathwayMetabolites.add(bpe.getName());
                    j++;
                }
            }
            Collections.sort(listPathwayMetabolites);
            String coverage = round((double) j / (double) path.getListOfInvolvedMetabolite().size() * (double) 100);
            list_pathwayEnrElement.add(new PathwayEnrichmentElement(pathEnrEntry.getKey().getName(),pathEnrEntry.getValue(),(double)itBonCorr.next(),(double)itBenHocCorr.next(),listPathwayMetabolites,j,coverage));
        }
        Collections.sort(list_pathwayEnrElement);
        for (int i=0;i< list_pathwayEnrElement.size();i++){
            f.write(list_pathwayEnrElement.get(i).toString());
        }
        if (f != null) {
            f.close();
        }
    }

    public class PathwayEnrichmentElement implements Comparable <PathwayEnrichmentElement> {

        public String pathName, mappedMetabolites, coverage;
        public double p_value, q_value_Bonf, q_value_BenHoc;
        public int nb_mapped;

        public PathwayEnrichmentElement(String pathName, double p_value, double q_value_Bonf, double q_value_BenHoc,
                                        List<String> mappedMetabolites, int nb_mapped, String coverage) {
            this.pathName = pathName;
            this.p_value = p_value;
            this.q_value_Bonf = q_value_Bonf;
            this.q_value_BenHoc = q_value_BenHoc;
            this.mappedMetabolites = String.join(";", mappedMetabolites);
            this.nb_mapped = nb_mapped;
            this.coverage = coverage;
        }

        public int compareTo(PathwayEnrichmentElement p) {
            return (this.pathName).compareToIgnoreCase(p.pathName);
        }

        public String toString() {
            return (this.pathName + "\t" + removeSciNot(this.p_value) + "\t" + removeSciNot(this.q_value_Bonf) + "\t" + removeSciNot(this.q_value_BenHoc) + "\t" + this.mappedMetabolites.toString() + "\t" + this.nb_mapped + "\t" + this.coverage + "\n");
        }
    }

    public void writeLog(String message) {
        System.out.println(message.replaceAll("\n", ""));
        text4outputFileInfo += message;
    }

    public String removeSciNot(double value) {
        String tmp = (String.valueOf(value));
        if (value < 1e-3) {
            String[] splitting = tmp.split("E");
            String power = splitting[1];
            if (power.charAt(1) == '0') {//remove the O after the E if exists
                power = power.substring(2, power.length());
            } else power = power.substring(1, power.length());
            String[] number = splitting[0].split("\\.");//obtain the integer and the decimal parts
            return "0." + (new String(new char[Integer.parseInt(power) - 1]).replace("\0", "0")) + number[0] + number[1];
        }
        return tmp;
    }

    public String round(double value) {
        return String.valueOf((double) Math.round(value * 100) / 100);
    }

}