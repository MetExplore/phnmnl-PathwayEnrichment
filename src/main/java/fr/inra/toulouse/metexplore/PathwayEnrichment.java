package fr.inra.toulouse.metexplore;

import parsebionet.biodata.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class PathwayEnrichment {
    public List<HashMap<BioPathway, Double>> list_pathwayEnr = new ArrayList<HashMap<BioPathway, Double>>(); //list of pathway containing mapped metabolites, p-value and corrections values
    public BioNetwork network;
    public Set <BioPhysicalEntity> list_mappedMetabolites;
    public String outFilePathEnr, text4outputFileInfo="";
    public Boolean ifGalaxy = false;
    
    public PathwayEnrichment (BioNetwork network, Set <BioPhysicalEntity> list_mappedMetabolites,
                              String outFilePathEnr, Boolean ifGalaxy) throws IOException {
        this.network=network;
        this.list_mappedMetabolites=list_mappedMetabolites;
        this.outFilePathEnr=outFilePathEnr;
        this.ifGalaxy=ifGalaxy;
        this.computeEnrichment();
    }

    public void computeEnrichment() throws IOException {

        System.out.println("Pathway enrichment in progress...");
        parsebionet.statistics.PathwayEnrichment pathEnr = new parsebionet.statistics.PathwayEnrichment(network, list_mappedMetabolites);
        HashMap<BioPathway, Double> pathEnrWhithPval = pathEnr.computeEnrichment(); //obtaining p-values for mapped pathway
        HashMap<BioPathway, Double> pathEnrBenHoc = pathEnr.benjaminiHochbergCorrection(pathEnrWhithPval);

        list_pathwayEnr.add(sortPathEnrByBenHocPath(pathEnrWhithPval, pathEnrBenHoc));//benjaminiHochberg function sorts biopath by pval, need to do the same here to join with it
        list_pathwayEnr.add(sortPathEnrByBenHocPath(pathEnr.bonferroniCorrection(pathEnrWhithPval), pathEnrBenHoc));
        list_pathwayEnr.add(pathEnrBenHoc);//same for Benjamini Hochberg

        writeLog(list_pathwayEnr.get(0).size() + " pathways are concerned among the network (on " + network.getPathwayList().size() + " in the network).");
        writeOutputPathEnr();
        writeOutputInfo();
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
        List <String> listPathwayMetabolitesID = new ArrayList<String>();

        f.write("Pathway_name\tFisher_p-value\tBonferroni_correction\tBenjamini-Hochberg_correction\tMapped_metabolites\tMapped_metabolites_ID\tNb. of mapped\tCoverage (%)");
        if (this.ifGalaxy)  f.write("\tNb. of unmapped in pathway\tNb. of unmapped in fingerprint\tNb. of remaining in network\n");
        else f.write("\n");

        HashMap<BioPathway, Double> result = list_pathwayEnr.get(0);//get pathway enrichment without corrections
        Iterator itBonCorr = list_pathwayEnr.get(1).values().iterator(); //Create a loop on Bonferonni values
        Iterator itBenHocCorr = list_pathwayEnr.get(2).values().iterator();//Same for Benjamini Hochberg

        for (Map.Entry<BioPathway, Double> pathEnrEntry : result.entrySet()) {//Loop on pathway enrichment without corrections
            listPathwayMetabolites = new ArrayList<String>();
            listPathwayMetabolitesID = new ArrayList<String>();
            BioPathway path = pathEnrEntry.getKey();

            int j = 0; //number of mapped metabolites contained in a BioPathway

            //Extracting metabolites from the mapping list contained in BioPathway
            for (BioPhysicalEntity bpe : list_mappedMetabolites) {
                if (path.getListOfInvolvedMetabolite().containsValue(bpe)) {
                    listPathwayMetabolites.add(bpe.getName());
                    listPathwayMetabolitesID.add(bpe.getId());
                    j++;
                }
            }
            //Collections.sort(listPathwayMetabolites);
            String coverage = round((double) j / (double) path.getListOfInvolvedMetabolite().size() * (double) 100);
            PathwayEnrichmentElement pathEnrElement = new PathwayEnrichmentElement(pathEnrEntry.getKey().getName(),pathEnrEntry.getValue(),(double)itBonCorr.next(),(double)itBenHocCorr.next(),listPathwayMetabolites,listPathwayMetabolitesID,j,coverage);
            if (this.ifGalaxy) pathEnrElement.settings4Galaxy(getFisherTestParameters(path, j));
            list_pathwayEnrElement.add(pathEnrElement);
        }
                    
        Collections.sort(list_pathwayEnrElement);
        for (int i=0;i< list_pathwayEnrElement.size();i++){
            f.write(list_pathwayEnrElement.get(i).toString());
        }
        if (f != null) {
            f.close();
        }
    }

    static Set<BioChemicalReaction> getReaction(Set<? extends BioEntity> BioEntitySet) {

        Set<BioChemicalReaction> reactionSet = new HashSet();
        Iterator it = BioEntitySet.iterator();

        while(it.hasNext()) {
            BioPhysicalEntity m = (BioPhysicalEntity)it.next();
            reactionSet.addAll(m.getReactionsAsProduct().values());
            reactionSet.addAll(m.getReactionsAsSubstrate().values());
        }
        return reactionSet;

    }

    public int[] getFisherTestParameters(BioPathway pathway, int nb_mapped) {
        Set<BioChemicalReaction> reactionSet=getReaction(list_mappedMetabolites);
        Collection<BioChemicalReaction> reactionInPathway = pathway.getReactions().values();
        int fisherTestParameters[] = new int[3];
        
        //unmapped metabolites in the fingerprint
        fisherTestParameters[0] = reactionSet.size() - nb_mapped;
        //unmapped metabolites in the pathway
        fisherTestParameters[1] = reactionInPathway.size() - nb_mapped;
        //remaining metabolites in the network
        fisherTestParameters[2] = network.getBiochemicalReactionList().size() - (nb_mapped + fisherTestParameters[0] + fisherTestParameters[1]);

        return fisherTestParameters;
    }

    public void writeOutputInfo() throws IOException {
        if (this.ifGalaxy) {//if "writing console output in a file" functionality is activated
            File f = new File("information.tsv");
            f.createNewFile();
            BufferedWriter b = new BufferedWriter(new FileWriter(f));
            b.write(text4outputFileInfo);
            b.close();
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

    public class PathwayEnrichmentElement implements Comparable <PathwayEnrichmentElement> {

        public String pathName, mappedMetabolites, mappedMetabolitesID, coverage;
        public double p_value, q_value_Bonf, q_value_BenHoc;
        public int nb_mapped, nb_unmappedInPathway, nb_unmappedInFingerprint, nb_remainingInNetwork;

        public PathwayEnrichmentElement(String pathName, double p_value, double q_value_Bonf, double q_value_BenHoc,
                                        List<String> mappedMetabolites, List<String> mappedMetabolitesID, int nb_mapped, String coverage) {
            this.pathName = pathName;
            this.p_value = p_value;
            this.q_value_Bonf = q_value_Bonf;
            this.q_value_BenHoc = q_value_BenHoc;
            this.mappedMetabolites = String.join(";", mappedMetabolites);
            this.mappedMetabolitesID = String.join(";", mappedMetabolitesID);
            this.nb_mapped = nb_mapped;
            this.coverage = coverage;
        }

        public void settings4Galaxy(int[] fisherTestParameters){
            this.nb_unmappedInPathway=fisherTestParameters[1];
            this.nb_unmappedInFingerprint=fisherTestParameters[0];
            this.nb_remainingInNetwork=fisherTestParameters[2];
        }

        public int compareTo(PathwayEnrichmentElement p) {
            return (this.pathName).compareToIgnoreCase(p.pathName);
        }

        public String toString() {
           String line  = this.pathName + "\t" + removeSciNot(this.p_value) + "\t" + removeSciNot(this.q_value_Bonf) + "\t" + removeSciNot(this.q_value_BenHoc) + "\t" + this.mappedMetabolites.toString() + "\t" + this.mappedMetabolitesID.toString() + "\t" + this.nb_mapped + "\t" + this.coverage;
            if (ifGalaxy) line += "\t" + this.nb_unmappedInPathway + "\t" + this.nb_unmappedInFingerprint + "\t" + this.nb_remainingInNetwork;
            return line + "\n";
        }
    }

}