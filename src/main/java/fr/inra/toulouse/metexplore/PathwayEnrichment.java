package fr.inra.toulouse.metexplore;

import parsebionet.biodata.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class PathwayEnrichment extends Omics{

    protected List<HashMap<BioPathway, Double>> list_pathwayEnr = new ArrayList<HashMap<BioPathway, Double>>(); //list of pathway containing mapped metabolites, p-value and corrections values
    protected String outFilePathEnr;

    public PathwayEnrichment (BioNetwork network, HashMap<String, String[]> list_fingerprint, Set <BioEntity> list_mappedEntities,
                              String outFilePathEnr, Boolean ifGalaxy , int bioEntityType) throws IOException {
        super(ifGalaxy, list_fingerprint, list_mappedEntities, network, bioEntityType);
        this.outFilePathEnr=outFilePathEnr;
        this.computeEnrichmentWithCorrections();
    }

    public void computeEnrichmentWithCorrections() throws IOException {
        System.out.println("Pathway enrichment in progress...");
        parsebionet.statistics.PathwayEnrichment pathEnr;
        BioEntity bpe = this.list_mappedEntities.iterator().next();
        if(bpe instanceof BioChemicalReaction){
            pathEnr = new parsebionet.statistics.PathwayEnrichment(this.network, this.list_mappedEntities);
        }else {
            pathEnr = new PathwayEnrichmentCalculation(this.network, this.list_mappedEntities,this.bioEntityType);
        }
        HashMap<BioPathway, Double> pathEnrWhithPval = pathEnr.computeEnrichment(); //obtaining p-values for mapped pathway
        HashMap<BioPathway, Double> pathEnrBenHoc = pathEnr.benjaminiHochbergCorrection(pathEnrWhithPval);

        this.list_pathwayEnr.add(sortPathEnrByBenHocPath(pathEnrWhithPval, pathEnrBenHoc));//benjaminiHochberg function sorts biopath by pval,
        // need to do the same here to join with it
        this.list_pathwayEnr.add(sortPathEnrByBenHocPath(pathEnr.bonferroniCorrection(pathEnrWhithPval), pathEnrBenHoc));
        this.list_pathwayEnr.add(pathEnrBenHoc);//same for Benjamini Hochberg
        writeLog(this.list_pathwayEnr.get(0).size() + " pathways are concerned among the network (on " + this.network.getPathwayList().size() + " in the network).");
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

        File fo2 = new File(this.outFilePathEnr);
        fo2.createNewFile();
        BufferedWriter f = new BufferedWriter(new FileWriter(fo2));
        List <PathwayEnrichmentElement> list_pathwayEnrElement = new ArrayList<PathwayEnrichmentElement>(); //list of pathway enrichment instantiation for sorting
        List <String> listPathwayMetabolites = new ArrayList<String>();
        List <String> listPathwayMetabolitesID = new ArrayList<String>();

        f.write("Pathway_name\tFisher_p-value\tBonferroni_correction\tBenjamini-Hochberg_correction\tMapped_metabolites\t" +
                "Mapped_metabolites_ID\tNb. of mapped\tCoverage (%)");
        if (this.ifGalaxy)  f.write("\tNb. of unmapped in pathway\tNb. of unmapped in fingerprint\tNb. of remaining in network\n");
        else f.write("\n");

        HashMap<BioPathway, Double> result = this.list_pathwayEnr.get(0);//get pathway enrichment without corrections
        Iterator itBonCorr = this.list_pathwayEnr.get(1).values().iterator(); //Create a loop on Bonferonni values
        Iterator itBenHocCorr = this.list_pathwayEnr.get(2).values().iterator();//Same for Benjamini Hochberg

        for (Map.Entry<BioPathway, Double> pathEnrEntry : result.entrySet()) {//Loop on pathway enrichment without corrections
            listPathwayMetabolites = new ArrayList<String>();
            listPathwayMetabolitesID = new ArrayList<String>();
            BioPathway path = pathEnrEntry.getKey();

            int j = 0; //number of mapped metabolites contained in a BioPathway

            //Extracting metabolites from the mapping list contained in BioPathway
            for (BioEntity bpe : this.list_mappedEntities) {
                if (path.getReactions().containsValue(bpe) || path.getGenes().contains(bpe) ||
                        path.getListOfInvolvedMetabolite().containsValue(bpe)) {
                    listPathwayMetabolites.add(bpe.getName());
                    listPathwayMetabolitesID.add(bpe.getId());
                    j++;
                }
            }
            //Collections.sort(listPathwayMetabolites);
            String coverage = this.writingComportment.round((double) j / (double) getEntityPerPathwaySize(path) * (double) 100);
            PathwayEnrichmentElement pathEnrElement = new PathwayEnrichmentElement(pathEnrEntry.getKey().getName(),pathEnrEntry.getValue(),
                    (double)itBonCorr.next(),(double)itBenHocCorr.next(),listPathwayMetabolites,listPathwayMetabolitesID,j,coverage, this.ifGalaxy);
            if (this.ifGalaxy) pathEnrElement.settings4Galaxy(this.getFisherTestParameters(path, j));
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

    public int getEntityPerPathwaySize(BioPathway pathway) {
        BioEntity bpe = this.list_mappedEntities.iterator().next();
        if (bpe instanceof BioChemicalReaction) {
            return pathway.getReactions().size();
        } else if (bpe instanceof BioGene) {
            return pathway.getGenes().size();
        } else {
            return pathway.getListOfInvolvedMetabolite().size();
        }
    }

    public int getEntityPerNetworkSize() {
        BioEntity bpe = this.list_mappedEntities.iterator().next();
        if (bpe instanceof BioChemicalReaction) {
            return this.network.getBiochemicalReactionList().size();
        } else if (bpe instanceof BioGene) {
            return this.network.getGeneList().size();
        } else {
            return this.network.getPhysicalEntityList().size();
        }
    }


    public int[] getFisherTestParameters(BioPathway pathway,  int nbMapped) {

        //nb of mapped in the pathway
        int a = nbMapped;
        //unmapped metabolites in the fingerprint
        int b = this.list_mappedEntities.size() - a;
        //unmapped metabolites in the pathway
        int c = this.getEntityPerPathwaySize(pathway) - a;
        //remaining metabolites in the network
        int d = this.getEntityPerNetworkSize() - (a + b + c);

        int fisherTestParameters[] = {b,c,d};
        return fisherTestParameters;
    }
}