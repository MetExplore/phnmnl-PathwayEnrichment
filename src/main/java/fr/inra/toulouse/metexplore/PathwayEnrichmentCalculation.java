package fr.inra.toulouse.metexplore;

import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioPathway;
import parsebionet.biodata.BioEntity;
import parsebionet.biodata.BioPhysicalEntity;

import java.util.*;

public class PathwayEnrichmentCalculation extends parsebionet.statistics.PathwayEnrichment {

    BioNetwork network;
    Set <BioEntity> list_mappedMetabolites;

    protected PathwayEnrichmentCalculation(BioNetwork network, Set <BioEntity> list_mappedMetabolites){
        super(network, list_mappedMetabolites);
        this.network = network;
        this.list_mappedMetabolites = list_mappedMetabolites;
    }

    @Override
    public double getPvalue(BioPathway pathway) throws IllegalArgumentException {
        if (!this.network.getPathwayList().values().contains(pathway)) {
            throw new IllegalArgumentException("pathway not in network");
        } else {
            int fisherTestParameters[] = this.getFisherTestParameters(pathway);
            return this.exactFisherOneTailed(fisherTestParameters[0], fisherTestParameters[1],
                    fisherTestParameters[2], fisherTestParameters[3]);
        }
    }

    public int[] getFisherTestParameters(BioPathway pathway) {
        Collection<BioPhysicalEntity> metaboliteInPathway = pathway.getListOfInvolvedMetabolite().values();

        //nb of mapped in the pathway
        int a = this.intersect(metaboliteInPathway).size();
        //unmapped metabolites in the fingerprint
        int b = this.list_mappedMetabolites.size() - a;
        //unmapped metabolites in the pathway
        int c = metaboliteInPathway.size() - a;
        //remaining metabolites in the network
        int d = this.network.getPhysicalEntityList().size() - (a + b + c);

        int fisherTestParameters[] = {a,b,c,d};
        return fisherTestParameters;
    }

    @Override
    public HashMap<BioPathway, Double> benjaminiHochbergCorrection(HashMap<BioPathway, Double> pvalues) {
        ArrayList<BioPathway> orderedPaths = this.sortPval(pvalues);
        HashMap<BioPathway, Double> adjPvalues = new HashMap();
        BioPathway p, p_next;
        double pval, adjPval;
        Boolean p_valEqual = false;

        for(int k = 0; k < orderedPaths.size(); ++k) {
            p = (BioPathway)orderedPaths.get(k);
            pval = ((Double)pvalues.get(p)).doubleValue();

            //case of probability equality
            if (k+1 < orderedPaths.size()){
                p_next = (BioPathway)orderedPaths.get(k+1);
                double pval_next = ((Double)pvalues.get(p_next)).doubleValue();
                if (pval_next == pval){
                    adjPval = pval * (double)pvalues.size() / ((double) k + (double)1.5);
                    adjPvalues.put(p, new Double(adjPval));
                    adjPvalues.put(p_next, new Double(adjPval));
                    p_valEqual = true;
                    continue;
                }
            }
            //default case
            if (!p_valEqual) {
                adjPval = pval * (double) pvalues.size() / (double) (k + 1);
                adjPvalues.put(p, new Double(adjPval));
                p_valEqual = false;
            }
        }
        return adjPvalues;
    }

    public HashSet<BioEntity> intersect(Collection<BioPhysicalEntity> set2) {
        HashSet<BioEntity> inter = new HashSet();
        for (BioEntity bpe: this.list_mappedMetabolites){
            if (set2.contains(bpe)) inter.add(bpe);
        }
        return inter;
    }
}