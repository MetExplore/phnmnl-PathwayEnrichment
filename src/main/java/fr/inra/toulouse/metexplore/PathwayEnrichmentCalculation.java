package fr.inra.toulouse.metexplore;

import parsebionet.biodata.*;

import java.util.*;

public class PathwayEnrichmentCalculation extends parsebionet.statistics.PathwayEnrichment {

    protected  BioNetwork network;
    protected  HashMap <BioEntity, String> list_mappedEntities;
    protected  OmicsMethods methods;
    protected Set<BioChemicalReaction> reactionSet;

    public PathwayEnrichmentCalculation(BioNetwork network, HashMap <BioEntity, String> list_mappedEntities, int bioEntityType){
        super(network, list_mappedEntities.keySet());
        this.network = network;
        this.list_mappedEntities = list_mappedEntities;
        this.methods = new OmicsMethods(list_mappedEntities,network, bioEntityType);
        this.setReactionSet(network,list_mappedEntities.keySet());
    }

    public void setReactionSet(BioNetwork network, Set<? extends BioEntity> BioEntitySet) {
        this.reactionSet = new HashSet();

        for(BioEntity e : BioEntitySet) {
            if (e instanceof BioChemicalReaction) {
                this.reactionSet.add((BioChemicalReaction)e);
            } else if (e instanceof BioPhysicalEntity) {
                BioPhysicalEntity m = (BioPhysicalEntity)e;
                this.reactionSet.addAll(m.getReactionsAsProduct().values());
                this.reactionSet.addAll(m.getReactionsAsSubstrate().values());
            } else if (e instanceof BioGene) {
                BioGene g = (BioGene) e;
                addReactionsFromGene(g);
            } else if (e instanceof BioProtein) {
                 BioProtein p = (BioProtein) e;
                 HashMap<String, BioGene> list_genes = p.getGeneList();
                 for (BioGene g : list_genes.values()){
                     addReactionsFromGene(g);
                 }
            }
        }
    }

    public void addReactionsFromGene(BioGene g){
        Set <String> list_reactions_ID = network.getReactionsFromGene(g.getId());
        for (String reac_ID : list_reactions_ID) {
            BioChemicalReaction r = network.getBiochemicalReactionList().get(reac_ID);
            this.reactionSet.add(r);
        }
    }


    @Override
    public double getPvalue(BioPathway pathway) throws IllegalArgumentException {
        if (!this.network.getPathwayList().values().contains(pathway)) {
            throw new IllegalArgumentException("pathway not in network");
        } else {
            int fisherTestParameters[] = this.methods.getFisherTestParameters(pathway);
            return this.exactFisherOneTailed(fisherTestParameters[0], fisherTestParameters[1],
                    fisherTestParameters[2], fisherTestParameters[3]);
        }
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
}