package fr.inra.toulouse.metexplore;

import parsebionet.biodata.*;

import java.util.*;

public class PathwayEnrichmentCalculation extends parsebionet.statistics.PathwayEnrichment {

    protected  BioNetwork network;
    protected  HashMap <BioEntity, String> list_mappedEntities;
    protected  OmicsMethods methods;
    protected  Set<BioChemicalReaction> reactionSet;

    public PathwayEnrichmentCalculation(BioNetwork network, HashMap <BioEntity, String> list_mappedEntities, int bioEntityType){
        super(network, list_mappedEntities.keySet());
        this.network = network;
        this.list_mappedEntities = list_mappedEntities;
        this.methods = new OmicsMethods(list_mappedEntities,network, bioEntityType);
        this.setReactionSet(list_mappedEntities.keySet());
    }

    public void setReactionSet(Set<? extends BioEntity> BioEntitySet) {
        this.reactionSet = new HashSet();
        for(BioEntity e : BioEntitySet) {
            if (e instanceof BioChemicalReaction) {
                this.reactionSet.add((BioChemicalReaction) e);
            /*}else if(this.methods.bioEntityTYpe == 4) {
                for (BioChemicalReaction r : network.getBiochemicalReactionList().values()) {
                    if (r.getEnzList().values().contains(e)) {
                        System.out.println("\t"  + ". " + r.getId());
                        for (BioPathway p : network.getPathwayList().values()) {
                            if(p.getReactions().values().contains(r)){
                                System.out.println("\t\t" + p.getId());
                            }
                        }
                        reactionSet.add(r);
                    }
                }*/
            //BUG: TODO: see in JSBML2BioNetwork why only one enzyme is associated to a reaction
            // (instead of multiple for protein and genes)
            }else if(this.methods.bioEntityTYpe == 4 || e instanceof BioProtein) {
                BioProtein p1 = (BioProtein) e;
                BioProtein p = network.getProteinList().get(e.getId());
                HashMap<String, BioGene> list_genes = p.getGeneList();
                for (BioGene g : list_genes.values()){
                    addReactions(network.getReactionsFromGene(g.getId()));
                }
            } else if (e instanceof BioGene) {
                BioGene g = (BioGene) e;
                addReactions(network.getReactionsFromGene(g.getId()));
            } else if (e instanceof BioPhysicalEntity) {
                if (this.methods.bioEntityTYpe == 1) {
                    BioPhysicalEntity m = (BioPhysicalEntity) e;
                    this.reactionSet.addAll(m.getReactionsAsProduct().values());
                    this.reactionSet.addAll(m.getReactionsAsSubstrate().values());
                }
            }
        }
    }

    public void addReactions(Set <String> list_reactions_ID){

        for (String reac_ID : list_reactions_ID) {
            BioChemicalReaction r = network.getBiochemicalReactionList().get(reac_ID);
            for (BioPathway p : network.getPathwayList().values()) {
                if(p.getReactions().values().contains(r)){
                }
            }
            this.reactionSet.add(r);
        }
    }

    @Override
    public HashMap<BioPathway, Double> computeEnrichment() {
        HashSet<BioPathway> paths = new HashSet();

       for ( BioChemicalReaction r : this.reactionSet){
            //System.out.println("Reac: " + r.getName() + ": " + r.getListOfGeneNames().size());
            paths.addAll(r.getPathwayList().values());
        }

       HashMap<BioPathway, Double> res = new HashMap();
       for (BioPathway p : paths){
            //System.out.println("PathName: " + p.getName() + ": " + p.getGenes().size());
            res.put(p, this.getPvalue(p));
        }
        return res;
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

        for(int k = 0; k < orderedPaths.size(); ++k) {
            p = (BioPathway)orderedPaths.get(k);
            pval = ((Double)pvalues.get(p)).doubleValue();
           // System.out.println("Pval : " + p.getName() + ": " + pval);

            //case of probability equality
            if (k+1 < orderedPaths.size()){
                p_next = (BioPathway)orderedPaths.get(k+1);
                double pval_next = ((Double)pvalues.get(p_next)).doubleValue();
                if (pval_next == pval) {
                    //System.out.println(true);
                    adjPval = pval * (double) pvalues.size() / ((double) k + (double) 1.5);
                    adjPvalues.put(p, new Double(adjPval));
                    adjPvalues.put(p_next, new Double(adjPval));
                    continue;
                }
            }

            //default case
            adjPval = pval * (double) pvalues.size() / (double) (k + 1);
            adjPvalues.put(p, new Double(adjPval));
        }
        return adjPvalues;
    }
}