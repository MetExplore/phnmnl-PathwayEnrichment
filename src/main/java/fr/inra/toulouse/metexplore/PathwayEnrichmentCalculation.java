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
        int i=0, j;
        for(BioEntity e : BioEntitySet) {
            i++;j=0;
            System.out.println("\n" + i + ". " + e.getId());
            if (e instanceof BioChemicalReaction) {
                this.reactionSet.add((BioChemicalReaction) e);
            }else if(this.methods.bioEntityTYpe == 4) {
                for (BioChemicalReaction r : network.getBiochemicalReactionList().values()) {
                    if (r.getEnzList().values().contains(e)) {
                        j++;
                        System.out.println("\t" + j + ". " + r.getId());
                        for (BioPathway p : network.getPathwayList().values()) {
                            if(p.getReactions().values().contains(r)){
                                System.out.println("\t\t" + p.getId());
                            }
                        }
                        reactionSet.add(r);
                    }
                }
            } else if (e instanceof BioGene) {
                System.out.println("Genes");
                BioGene g = (BioGene) e;
                addReactions(network.getReactionsFromGene(g.getId()));
            } else if (e instanceof BioProtein) {
                System.out.println("Protein");
                BioProtein p = (BioProtein) e;
                HashMap<String, BioGene> list_genes = p.getGeneList();
                for (BioGene g : list_genes.values()){
                    addReactions(network.getReactionsFromGene(g.getId()));
                }
            } else if (e instanceof BioPhysicalEntity) {
                if (this.methods.bioEntityTYpe == 1) {
                    System.out.println("Metabolite");
                    BioPhysicalEntity m = (BioPhysicalEntity) e;
                    this.reactionSet.addAll(m.getReactionsAsProduct().values());
                    this.reactionSet.addAll(m.getReactionsAsSubstrate().values());
                }
            }
        }
        System.out.println("ListReacSize: " + this.reactionSet.size());
    }

    public void addReactions(Set <String> list_reactions_ID){

        //System.out.println("sizeReac: " + list_reactions_ID.size());
        int j = 0;
        for (String reac_ID : list_reactions_ID) {
            j++;
            BioChemicalReaction r = network.getBiochemicalReactionList().get(reac_ID);
            //System.out.println("Reac: " + r.getName());
            System.out.println("\t" + j + ". " + r.getId());
            for (BioPathway p : network.getPathwayList().values()) {
                if(p.getReactions().values().contains(r)){
                    System.out.println("\t\t" + p.getId());
                }
            }
            this.reactionSet.add(r);
        }
    }

    @Override
    public HashMap<BioPathway, Double> computeEnrichment() {
        HashSet<BioPathway> paths = new HashSet();

        System.out.println("SizeCompute: " + this.reactionSet.size());
       for ( BioChemicalReaction r : this.reactionSet){
            //System.out.println("Reac: " + r.getName() + ": " + r.getListOfGeneNames().size());
            paths.addAll(r.getPathwayList().values());
        }

        System.out.println("PathSize: " + paths.size());

       HashMap<BioPathway, Double> res = new HashMap();
       for (BioPathway p : paths){
            //System.out.println("PathName: " + p.getName() + ": " + p.getGenes().size());
            res.put(p, this.getPvalue(p));
        }
        System.out.println("ResSIze: "+ res.size());
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