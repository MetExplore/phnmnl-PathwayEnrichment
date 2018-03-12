package fr.inra.toulouse.metexplore;

import parsebionet.biodata.*;

import java.util.*;

public class OmicsMethods {

    protected HashMap<BioEntity, String> list_mappedEntities; //list of mapped metabolites used for analysis
    //Set type is used to avoid metabolites duplicates (no need to used Set now, could be refactored)
    protected BioNetwork network;
    protected int bioEntityTYpe;

    public OmicsMethods(HashMap<BioEntity, String> list_mappedEntities, BioNetwork network, int bioEntityType){
        this.list_mappedEntities = list_mappedEntities;
        this.network = network;
        this.bioEntityTYpe=bioEntityType;
    }

    public Collection getEntitySetInPathway(BioPathway pathway) {
       switch (this.bioEntityTYpe){
           case 1:
               return pathway.getListOfInvolvedMetabolite().values();
           case 2:
               return pathway.getReactions().values();
           case 4:
                ;
           case 5:
               Set <BioProtein> proteins = new HashSet<BioProtein>();
               for (BioGene g : pathway.getGenes()){
                   proteins.addAll(g.getProteinList().values());
               }
               return proteins;
           case 6:
            return pathway.getGenes();
       }
        return null;
    }

    public HashMap getEntitySetInNetwork() {
        switch (this.bioEntityTYpe){
            case 1:
                return this.network.getPhysicalEntityList();
            case 2:
                return this.network.getBiochemicalReactionList();
            case 3:
                return this.network.getPathwayList();
            case 4:
                return this.network.getEnzymeList();
            case 5:
                return this.network.getProteinList();
            case 6:
            return this.network.getGeneList();
        }
        return null;
    }
    public int[] getFisherTestParameters(BioPathway pathway) {
        System.out.println("Hi again !");
        Collection entityInPathway = this.getEntitySetInPathway(pathway);
        //nb of mapped in the pathway
        int a = this.intersect(entityInPathway).size();
        System.out.println(pathway.getName() + ": size: " + a);
        //unmapped metabolites in the fingerprint
        int b = this.list_mappedEntities.size() - a;
        //unmapped metabolites in the pathway
        int c = entityInPathway.size() - a;
        //remaining metabolites in the network
        int d = this.getEntitySetInNetwork().size() - (a + b + c);

        int fisherTestParameters[] = {a,b,c,d};
        return fisherTestParameters;
    }

    public HashSet<BioEntity> intersect(Collection<BioEntity> set2) {
        HashSet<BioEntity> inter = new HashSet();
        for (BioEntity bpe: this.list_mappedEntities.keySet()){
            if (set2.contains(bpe)) inter.add(bpe);
        }
        return inter;
    }
}
