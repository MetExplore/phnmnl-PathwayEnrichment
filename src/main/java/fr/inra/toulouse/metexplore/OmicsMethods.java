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

    public HashSet<BioEntity> intersect(Collection<BioEntity> set1) {
        return intersect(set1,this.list_mappedEntities.keySet());
    }

    public HashSet<BioEntity> intersect(Collection<BioEntity> set1, Collection<BioEntity> set2) {
        HashSet<BioEntity> inter = new HashSet();
        for (BioEntity bpe : set2){
        if (set1.contains(bpe)) inter.add(bpe);
        }
        return inter;
    }
}
