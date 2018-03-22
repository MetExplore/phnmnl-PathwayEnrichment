package fr.inra.toulouse.metexplore;

import parsebionet.biodata.*;
import java.util.*;

public interface OmicsMethods {

     default HashMap getEntitySetInNetwork(BioNetwork network, int bioEntityType) {
        switch (bioEntityType){
            case 1:
                return network.getPhysicalEntityList();
            case 2:
                return network.getBiochemicalReactionList();
            case 3:
                return network.getPathwayList();
            case 4:
                return network.getEnzymeList();
            case 5:
                return network.getProteinList();
            case 6:
            return network.getGeneList();
        }
        return null;
    }

     default HashSet<BioEntity> intersect(Collection<? extends BioEntity> set1, Collection<BioEntity> set2) {
        HashSet<BioEntity> inter = new HashSet<>();
        for (BioEntity bpe : set2){
            if (set1.contains(bpe)) inter.add(bpe);
        }
        return inter;
    }

     default String getTypeOfEntity(int entityType){
        switch (entityType) {
            case 1:
                return "Metabolite";
            case 2:
                return "Reaction";
            case 3:
                return "Pathway";
            case 4:
                return "Enzyme";
            case 5:
                return "Protein";
            case 6:
                return "Gene";
        }
        return null;
    }
}
