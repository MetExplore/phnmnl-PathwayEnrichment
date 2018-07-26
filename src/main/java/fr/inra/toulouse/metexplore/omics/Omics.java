/*******************************************************************************
 * Copyright INRA
 *
 *  Contact: ludovic.cottret@toulouse.inra.fr
 *
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *  In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *  The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *******************************************************************************/

package fr.inra.toulouse.metexplore.omics;

import fr.inra.toulouse.metexplore.io.WritingComportment;
import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioEntity;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

public abstract class Omics implements WritingComportment{
    protected ArrayList<String[]> list_fingerprint;//input file after formatting and filtering

    public HashMap<BioEntity, String> getList_mappedEntities() {
        return list_mappedEntities;
    }

    protected HashMap<BioEntity, String> list_mappedEntities; //list of mapped metabolites used for analysis
    protected BioNetwork network;

    protected int entityType2Map;
    protected String typeOfMappedEntity, galaxy, logContent;

    public String getLogContent() {
        return logContent;
    }

    public Omics (String logContent, String galaxy, ArrayList<String[]> list_fingerprint,
                       HashMap<BioEntity, String> list_mappedEntities, BioNetwork network, int entityType2Map){
        this.logContent = logContent;
        this.galaxy = galaxy;
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = list_mappedEntities;
        this.network = network;
        this.entityType2Map=entityType2Map;
        this.typeOfMappedEntity = getTypeOfEntity(entityType2Map).toLowerCase();
    }

    public Omics (String logContent, String galaxy, ArrayList<String[]> list_fingerprint,
                  BioNetwork network, int entityType2Map){
        this.logContent = logContent;
        this.galaxy = galaxy;
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = new HashMap<>();
        this.network = network;
        this.entityType2Map=entityType2Map;
        this.typeOfMappedEntity = getTypeOfEntity(entityType2Map).toLowerCase();
    }

    public HashMap getEntitySetInNetwork(int bioEntityType) {
        switch (bioEntityType){
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

    public HashMap getEntitySetInNetwork() {
       return getEntitySetInNetwork(this.entityType2Map);
    }

    public HashSet<BioEntity> intersect(Collection<BioEntity> set2) {
        HashSet<BioEntity> inter = new HashSet<>();
        for (BioEntity bpe : set2){
            if (this.list_mappedEntities.keySet().contains(bpe)) inter.add(bpe);
        }
        return inter;
    }

    public String getTypeOfEntity(int entityType){
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

    public void writeLog(String warning){
        this.logContent = writeLog(this.logContent,warning);
    }
}