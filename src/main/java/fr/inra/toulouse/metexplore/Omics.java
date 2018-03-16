package fr.inra.toulouse.metexplore;

import parsebionet.biodata.BioGene;
import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioEntity;
import parsebionet.biodata.BioProtein;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public abstract class Omics {
    protected ArrayList<String[]> list_fingerprint;//input file after formatting and filtering
    protected HashMap<BioEntity, String> list_mappedEntities; //list of mapped metabolites used for analysis
    protected BioNetwork network;
    protected WritingComportment write;
    protected OmicsMethods omics;
    protected int bioEntityType;
    protected String typeOfMappedEntity;


    public Omics (String galaxy, ArrayList<String[]> list_fingerprint,
                       HashMap<BioEntity, String> list_mappedEntities, BioNetwork network, int bioEntityType){
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = list_mappedEntities;
        this.network = network;
        this.write = new WritingComportment(galaxy);
        this.bioEntityType=bioEntityType;
        this.omics = new OmicsMethods(list_mappedEntities,network,bioEntityType);
        this.typeOfMappedEntity = omics.getTypeOfEntity(bioEntityType).toLowerCase();
    }

    public Omics (String galaxy, ArrayList<String[]> list_fingerprint,
                  BioNetwork network, int bioEntityType){
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = new HashMap<BioEntity, String>();
        this.network = network;
        this.write = new WritingComportment(galaxy);
        this.bioEntityType=bioEntityType;
        this.omics = new OmicsMethods(list_mappedEntities,network,bioEntityType);
        this.typeOfMappedEntity = omics.getTypeOfEntity(bioEntityType).toLowerCase();
    }

    public void parseProtList4Genes(){
        HashMap<String,BioProtein> protList;
        for (BioGene g : network.getGeneList().values()) {
            protList = new HashMap();
            for (BioProtein p : network.getProteinList().values()) {
                if(p.getGeneList().values().contains(g)) protList.put(p.getId(),p);
            }
            g.setProteinList(protList);
        }
    }
}