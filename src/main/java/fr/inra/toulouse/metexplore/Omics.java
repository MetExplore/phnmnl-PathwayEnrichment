package fr.inra.toulouse.metexplore;

import parsebionet.biodata.BioNetwork;
import parsebionet.biodata.BioEntity;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public abstract class Omics {

    protected Boolean ifGalaxy;
    protected String text4outputFileInfo;
    protected HashMap<String, String[]> list_fingerprint;//input file after formatting and filtering
    protected Set<BioEntity> list_mappedEntities; //list of mapped metabolites used for analysis
    //Set type is used to avoid metabolites duplicates (no need to used Set now, could be refactored)
    protected BioNetwork network;
    protected WritingComportment writingComportment;
    protected OmicsMethods methods;
    protected int bioEntityType;
    protected File log = new File("information.txt");
    static int nbInstance=0;


    public Omics (Boolean ifGalaxy, HashMap<String, String[]> list_fingerprint,
                       Set<BioEntity> list_mappedEntities, BioNetwork network, int bioEntityType){
        nbInstance++;
        this.ifGalaxy =ifGalaxy;
        this.text4outputFileInfo="";
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = list_mappedEntities;
        this.network = network;
        this.writingComportment = new WritingComportment();
        this.bioEntityType=bioEntityType;
        this.methods = new OmicsMethods(list_mappedEntities,network,bioEntityType);
    }

    public Omics (Boolean ifGalaxy, HashMap<String, String[]> list_fingerprint,
                  BioNetwork network, int bioEntityType){
        nbInstance++;
        this.ifGalaxy =ifGalaxy;
        this.text4outputFileInfo="";
        this.list_fingerprint = list_fingerprint;
        this.list_mappedEntities = new HashSet<BioEntity>();
        this.network = network;
        this.writingComportment = new WritingComportment();
        this.bioEntityType=bioEntityType;
        this.methods = new OmicsMethods(list_mappedEntities,network,bioEntityType);
    }

    public void writeOutputInfo() throws IOException {
        if (this.ifGalaxy) {//if "writing console output in a file" functionality is activated
            if(nbInstance==1){
                //avoid to erase the file if Omics class is called more than one time
                if(log.isFile()) {
                    //write a new file if already exists
                    log.delete();
                }
                log.createNewFile();
            }
            BufferedWriter b = new BufferedWriter(new FileWriter(log, true));
            b.write(text4outputFileInfo);
            b.close();
        }
    }

    public void writeLog(String message) {
        System.out.println(message.replaceAll("\n", ""));
        text4outputFileInfo += message;
    }
}