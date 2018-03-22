package fr.inra.toulouse.metexplore.launcher;

import java.io.IOException;

import fr.inra.toulouse.metexplore.io.WritingComportment;
import fr.inra.toulouse.metexplore.omics.Omics;
import fr.inra.toulouse.metexplore.omics.PathwayEnrichment;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.kohsuke.args4j.CmdLineException;

public class Launcher_PathEnr extends Launcher_mapping implements WritingComportment{

    @Option(name = "-o3", aliases = "--outPath", usage = "Output file name for pathway enrichment result (by default: pathwayEnrichment.tsv).")
    protected String outFilePathEnr = "pathwayEnrichment.tsv";

    @Option(name = "-tEnr", aliases = "--typeEnr", usage = "1 for metabolites, 2 for reactions, 3 for pathway, 4 for enzyme, 5 for protein, 6 for gene (by default: pathways).")
    protected int entityType2Enrich = 3;

    @SuppressWarnings("deprecation")
    public void printInfo(CmdLineParser parser, String[] args) throws CmdLineException {

        super.printInfo(parser, args);

        if (this.entityType2Enrich < 1 || this.entityType2Enrich > 6) {
            throw new CmdLineException("Type of enriched entity be between 1 and 6.");
        }
    }

    public Omics analyse(CmdLineParser parser, String[] args) throws IOException{

        Omics mapping = super.analyse(parser, args);

        PathwayEnrichment enr = new PathwayEnrichment(logContent, network, fingerprint.getList_entities(), mapping.getList_mappedEntities(),
                    this.outFilePathEnr, this.galaxyFile, this.entityType2Map, this.entityType2Enrich);

        this.logContent = enr.getLogContent();
        return enr;
    }

    public static void main(String[] args) {
        exec(new Launcher_PathEnr(), args);
    }

}