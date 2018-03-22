package fr.inra.toulouse.metexplore.launcher;

import fr.inra.toulouse.metexplore.io.Fingerprint;
import fr.inra.toulouse.metexplore.io.JSBML2Bionetwork;
import fr.inra.toulouse.metexplore.omics.Mapping;
import fr.inra.toulouse.metexplore.omics.Omics;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import parsebionet.biodata.BioNetwork;
import java.io.IOException;

public class Launcher_mapping extends Launcher_Fingerprint {

    @Option(name = "-o2", aliases = "--outMap", usage = "Output file name for mapping result (by default: mapping.tsv).")
    protected String outFileMapping = "mapping.tsv";

    @Option(name = "-s", aliases = "--sbml", usage = "SBML file name (by default: Recon v2.02).")
    protected String sbml = "data/recon2.02_without_compartment.xml";

    /*****MAPPING PARAMETERS*****/

    @Option(name = "-t", aliases = "--type", usage = "1 for metabolites, 2 for reactions, 3 for pathway, 4 for enzyme, 5 for protein, 6 for gene (by default: metabolites).")
    protected int entityType2Map = 1;

    @Option(name = "-name", usage = "Activate this option for a name mapping .")
    protected Boolean nameMapping = false;

    @Option(name = "-prec", aliases = "--precision", usage = "Indicate the allowed error in ppm (used in mass mapping).")
    protected int weightPrecision = 2;

    protected BioNetwork network;
    protected Fingerprint fingerprint;

    @SuppressWarnings("deprecation")
    public void printInfo(CmdLineParser parser, String[] args) throws CmdLineException {

        super.printInfo(parser, args);

        if (this.entityType2Map < 1 || this.entityType2Map > 6) {
            throw new CmdLineException("Type of mapped entity must be between 1 and 6.");
        }
    }

    public Omics analyse(CmdLineParser parser, String[] args) throws IOException {

        this.fingerprint = (Fingerprint) super.analyse(parser, args);
        this.network = (new JSBML2Bionetwork(this.sbml)).getBioNetwork();
        Mapping map = new Mapping(logContent, network, fingerprint.getList_entities(), this.tab_inchiLayers,
                this.nameMapping, this.weightPrecision, this.outFileMapping, this.galaxyFile, this.entityType2Map);
        this.logContent = map.getLogContent();
        return map;
    }

    public static void main(String[] args) {
        exec(new Launcher_mapping(), args);
    }
}
