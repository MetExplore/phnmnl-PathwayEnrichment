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

package fr.inra.toulouse.metexplore.launcher;

import java.io.IOException;

import fr.inra.toulouse.metexplore.io.WritingComportment;
import fr.inra.toulouse.metexplore.omics.Omics;
import fr.inra.toulouse.metexplore.omics.PathwayEnrichment;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.kohsuke.args4j.CmdLineException;

public class Launcher_PathEnr extends Launcher_mapping implements WritingComportment{

    @Option(name = "-o3", aliases = "-outPath", usage = "Output file name for pathway enrichment result (by default: pathwayEnrichment.tsv).")
    protected String outFilePathEnr = "pathwayEnrichment.tsv";

    @Option(name = "-tEnr", aliases = "-typeEnr", usage = "1 for metabolites, 2 for reactions, 3 for pathway, 4 for enzyme, 5 for protein, 6 for gene (by default: pathways).")
    protected int entityType2Enrich = 3;

    @SuppressWarnings("deprecation")
    public void testParameters(String[] args) throws CmdLineException2 {

        if (this.entityType2Enrich < 1 || this.entityType2Enrich > 6) {
            throw new CmdLineException2(warn_type[0] + "enriched" + warn_type[1]);
        }

        super.testParameters(args);
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