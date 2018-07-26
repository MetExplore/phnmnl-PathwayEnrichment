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

import fr.inra.toulouse.metexplore.io.WritingComportment;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.File;
import java.io.IOException;

import static java.lang.System.exit;

public abstract class Launcher implements WritingComportment{
    @Option(name = "-h", aliases = "-help", usage = "Prints this help.")
    protected boolean phelp = false;

    @Option(name = "-v", aliases = "-version", usage = "Prints the current version of the program.")
    protected boolean version;

    @Option(name = "-gal", aliases = "-galaxy", usage = "For galaxy compliance: formatting pathway output and creating a new one containing log information.")
    protected String galaxyFile = "";

    protected static File logFile;
    protected String logContent;

    protected static long startTime;

    public static void timeCalculation(long elapsedTime) {
        long min = elapsedTime / 60000000000L;
        long sec = elapsedTime / 1000000000L - (min * 60L);
        System.out.println("Time to run the process : " + min + "min " + sec + "s");

    }

    public void printError(CmdLineParser parser, CmdLineException e){
        writeLog("[FATAL] " + e.getMessage());
        System.err.println("Options:");
        parser.printUsage(System.err);
        sysExit(this.logContent,"", this.galaxyFile,1);
    }

    public void printInfo(CmdLineParser parser){

        if (this.phelp) {
            System.out.println("Options:");
            parser.printUsage(System.out);
            exit(0);
        }

        if (this.version) {
            System.out.println("Version: 2.0");
            exit(0);
        }
    }

    public void writeLog(String warning){
        this.logContent = writeLog(this.logContent,warning);
    }
}