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

package fr.inra.toulouse.metexplore.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import static java.lang.System.exit;

public interface WritingComportment {

    default File createFile (String galaxyFile) {
        File log = null;
        if (galaxyFile != "") {
            log = new File(galaxyFile);
            if (log.isFile()) {
                //write a new file if already exists
                log.delete();
            }
            try {
                log.createNewFile();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return log;
    }

    default void writeOutput(String text4outputFileInfo, File galaxyFile) throws IOException {
        if (!galaxyFile.equals("")) {//if "writing console log in a file" functionality is activated
            BufferedWriter b = new BufferedWriter(new FileWriter(galaxyFile, true));
            b.write(text4outputFileInfo);
            b.close();
        }
    }

    default String writeLog(String log, String message) {
        //Avoid to skip a line
        if (message.contains("[WARNING]") || message.contains("[FATAL]")) {
            System.err.println(message.replaceAll("\n$", ""));
        }else{
            System.out.println(message.replaceAll("\n$", ""));
        }
        log += message;
        return log;
    }

    default void sysExit(String log, String message, String galaxy, int exitCode){

        writeLog(log,message);
        try{
            if(!galaxy.isEmpty()) {
                File galaxyFile = new File(galaxy);
                writeOutput(log + message, galaxyFile);
            }
        }catch (IOException e) {
            e.printStackTrace();
        }
        exit(exitCode);
    }

    default String removeSciNot(double value) {
        String tmp = (String.valueOf(value));
        if (value < 1e-3) {
            String[] splitting = tmp.split("E");
            String power = splitting[1];
            if (power.charAt(1) == '0') {//remove the O after the E if exists
                power = power.substring(2, power.length());
            } else power = power.substring(1, power.length());
            String[] number = splitting[0].split("\\.");//obtain the integer and the decimal parts
            return "0." + (new String(new char[Integer.parseInt(power) - 1]).replace("\0", "0")) + number[0] + number[1];
        }
        return tmp;
    }

    default String round(double value) {
        return String.valueOf((double) Math.round(value * 100) / 100);
    }

    default String calculPercent(int numerator, int denominator) {
        return round(((double) numerator / (double) denominator) * 100);
    }
}
