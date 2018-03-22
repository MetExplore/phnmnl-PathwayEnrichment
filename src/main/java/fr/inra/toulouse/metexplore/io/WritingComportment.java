package fr.inra.toulouse.metexplore.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

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
        System.out.println(message.replaceAll("\n$", ""));
        log += message;
        return log;
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
