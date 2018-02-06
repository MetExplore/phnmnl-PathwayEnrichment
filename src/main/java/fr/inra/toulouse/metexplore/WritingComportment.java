package fr.inra.toulouse.metexplore;

public class WritingComportment {

    public String removeSciNot(double value) {
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

    public String round(double value) {
        return String.valueOf((double) Math.round(value * 100) / 100);
    }
}
