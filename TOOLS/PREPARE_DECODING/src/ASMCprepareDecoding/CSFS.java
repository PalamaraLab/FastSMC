//    This file is part of ASMC, developed by Pier Francesco Palamara.
//
//    ASMC is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ASMC is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ASMC.  If not, see <https://www.gnu.org/licenses/>.


package ASMCprepareDecoding;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;

    public final class CSFS {

    ArraySpectrum arraySpectrum;
    // sequence emission
    TreeMap<Double, CSFSEntry> CSFS = new TreeMap<Double, CSFSEntry>();
    TreeMap<Double, CSFSEntry> foldedCSFS = new TreeMap<Double, CSFSEntry>();

    // ascertained emission
    TreeMap<Double, CSFSEntry> ascertainedCSFS = new TreeMap<Double, CSFSEntry>();
    TreeMap<Double, CSFSEntry> foldedAscertainedCSFS = new TreeMap<Double, CSFSEntry>();
    double[][] compressedAscertainedEmissionTable;

    double[] arraySamplingFactors;
    int samples;

    private CSFS(TreeMap<Double, CSFSEntry> CSFS) {
        this.CSFS = CSFS;
        if (CSFS.size() > 0) {
            this.foldedCSFS = this.foldCSFS(CSFS);
        }
    }

    public static CSFS loadFromFile(String fileName) throws IOException {
        BufferedReader br = Utils.openFile(fileName);
        String line;
        TreeMap<Double, CSFSEntry> parsedCSFS = new TreeMap<Double, CSFSEntry>();
        while ((line = br.readLine()) != null) {
            String[] splitString = line.split("\\s+");
            if (splitString[0].compareToIgnoreCase("Time:") == 0) {
                // this is a new CSFS line entry. Expect to read time vector, size vector, Mu, Samples, Intervals, then 3 CSFS lines.
                // parse time line
                ArrayList<Double> timeVector = new ArrayList<Double>();
                for (int i = 1; i < splitString.length; i++) {
                    double thisTime = Double.parseDouble(splitString[i]);
                    timeVector.add(thisTime);
                }
                // parse size line
                line = br.readLine();
                splitString = line.split("\\s+");
                if (splitString[0].compareToIgnoreCase("Size:") != 0) {
                    Utils.exit("Parsed line \"" + line + "\" after Time line in " + fileName + ".");
                }
                ArrayList<Double> sizeVector = new ArrayList<Double>();
                for (int i = 1; i < splitString.length; i++) {
                    double thisSize = Double.parseDouble(splitString[i]);
                    sizeVector.add(thisSize);
                }
                // parse mu
                line = br.readLine();
                splitString = line.split("\\s+");
                if (splitString[0].compareToIgnoreCase("Mu:") != 0) {
                    Utils.exit("Parsed line \"" + line + "\" after Size line in " + fileName + ".");
                }
                double mu = Double.parseDouble(splitString[1]);
                // parse samples
                line = br.readLine();
                splitString = line.split("\\s+");
                if (splitString[0].compareToIgnoreCase("Samples:") != 0) {
                    Utils.exit("Parsed line \"" + line + "\" after Mu line in " + fileName + ".");
                }
                int samples = Integer.parseInt(splitString[1]);
                // parse samples
                line = br.readLine();
                splitString = line.split("\\s+");
                if (splitString[0].compareToIgnoreCase("Interval:") != 0) {
                    Utils.exit("Parsed line \"" + line + "\" after Samples line in " + fileName + ".");
                }
                double from = Double.parseDouble(splitString[1]);
                double to = Double.parseDouble(splitString[2]);
                // parse CSFS
                double[][] CSFS = new double[3][samples - 1];
                for (int dist = 0; dist < 3; dist++) {
                    line = br.readLine();
                    splitString = line.split("\\s+");
                    for (int undist = 0; undist < splitString.length; undist++) {
                        double thisEntry = Double.parseDouble(splitString[undist]);
                        CSFS[dist][undist] = thisEntry;
                    }
                }
                CSFSEntry thisEntry = new CSFSEntry(timeVector, sizeVector, mu, from, to, samples, CSFS);
                parsedCSFS.put(from, thisEntry);
            } else {
                Utils.exit("Badly formatted CSFS file. Parsed line \"" + line + "\", which does not contain a Time definition. ");
            }
        }
        Utils.print("Read " + parsedCSFS.size() + " CSFS entries.");
        return new CSFS(parsedCSFS);
    }

    public boolean verify(ArrayList<Double> timeVectorOriginal, ArrayList<Double> sizeVectorOriginal, double mu, int samples, ArrayList<Double> discretizationOriginal) {
        try {
            ArrayList<Double> timeVector = (ArrayList<Double>) timeVectorOriginal.clone();
            timeVector.remove(timeVector.size() - 1);
            ArrayList<Double> sizeVector = (ArrayList<Double>) sizeVectorOriginal.clone();
            sizeVector.remove(sizeVector.size() - 1);
            ArrayList<Double> discretization = (ArrayList<Double>) discretizationOriginal.clone();
            discretization.remove(discretization.size() - 1);
            for (double from : discretization) {
                if (!CSFS.containsKey(from)) {
                    Utils.warning("CSFS does not contain interval " + from + ".");
                    return false;
                }
                CSFSEntry thisEntry = CSFS.get(from);
                if (thisEntry.mu != mu) {
                    Utils.warning("CSFS entry " + from + " has different mu: " + thisEntry.mu + ".");
                    return false;
                }
                if (!compareArrays(thisEntry.timeVector, timeVector)) {
                    Utils.print(thisEntry.timeVector);
                    Utils.print(timeVector);
                    Utils.warning("CSFS entry " + from + " has different time vector.");
                    return false;
                }
                if (!compareArrays(thisEntry.sizeVector, sizeVector)) {
                    Utils.warning("CSFS entry " + from + " has different size vector.");
                    return false;
                }
                if (thisEntry.samples != samples) {
                    if (samples == Integer.MAX_VALUE) {
                        samples = thisEntry.samples;
                    } else {
                        Utils.warning("CSFS entry " + from + " has different samples (want: " + samples + ", found: " + thisEntry.samples + ")");
                        return false;
                    }
                }
            }
            return true;
        } catch (Exception e) {
            Utils.warning("Something went wrong.");
            return false;
        }

    }

    private static boolean compareArrays(ArrayList<Double> a1, ArrayList<Double> a2) {
        if (a1.size() != a2.size()) {
            return false;
        }
        for (int i = 0; i < a1.size(); i++) {
            if (!a1.get(i).equals(a2.get(i))) {
                return false;
            }
        }
        return true;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (double from : CSFS.keySet()) {
            sb.append(CSFS.get(from));
        }
        return sb.toString();
    }

    public void fixAscertainment(Data data, int samples, Transition transition) {
        this.computeArraySamplingFactors(data, samples, transition);
        // CSFS is loaded here, but fixed later.
        for (double d : CSFS.keySet()) {
            ascertainedCSFS.put(d, (CSFSEntry) new CSFSEntry(CSFS.get(d)));
        }
        this.applyFactors();
        this.foldedAscertainedCSFS = this.foldCSFS(ascertainedCSFS);
        this.compressedAscertainedEmissionTable = compressCSFS(foldedAscertainedCSFS);
    }

    public static double[][] computeClassicEmission(double[] expectedTimes, double mu) {
        double[][] emission = new double[2][expectedTimes.length];
        for (int timeInterval = 0; timeInterval < expectedTimes.length; timeInterval++) {
            emission[0][timeInterval] = Math.exp(-2 * expectedTimes[timeInterval] * mu);
            emission[1][timeInterval] = 1 - emission[0][timeInterval];
        }
        return emission;
    }

    private void computeArraySamplingFactors(Data data, int samples, Transition transition) {
        this.samples = samples;
        double[] coalDist = transition.getCoalDist();
        double[] AFS = new double[samples];
        // the first entry of the CSFS may not be zero, since it's a shared doubleton
        int counter = 0;
        for (double from : CSFS.keySet()) {
            double coalescentProbabilityThisTimeSlice = coalDist[counter];
            CSFSEntry thisCsfsDM = CSFS.get(from);
            for (int row = 0; row < 3; row++) {
                for (int column = 0; column < samples - 1; column++) {
                    int pos = row + column;
                    if (pos > samples / 2) {
                        pos = samples - pos;
                    }
                    AFS[pos] += coalescentProbabilityThisTimeSlice * thisCsfsDM.CSFS[row][column];
                }
            }
            counter++;
        }

        // normalize spectrum
        AFS[0] = 0.;
        double norm = 0.;
        for (int i = 0; i < samples; i++) {
            norm += AFS[i];
        }
        for (int i = 0; i < samples; i++) {
            AFS[i] /= norm;
        }

        // fold AFS
        int halfTotal = samples / 2;
        for (int i = halfTotal + 1; i < samples; i++) {
            AFS[samples - i] += AFS[i];
            AFS[i] = 0;
        }
        // normalize spectrum
        norm = 0.;
        for (int i = 0; i < samples; i++) {
            norm += AFS[i];
        }
        for (int i = 0; i < samples; i++) {
            AFS[i] /= norm;
        }

        // foldedAFS contains probability a site has MAF i given the site is polymorphic in the sequence data
        double[] foldedAFS = new double[halfTotal + 1];
        for (int i = 0; i <= halfTotal; i++) {
            foldedAFS[i] = AFS[i];
        }

        // now get foldedAFS_array, the probability a site has MAF i given it is polymorphic in the sample (array)
        this.arraySpectrum = new ArraySpectrum(data, samples);
        double[] foldedAFS_array = arraySpectrum.spectrum;
        for (int i = 0; i < foldedAFS_array.length; i++) {
        }
        double[] samplingFactors = new double[halfTotal + 1];

        for (int i = 1; i < foldedAFS_array.length; i++) {
            samplingFactors[i] = foldedAFS_array[i] / foldedAFS[i];
        }
        arraySamplingFactors = samplingFactors;
    }

    private void applyFactors() {
        // apply sampling factors and renormalize
        // note that the first entry of the CSFS may not be zero, since it's a shared doubleton
        for (double from : ascertainedCSFS.keySet()) {
            double[][] thisCSFS = ascertainedCSFS.get(from).CSFS;
            thisCSFS[0][0] = 0.;
            double norm = 0.;
            for (int row = 0; row < 3; row++) {
                for (int column = 0; column < samples - 1; column++) {
                    // if the spectrum is folded, this emission is mapped to this position
                    int pos = row + column;
                    if (pos > samples / 2) {
                        pos = samples - pos;
                    }
                    // and if we're looking at array data, this MAF is adjusted using this factor
                    double factor = arraySamplingFactors[pos];
                    // apply factor
                    thisCSFS[row][column] *= factor;
                    // sum value to renomralize to 1 later on
                    norm += thisCSFS[row][column];
                }
            }
            norm /= 1 - arraySpectrum.monomorphic;
            for (int row = 0; row < 3; row++) {
                for (int column = 0; column < samples - 1; column++) {
                    thisCSFS[row][column] /= norm;
                }
            }
            thisCSFS[0][0] = arraySpectrum.monomorphic;
            ascertainedCSFS.get(from).CSFS = thisCSFS;
        }
    }

    public TreeMap<Double, CSFSEntry> foldCSFS(TreeMap<Double, CSFSEntry> CSFS) {
        TreeMap<Double, CSFSEntry> foldedCSFS = new TreeMap<Double, CSFSEntry>();
        int samples = CSFS.firstEntry().getValue().samples;
        int undistinguished = samples - 2;
        for (double from : CSFS.keySet()) {
            CSFSEntry foldedEntry = new CSFSEntry(CSFS.get(from));
            double[][] thisCsfs_double = foldedEntry.CSFS;
            // code to fold the spectrum
            if (samples % 2 != 0) {
                Utils.exit("ConditionalSFS called with odd number of samples.");
            }
            int half = samples / 2;
            double[][] thisCsfs_double_folded = new double[2][half + 1];
            for (int row = 0; row < 3; row++) {
                for (int column = 0; column < undistinguished + 1; column++) {
                    Integer[] coord = new Integer[]{row, column};
                    Integer[] foldedCoord = getFoldedObservationFromUnfolded(coord, samples);
                    thisCsfs_double_folded[foldedCoord[0]][foldedCoord[1]] += thisCsfs_double[row][column];
                }
            }
            foldedEntry.CSFS = thisCsfs_double_folded;
            foldedCSFS.put(from, foldedEntry);
        }
        return foldedCSFS;
    }

    private static Integer[] getFoldedObservationFromUnfolded(Integer[] unfolded, int totalSamples) {
        int dist = unfolded[0];
        int undist = unfolded[1];
        if (totalSamples % 2 != 0) {
            Utils.exit("Function getFoldedObservationFromUnfolded was called with odd total sample size. Only diploid samples are supported at the moment.");
        }
        int halfTotal = totalSamples / 2;
        if (undist + dist > halfTotal) {
            // flip
            undist = (totalSamples - 2 - undist);
        }
        if (dist == 2) {
            dist = 0;
        }
        Integer[] ret = new Integer[]{dist, undist};
        return ret;
    }

    public double[][] compressCSFS(TreeMap<Double, CSFSEntry> CSFS) {
        double[][] compressed = new double[2][CSFS.size()];
        int timeInterval = 0;
        for (double from : CSFS.keySet()) {
            double[][] thisCSFS = CSFS.get(from).CSFS;
            for (int k = 0; k < thisCSFS[0].length; k++) {
                compressed[0][timeInterval] += thisCSFS[0][k];
                compressed[1][timeInterval] += thisCSFS[1][k];
            }
            timeInterval++;
        }
        return compressed;
    }

//    public static void unitTest() throws IOException {
//        Data data = new Data("CHR22_100_1/out.1234.CEU.full.n100.chr22.len1.array");
//        double[] discretization = {0.0, 1118.2, 1472.2, 1849.7, 2497.0, 3963.8, 9120.8, 15832.9, 24139.9, 34891.6, Double.POSITIVE_INFINITY};
//        Transition transition = new Transition(Transition.EUtime_array, Transition.EUsize_array, discretization, Transition.TransitionType.CSC);
//        CSFS csfs = loadFromFile("CHR22_100_1/out.1234.CEU.full.n100.chr22.len1.array.csfs");
//        csfs.computeArraySamplingFactors(data, 100, transition);
//        csfs.applyFactors();
//        csfs.compressedAscertainedEmissionTable = csfs.compressCSFS(csfs.CSFS);
//        Utils.printMatrixOfDoubles(csfs.compressedAscertainedEmissionTable);
//    }
}

class CSFSEntry {

    ArrayList<Double> timeVector;
    ArrayList<Double> sizeVector;
    double mu;
    double from;
    double to;
    int samples;
    double[][] CSFS;

    // copy constructor
    CSFSEntry(CSFSEntry other) {
        this.timeVector = new ArrayList<Double>();
        for (int i = 0; i < other.timeVector.size(); i++) {
            this.timeVector.add(other.timeVector.get(i));
        }
        this.sizeVector = new ArrayList<Double>();
        for (int i = 0; i < other.sizeVector.size(); i++) {
            this.sizeVector.add(other.sizeVector.get(i));
        }
        this.mu = other.mu;
        this.from = other.from;
        this.to = other.to;
        this.samples = other.samples;
        this.CSFS = new double[other.CSFS.length][other.CSFS[0].length];
        for (int i = 0; i < other.CSFS.length; i++) {
            for (int j = 0; j < other.CSFS[0].length; j++) {
                this.CSFS[i][j] = other.CSFS[i][j];
            }
        }
    }

    // constructor
    CSFSEntry(ArrayList<Double> timeVector, ArrayList<Double> sizeVector,
            double mu, double from, double to, int samples, double[][] CSFS) {
        if (timeVector.size() != sizeVector.size() || CSFS.length != 3 || CSFS[0].length != (samples - 1) || from >= to) {
            Utils.print(timeVector.size());
            Utils.print(sizeVector.size());
            Utils.print(CSFS.length);
            Utils.print(CSFS[0].length);
            Utils.print(from);
            Utils.print(to);
            Utils.exit("Malformed CSFS entry.");
        }
        this.timeVector = timeVector;
        this.sizeVector = sizeVector;
        this.mu = mu;
        this.samples = samples;
        this.from = from;
        this.to = to;
        this.CSFS = CSFS;
    }

    // used to write to file
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Time:\t");
        for (int i = 0; i < timeVector.size(); i++) {
            if (i > 0) {
                sb.append(" ");
            }
            sb.append(timeVector.get(i));
        }
        sb.append("\n");
        sb.append("Size:\t");
        for (int i = 0; i < sizeVector.size(); i++) {
            if (i > 0) {
                sb.append(" ");
            }
            sb.append(sizeVector.get(i));
        }
        sb.append("\n");
        sb.append("Mu:\t").append(mu).append("\n");
        sb.append("Samples:\t").append(samples).append("\n");
        sb.append("Interval:\t").append(from).append(" ").append(to).append("\n");
        for (int dist = 0; dist < CSFS.length; dist++) {
            for (int undist = 0; undist < CSFS[0].length; undist++) {
                if (undist > 0) {
                    sb.append(" ");
                }
                sb.append(CSFS[dist][undist]);
            }
            sb.append("\n");
        }
        return sb.toString();
    }

}
