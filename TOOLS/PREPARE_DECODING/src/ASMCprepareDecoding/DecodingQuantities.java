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

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;
import static ASMCprepareDecoding.CSFS.loadFromFile;
import org.jblas.DoubleMatrix;

public class DecodingQuantities {

    // these are the quantities loaded from file, which are required for the linear-time decoding
    Transition.TransitionType transitionType;
    DoubleMatrix initialStateProb;
    double[] timeVector;
    double[] sizeVector;
    double[] discretization;
    double[] expectedTimes;
    ArrayList<Double> geneticDistancesList = new ArrayList<Double>();
    ArrayList<Integer> physDistancesList = new ArrayList<Integer>();
    ArrayList<DoubleMatrix> DvectorMap = new ArrayList<DoubleMatrix>();
    ArrayList<DoubleMatrix> BvectorMap = new ArrayList<DoubleMatrix>();
    ArrayList<DoubleMatrix> UvectorMap = new ArrayList<DoubleMatrix>();
    ArrayList<DoubleMatrix> rowRatiosVectorMap = new ArrayList<DoubleMatrix>();
    ArrayList<DoubleMatrix> HomozygousEmissionMap = new ArrayList<DoubleMatrix>();
    DoubleMatrix columnRatios;

    private int states, CSFSSamples;
    double mu;

    // sequence emission
    TreeMap<Double, CSFSEntry> CSFSmap;
    TreeMap<Double, CSFSEntry> foldedCSFSmap;
    private double[][] classicEmissionTable;

    // ascertained emission
    TreeMap<Double, CSFSEntry> ascertainedCSFSmap;
    TreeMap<Double, CSFSEntry> foldedAscertainedCSFSap;
    private double[][] compressedEmissionTable;

    // for a number between 10^(n) and 10^(n+1), proceed in steps of size max(1, 10^(n-precision))
    final static int precision = 2;
    final static double minGenetic = 1e-10;

    // start/end for physical and genetic distances (in Morgans)
    private static final double startGen = 1e-10;
    private static final double maxGen = 0.3; // 30 centiMorgans
    private static final int startPhys = 1;
    private static final int maxPhys = 100000000; // 100 Mb

    // for each genetic distance, compute decoding quantities
    public DecodingQuantities(CSFS csfs, Transition transition, double mu) throws IOException {
        this.timeVector = transition.timeVector;
        this.mu = mu;
        this.sizeVector = transition.sizeVector;
        this.CSFSSamples = csfs.samples;
        this.discretization = transition.discretization;
        this.expectedTimes = transition.expectedTimes;
        this.states = transition.discretization.length - 1;
        int lastPercentage = -1;
        columnRatios = transition.columnRatios;

        // compute transition quantities
        // putting 0 in. This is not really used, depending on the value of smallNumber in the decoding side.
        computeTransitionQuantitiesForOnePosition(0.0, transition);
        geneticDistancesList.add(0.0);
        double genDist = startGen;
        while (genDist < maxGen) {
            geneticDistancesList.add(genDist);
            genDist = DecodingQuantities.nextGen(genDist);
        }
        // start at 1, 0 is done already
        for (int i = 1; i < geneticDistancesList.size(); i++) {
            genDist = geneticDistancesList.get(i);
            // compute emission and transition quantities in linear time
            computeTransitionQuantitiesForOnePosition(genDist, transition);
            int percentage = (int) Math.round(100 * i / (float) geneticDistancesList.size());
            if (percentage != lastPercentage) {
                Utils.printf("\rGenetic distances progress: " + percentage + "%");
            }
            lastPercentage = percentage;
        }
        Utils.print("");

        // compute homozygous emissions
        int phys = startPhys;
        while (phys < maxPhys) {
            this.physDistancesList.add(phys);
            phys = DecodingQuantities.nextPhys(phys);
        }
        for (int i = 0; i < physDistancesList.size(); i++) {
            HomozygousEmissionMap.add(new DoubleMatrix(CSFS.computeClassicEmission(this.expectedTimes, physDistancesList.get(i) * this.mu)));
            int percentage = (int) Math.round(100 * i / (float) physDistancesList.size());
            if (percentage != lastPercentage) {
                Utils.printf("\rPhysical distances progress: " + percentage + "%");
            }
            lastPercentage = percentage;
        }
        Utils.print("");

        // get initial state distribution from coalescent distribution
        double lastProb = 0;
        DoubleMatrix cumProb = new DoubleMatrix(this.states);
        initialStateProb = new DoubleMatrix(this.states);
        for (int i = 0; i < this.states; i++) {
            double timeT = transition.discretization[i + 1];
            cumProb.put(i, transition.cumulativeCoalesceFromStoT(0, timeT));
            initialStateProb.put(i, cumProb.get(i) - lastProb);
            lastProb = cumProb.get(i);
        }

        this.transitionType = transition.type;

        this.CSFSmap = csfs.CSFS;
        this.foldedCSFSmap = csfs.foldedCSFS;
        this.classicEmissionTable = CSFS.computeClassicEmission(this.expectedTimes, mu);

        this.ascertainedCSFSmap = csfs.ascertainedCSFS;
        this.foldedAscertainedCSFSap = csfs.foldedAscertainedCSFS;
        this.compressedEmissionTable = csfs.compressedAscertainedEmissionTable;
    }

    public static int nextPhys(int phys) {
        if (phys < 0) {
            Utils.exit("Int overflow: " + phys);
        }
        int log10 = (int) Math.max(0, Math.floor(Math.log10(phys)) - precision);
        int factor = (int) Math.pow(10, log10);
        int next = Math.round(phys / (float) factor + 1) * factor;
        return next;
    }

    public static int roundPhysical(int phys) {
        if (phys < 0) {
            Utils.exit("Int overflow: " + phys);
        }
        int log10 = (int) Math.max(0, Math.floor(Math.log10(phys)) - precision);
        int factor = (int) Math.pow(10, log10);
        int rounded = Math.round(phys / (float) factor) * factor;
        return rounded;
    }

    public static double nextGen(double gen) {
        double gene1e10 = gen * 1e10;
        int log10 = (int) Math.max(0, Math.floor(Math.log10(gene1e10)) - precision);
        double factor = Math.pow(10, log10);
        double next = (Math.round(gene1e10 / factor) + 1) * factor;
        return next / 1e10;
    }

    public static double roundMorgans(double gen) {
        double gene1e10 = gen * 1e10;
        int log10 = (int) Math.max(0, Math.floor(Math.log10(gene1e10)) - precision);
        double factor = Math.pow(10, log10);
        double rounded = Math.round(gene1e10 / factor) * factor;
        return Math.max(minGenetic, rounded / 1e10);
    }

    // compute quantities for a give genetic distance
    private void computeTransitionQuantitiesForOnePosition(double genDist, Transition transition) {
        ArrayList<DoubleMatrix> res = transition.getLinearTimeDecodingQuantitiesAndMatrixGivenDistance(genDist);
        // index 0 is transition matrix
        DoubleMatrix D = res.get(0);
        DoubleMatrix B = res.get(1);
        DoubleMatrix U = res.get(2);
        DoubleMatrix RR = res.get(3);
        DvectorMap.add(D);
        BvectorMap.add(B);
        UvectorMap.add(U);
        rowRatiosVectorMap.add(RR);
    }

    // save all compute quantities to file
    public void saveDecodingQuantities(String outputFileRoot) throws IOException {
        BufferedWriter outBuff = Utils.openGzipFileForWriting(outputFileRoot + ".decodingQuantities.gz");
        // write model parameters
        outBuff.write("TransitionType\n" + this.transitionType + "\n\n");
        outBuff.write("States\n" + this.states + "\n\n");
        outBuff.write("CSFSSamples\n" + this.CSFSSamples + "\n\n");
        outBuff.write("TimeVector\n");
        outBuff.write(Utils.doubleMatrixToString(timeVector));
        outBuff.write("\n");
        outBuff.write("SizeVector\n");
        outBuff.write(Utils.doubleMatrixToString(sizeVector));
        outBuff.write("\n");
        outBuff.write("Discretization\n");
        outBuff.write(Utils.doubleMatrixToString(discretization));
        outBuff.write("\n");
        outBuff.write("ExpectedTimes\n");
        outBuff.write(Utils.doubleMatrixToString(expectedTimes));
        // write sequence Emissions
        outBuff.write("\n");
        for (int undistinguished = 0; undistinguished < this.CSFSmap.firstEntry().getValue().CSFS[0].length; undistinguished++) {
            outBuff.write("CSFS\t" + undistinguished + "\n");
            for (int distinguished = 0; distinguished < 3; distinguished++) {
                for (double from : this.CSFSmap.keySet()) {
                    outBuff.write(this.CSFSmap.get(from).CSFS[distinguished][undistinguished] + "\t");
                }
                outBuff.write("\n");
            }
        }
        outBuff.write("\n");
        for (int undistinguished = 0; undistinguished < this.foldedCSFSmap.firstEntry().getValue().CSFS[0].length; undistinguished++) {
            outBuff.write("FoldedCSFS\t" + undistinguished + "\n");
            for (int distinguished = 0; distinguished < 2; distinguished++) {
                for (double from : this.foldedCSFSmap.keySet()) {
                    outBuff.write(this.foldedCSFSmap.get(from).CSFS[distinguished][undistinguished] + "\t");
                }
                outBuff.write("\n");
            }
        }
        outBuff.write("\n");
        outBuff.write("ClassicEmission\n");
        outBuff.write(Utils.doubleMatrixToString(classicEmissionTable));
        // write ascertained Emissions
        outBuff.write("\n");
        for (int undistinguished = 0; undistinguished < this.ascertainedCSFSmap.firstEntry().getValue().CSFS[0].length; undistinguished++) {
            outBuff.write("AscertainedCSFS\t" + undistinguished + "\n");
            for (int distinguished = 0; distinguished < 3; distinguished++) {
                for (double from : this.ascertainedCSFSmap.keySet()) {
                    outBuff.write(this.ascertainedCSFSmap.get(from).CSFS[distinguished][undistinguished] + "\t");
                }
                outBuff.write("\n");
            }
        }
        outBuff.write("\n");
        for (int undistinguished = 0; undistinguished < this.foldedAscertainedCSFSap.firstEntry().getValue().CSFS[0].length; undistinguished++) {
            outBuff.write("FoldedAscertainedCSFS\t" + undistinguished + "\n");
            for (int distinguished = 0; distinguished < 2; distinguished++) {
                for (double from : this.foldedAscertainedCSFSap.keySet()) {
                    outBuff.write(this.foldedAscertainedCSFSap.get(from).CSFS[distinguished][undistinguished] + "\t");
                }
                outBuff.write("\n");
            }
        }
        outBuff.write("\n");
        outBuff.write("CompressedAscertainedEmission\n");
        outBuff.write(Utils.doubleMatrixToString(compressedEmissionTable));
        // write initial state distribution
        outBuff.write("\n");
        outBuff.write("initialStateProb\n");
        outBuff.write(Utils.doubleMatrixToString(initialStateProb.transpose()));
        // write column ratios
        outBuff.write("\n");
        outBuff.write("ColumnRatios\n");
        outBuff.write(Utils.doubleMatrixToString(columnRatios.transpose()));
        // write row ratios
        outBuff.write("\n");
        outBuff.write("RowRatios\n");
        for (int i = 0; i < rowRatiosVectorMap.size(); i++) {
            DoubleMatrix thisRR = rowRatiosVectorMap.get(i);
            outBuff.write(geneticDistancesList.get(i) + "\t" + Utils.doubleMatrixToString(thisRR.transpose()));
        }
        // write U vectors
        outBuff.write("\n");
        outBuff.write("Uvectors\n");
        for (int i = 0; i < UvectorMap.size(); i++) {
            DoubleMatrix thisU = UvectorMap.get(i);
            outBuff.write(geneticDistancesList.get(i) + "\t" + Utils.doubleMatrixToString(thisU.transpose()));
        }
        // write B vectors
        outBuff.write("\n");
        outBuff.write("Bvectors\n");
        for (int i = 0; i < BvectorMap.size(); i++) {
            DoubleMatrix thisB = BvectorMap.get(i);
            outBuff.write(geneticDistancesList.get(i) + "\t" + Utils.doubleMatrixToString(thisB.transpose()));
        }
        // write D vectors
        outBuff.write("\n");
        outBuff.write("Dvectors\n");
        for (int i = 0; i < DvectorMap.size(); i++) {
            DoubleMatrix thisD = DvectorMap.get(i);
            outBuff.write(geneticDistancesList.get(i) + "\t" + Utils.doubleMatrixToString(thisD.transpose()));
        }
        // write homozygous emissions
        outBuff.write("\n");
        outBuff.write("HomozygousEmissions\n");
        for (int i = 0; i < physDistancesList.size(); i++) {
            int phys = physDistancesList.get(i);
            outBuff.write(phys + "\t" + Utils.doubleMatrixToString(HomozygousEmissionMap.get(i).getRow(0)));
        }
        outBuff.close();
    }

//    // unit test
//    public static void unitTest() throws IOException {
//        Data data = new Data("CHR22_100_1/out.1234.CEU.full.n100.chr22.len1.array");
//        double[] discretization = {0.0, 1118.2, 1472.2, 1849.7, 2497.0, 3963.8, 9120.8, 15832.9, 24139.9, 34891.6, Double.POSITIVE_INFINITY};
//        Transition transition = new Transition(Transition.EUtime_array, Transition.EUsize_array, discretization, Transition.TransitionType.CSC);
//        CSFS csfs = loadFromFile("CHR22_100_1/out.1234.CEU.full.n100.chr22.len1.array.csfs");
//        csfs.fixAscertainment(data, 100, transition);
//        double mutRate = 1.65E-8;
//        DecodingQuantities decodingQuantities = new DecodingQuantities(csfs, transition, mutRate);
//        decodingQuantities.saveDecodingQuantities("CHR22_100_1/out.1234.CEU.full.n100.chr22.len1.array.temp");
//    }

}
