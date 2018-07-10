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

import java.util.ArrayList;
import java.util.Arrays;
import org.jblas.DoubleMatrix;
import static org.jblas.MatrixFunctions.expm;
import org.jblas.util.Logger;

public class Transition {

    enum TransitionType {

        SMC, SMC1, CSC
    }

    double[] timeVector;
    double[] sizeVector;
    double[] discretization;
    double[] timeVectorPlusInfinity;
    double[] expectedTimes;
    TransitionType type;
    public static int states;

    // coalescent arrays that will be computed once and depend on the demography
    double[] probNotCoalesceBetweenExpectedTimes;
    double[] probNotCoalesceBetweenTimeIntervals;
    double[] probCoalesceBetweenExpectedTimesAndUpperLimit;
    DoubleMatrix columnRatios;

    private static final DoubleMatrix fourState_timeInf = new DoubleMatrix(new double[][]{{0., 0., 0., 1.}, {0., 0., 0., 1.}, {0., 0., 0., 1.}, {0., 0., 0., 1.}});
    private static final DoubleMatrix threeState_timeInf = new DoubleMatrix(new double[][]{{0., 0., 1.}, {0., 0., 1.}, {0., 0., 1.}});
    private static final DoubleMatrix fourState_identity = new DoubleMatrix(new double[][]{{1., 0., 0., 0.}, {0., 1., 0., 0.}, {0., 0., 1., 0.}, {0., 0., 0., 1.}});
    private static final DoubleMatrix threeState_identity = new DoubleMatrix(new double[][]{{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}});
    
    public static final double[] EUsize_array = new double[]{145041., 129827., 116209., 104020., 93109., 83342., 74600., 66775., 59771., 53501., 47892., 44915., 43684., 42486., 41321., 40188., 39086., 38014., 36972., 35958., 34972., 34013., 33080., 32173., 31291., 30433., 29598., 28787., 27997., 27230., 26483., 25757., 25050., 24364., 23695., 23046., 22414., 21799., 21201., 20620., 20055., 19505., 18970., 18450., 17944., 17452., 16973., 16508., 16055., 15615., 15186., 14770., 14365., 13971., 13588., 13215., 12853., 12501., 12158., 11824., 11500., 11185., 10878., 10580., 10289., 10007., 9733., 9466., 9206., 8954., 8708., 8470., 7871., 7605., 6890., 6080., 5191., 4225., 3180., 2843., 2604., 2830., 2940., 2694., 5002., 6138., 5562., 6093., 7825., 9864., 13076., 15823., 18475., 20949., 23350., 25032., 26366., 26381., 26035., 24927., 23896., 22998., 23071., 23312., 23339., 23145., 22510., 21536., 20189., 18401., 15099., 15628., 15628.};
    public static final double[] EUtime_array = new double[]{0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280., 290., 300., 310., 320., 330., 340., 350., 360., 370., 380., 390., 400., 410., 420., 430., 440., 450., 460., 470., 480., 490., 500., 510., 520., 530., 540., 550., 560., 570., 580., 590., 600., 610., 620., 630., 640., 650., 660., 670., 680., 690., 700., 710., 720., 724., 795., 874., 960., 1055., 1160., 1275., 1402., 1541., 1694., 1862., 2047., 2250., 2474., 2720., 2990., 3287., 3614., 3973., 4368., 4802., 5280., 5805., 6382., 7017., 7715., 8482., 9325., 10252., 11271., 12391., 13623., 14977., 16466., 18102., 19901., 21879., 24053., 26443., Double.POSITIVE_INFINITY};
    public static final ArrayList<Double> EUsize = new ArrayList<Double>(Arrays.asList(145041., 129827., 116209., 104020., 93109., 83342., 74600., 66775., 59771., 53501., 47892., 44915., 43684., 42486., 41321., 40188., 39086., 38014., 36972., 35958., 34972., 34013., 33080., 32173., 31291., 30433., 29598., 28787., 27997., 27230., 26483., 25757., 25050., 24364., 23695., 23046., 22414., 21799., 21201., 20620., 20055., 19505., 18970., 18450., 17944., 17452., 16973., 16508., 16055., 15615., 15186., 14770., 14365., 13971., 13588., 13215., 12853., 12501., 12158., 11824., 11500., 11185., 10878., 10580., 10289., 10007., 9733., 9466., 9206., 8954., 8708., 8470., 7871., 7605., 6890., 6080., 5191., 4225., 3180., 2843., 2604., 2830., 2940., 2694., 5002., 6138., 5562., 6093., 7825., 9864., 13076., 15823., 18475., 20949., 23350., 25032., 26366., 26381., 26035., 24927., 23896., 22998., 23071., 23312., 23339., 23145., 22510., 21536., 20189., 18401., 15099., 15628., 15628.));
    public static final ArrayList<Double> EUtime = new ArrayList<Double>(Arrays.asList(0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280., 290., 300., 310., 320., 330., 340., 350., 360., 370., 380., 390., 400., 410., 420., 430., 440., 450., 460., 470., 480., 490., 500., 510., 520., 530., 540., 550., 560., 570., 580., 590., 600., 610., 620., 630., 640., 650., 660., 670., 680., 690., 700., 710., 720., 724., 795., 874., 960., 1055., 1160., 1275., 1402., 1541., 1694., 1862., 2047., 2250., 2474., 2720., 2990., 3287., 3614., 3973., 4368., 4802., 5280., 5805., 6382., 7017., 7715., 8482., 9325., 10252., 11271., 12391., 13623., 14977., 16466., 18102., 19901., 21879., 24053., 26443., Double.POSITIVE_INFINITY));

    // unoptimized
    public static ArrayList<Double> getTimeExponentialQuantiles(int numQuantiles, ArrayList<Double> timeVector, ArrayList<Double> sizeFromVector) {
        double slice = 1. / numQuantiles;
        double nextQuant = slice;
        double timeStep = 0.1;
        ArrayList<Double> quantiles = new ArrayList<Double>(numQuantiles);
        quantiles.add(0.);
        double pNotCoal = 1.;
        for (int i = 0; i < timeVector.size() - 1; i++) {
            double fromTime = timeVector.get(i);
            double toTime = timeVector.get(i + 1);
            double size = sizeFromVector.get(i);
            double notCoalRate = 1 - timeStep / size;
            for (double t = fromTime; t < toTime; t += timeStep) {
                pNotCoal *= notCoalRate;
                if (1 - pNotCoal > nextQuant) {
                    nextQuant += slice;
                    double rounded = Math.round(t * 1000.) / 1000.;
                    quantiles.add(rounded);
                    if (nextQuant >= 1.0 - 1E-10) {
                        return quantiles;
                    }
                }
            }
        }
        return quantiles;
    }

    // unoptimized
    public static ArrayList<Double> getTimeErlangQuantiles(int numQuantiles, ArrayList<Double> timeVector, ArrayList<Double> sizeFromVector) {
        double slice = 1. / numQuantiles;
        double nextQuant = slice;
        double timeStep = 0.1;
        ArrayList<Double> quantiles = new ArrayList<Double>(numQuantiles);
        quantiles.add(0.);
        double normalizer = 0.;
        double pNotCoal = 1.;
        final double MAX_T = sizeFromVector.get(sizeFromVector.size() - 1) * 20; //  10 times the ancestral size
        for (int i = 0; i < timeVector.size() - 1; i++) {
            double fromTime = timeVector.get(i);
            double toTime = timeVector.get(i + 1);
            double size = sizeFromVector.get(i);
            double coalRate = timeStep / size;
            double notCoalRate = 1 - coalRate;
            for (double t = fromTime; t < toTime && t < MAX_T; t += timeStep) {
                pNotCoal *= notCoalRate;
                normalizer += t * coalRate * pNotCoal;
            }
        }
        pNotCoal = 1.;
        double normalizedCumulative = 0.;
        for (int i = 0; i < timeVector.size() - 1; i++) {
            double fromTime = timeVector.get(i);
            double toTime = timeVector.get(i + 1);
            double size = sizeFromVector.get(i);
            double coalRate = timeStep / size;
            double notCoalRate = 1 - coalRate;
            for (double t = fromTime; t < toTime && t < MAX_T; t += timeStep) {
                pNotCoal *= notCoalRate;
                normalizedCumulative += t * coalRate * pNotCoal / normalizer;
                if (normalizedCumulative >= nextQuant) {
                    nextQuant += slice;
                    double rounded = Math.round(t * 1000.) / 1000.;
                    quantiles.add(rounded);
                    if (nextQuant >= 1.0) {
                        return quantiles;
                    }
                }
            }
        }
        return quantiles;
    }

    private static int logLevel = 100;

    public Transition(double[] timeVector, double[] sizeVector, double[] discretization, TransitionType type) {
        Logger.getLogger().setLevel(logLevel);
        this.timeVector = timeVector;
        this.sizeVector = sizeVector;
        this.discretization = discretization;
        this.type = type;
        // adjusting for expectedTimeFromStoT to work. Will fix
        timeVectorPlusInfinity = new double[timeVector.length + 1];
        for (int i = 0; i < timeVector.length; i++) {
            timeVectorPlusInfinity[i] = timeVector[i];
        }
        timeVectorPlusInfinity[timeVectorPlusInfinity.length - 1] = Double.POSITIVE_INFINITY;
        expectedTimes = expectedIntervalTimesPiecewise();
        states = discretization.length - 1;
        computeCoalescentVectors();
    }

    public double[] getExpectedTimes() {
        return this.expectedTimes;
    }

    ArrayList<DoubleMatrix> getLinearTimeDecodingQuantitiesAndMatrixGivenDistance(double rho) {
        Pair<ArrayList<DoubleMatrix>, ArrayList<DoubleMatrix>> omegas = getOmegas(rho, type);
        ArrayList<DoubleMatrix> omegasAtBoundaries = omegas.getKey();
        ArrayList<DoubleMatrix> omegasAtExpectedTimes = omegas.getValue();
        DoubleMatrix D = new DoubleMatrix(states);
        DoubleMatrix B = new DoubleMatrix(states - 1);
        DoubleMatrix U = new DoubleMatrix(states - 1);
        DoubleMatrix RR = new DoubleMatrix(states - 1);

        // others should be computed only for i > 0
        DoubleMatrix omegaAtBoundaries = omegasAtBoundaries.get(0);
        DoubleMatrix omegaAtExpectedTimes = omegasAtExpectedTimes.get(0);
        double diagonal = (type == TransitionType.CSC)
                ? omegaAtExpectedTimes.get(0, 0) + probCoalesceBetweenExpectedTimesAndUpperLimit[0] * (omegaAtExpectedTimes.get(0, 1) + omegaAtExpectedTimes.get(0, 2)) + omegaAtExpectedTimes.get(0, 3) - omegaAtBoundaries.get(0, 3)
                : omegaAtExpectedTimes.get(0, 0) + probCoalesceBetweenExpectedTimesAndUpperLimit[0] * omegaAtExpectedTimes.get(0, 1) + omegaAtExpectedTimes.get(0, 2) - omegaAtBoundaries.get(0, 2);
        D.put(0, diagonal);

        // now compute all for each i
        for (int i = 1; i < states; i++) {
            omegaAtBoundaries = omegasAtBoundaries.get(i);
            omegaAtExpectedTimes = omegasAtExpectedTimes.get(i);
            diagonal = (type == TransitionType.CSC)
                    ? omegaAtExpectedTimes.get(0, 0) + probCoalesceBetweenExpectedTimesAndUpperLimit[i] * (omegaAtExpectedTimes.get(0, 1) + omegaAtExpectedTimes.get(0, 2)) + omegaAtExpectedTimes.get(0, 3) - omegaAtBoundaries.get(0, 3)
                    : omegaAtExpectedTimes.get(0, 0) + probCoalesceBetweenExpectedTimesAndUpperLimit[i] * omegaAtExpectedTimes.get(0, 1) + omegaAtExpectedTimes.get(0, 2) - omegaAtBoundaries.get(0, 2);
            D.put(i, diagonal);
            double thisB = (type == TransitionType.CSC)
                    ? omegasAtBoundaries.get(i).get(0, 3) - omegasAtBoundaries.get(i - 1).get(0, 3)
                    : omegasAtBoundaries.get(i).get(0, 2) - omegasAtBoundaries.get(i - 1).get(0, 2);
            B.put(i - 1, thisB);
        }
        // do U and RR up to states - 2
        for (int i = 0; i < states - 2; i++) {
            double omegaSi = (type == TransitionType.CSC)
                    ? omegasAtExpectedTimes.get(i).get(0, 1) + omegasAtExpectedTimes.get(i).get(0, 2)
                    : omegasAtExpectedTimes.get(i).get(0, 1);
            double omegaSiplus1 = (type == TransitionType.CSC)
                    ? omegasAtExpectedTimes.get(i + 1).get(0, 1) + omegasAtExpectedTimes.get(i + 1).get(0, 2)
                    : omegasAtExpectedTimes.get(i + 1).get(0, 1);
            double thisU = omegaSi * (1 - probCoalesceBetweenExpectedTimesAndUpperLimit[i]) * (1 - probNotCoalesceBetweenTimeIntervals[i + 1]);
            U.put(i, thisU);
            // rho == 0 --> transition is identity, ratios are 0/0. Avoid NaN by setting to 1.
            double thisRR = (rho == 0) ? 1. : omegaSi * probNotCoalesceBetweenExpectedTimes[i] / omegaSiplus1;
            RR.put(i, thisRR);
        }
        // do last for U
        double omegaSi = (type == TransitionType.CSC)
                ? omegasAtExpectedTimes.get(states - 2).get(0, 1) + omegasAtExpectedTimes.get(states - 2).get(0, 2)
                : omegasAtExpectedTimes.get(states - 2).get(0, 1);
        double thisU = omegaSi * (1 - probCoalesceBetweenExpectedTimesAndUpperLimit[states - 2]) * (1 - probNotCoalesceBetweenTimeIntervals[states - 1]);
        U.put(states - 2, thisU);

        ArrayList<DoubleMatrix> res = new ArrayList<DoubleMatrix>();
        res.add(D);
        res.add(B);
        res.add(U);
        res.add(RR);
        return res;
    }

    // note: can compute these in linear time instead of building transition matrix (quadratic)
    public static ArrayList<DoubleMatrix> getLinearTimeDecodingQuantitiesGivenTransition(DoubleMatrix T) {
        int N = T.columns;
        DoubleMatrix D = T.diag();
        DoubleMatrix B = new DoubleMatrix(N - 1);
        for (int i = 0; i < N - 1; i++) {
            // B is below
            B.put(i, T.get(i + 1, i));
        }
        DoubleMatrix U = new DoubleMatrix(N - 1);
        for (int i = 0; i < N - 1; i++) {
            // U is above
            U.put(i, T.get(i, i + 1));
        }
        DoubleMatrix RR = new DoubleMatrix(N - 1);
        for (int i = 0; i < N - 2; i++) {
            // ratio of columns
            if (T.get(i, N - 1) == T.get(i + 1, N - 1)) {
                // avoids 0/0 for totally linked sites in which T = identity
                RR.put(i, 1.);
            } else {
                RR.put(i, T.get(i, N - 1) / T.get(i + 1, N - 1));
            }
        }
        ArrayList<DoubleMatrix> res = new ArrayList<DoubleMatrix>();
        res.add(D);
        res.add(B);
        res.add(U);
        res.add(RR);
        return res;
    }

    private DoubleMatrix transtionMatrix(double r) {
        double[] expIntervals = expectedIntervalTimesPiecewise();
        DoubleMatrix transitionMatrix = new DoubleMatrix(expIntervals.length, expIntervals.length);
        for (int i = 0; i < discretization.length - 1; i++) {
            double timeS = expIntervals[i];
            for (int j = 0; j < discretization.length - 1; j++) {
                double fromTime = discretization[j];
                double toTime = discretization[j + 1];
                transitionMatrix.put(i, j, getTransitionFromStoInterval(r, timeS, fromTime, toTime, type));
            }
        }
        return transitionMatrix;
    }

    // can change this for other models (e.g. piecewise exponential)
    private DoubleMatrix getExponentiatedTransitionMatrix(double N, double r, double time, TransitionType type) {
        double rho = 2 * r * time;
        double eta = 1 / N * time;
        switch (type) {
            case SMC:
                return expm(new DoubleMatrix(new double[][]{{-rho, rho, 0}, {0, -eta, eta}, {0, 0, 0}}));
            case SMC1:
                return expm(new DoubleMatrix(new double[][]{{-rho, rho, 0}, {eta, -2 * eta, eta}, {0, 0, 0}}));
            case CSC:
                return expm(new DoubleMatrix(new double[][]{{-rho, rho, 0, 0}, {eta, -(2 * eta + rho / 2), rho / 2, eta}, {0, 4 * eta, -5 * eta, eta}, {0, 0, 0, 0}}));
            default:
                Utils.exit("Unknown transition matrix requested.");
                return null;
        }
    }

    private double getTransitionFromStoInterval(double r, double timeS, double fromTime, double toTime, TransitionType type) {
        double fromCum = getCumulativeTransitionPobability(r, timeS, fromTime, type);
        double toCum = getCumulativeTransitionPobability(r, timeS, toTime, type);
        double prob = toCum - fromCum;
        return prob;
    }

    public double[] expectedIntervalTimesPiecewise() {
        double[] expectedTimes = new double[discretization.length - 1];
        for (int i = 0; i < discretization.length - 1; i++) {
            double timeFrom = discretization[i];
            double timeTo = discretization[i + 1];
            double expectedThisBin = expectedTimeFromStoT(timeFrom, timeTo);
            expectedTimes[i] = expectedThisBin;
        }
        return expectedTimes;
    }

    private double expectedTimeFromStoT(double timeS, double timeT) {
        int indexFrom = findIntervalForTime(timeS);
        int indexTo = findIntervalForTime(timeT);
        double expected = 0.;
        double rate = 0;
        for (int i = indexFrom; i < indexTo + 1; i++) {
            double time0 = Math.max(timeS, timeVectorPlusInfinity[i]);
            double time1 = Math.min(timeT, timeVectorPlusInfinity[i + 1]);
            double N = sizeVector[i];
            double T = time1 - time0;
            if (time0 == time1) {
                continue;
            }
            double expectedThisPiece = (time1 == Double.POSITIVE_INFINITY)
                    ? Math.exp((timeS - time0) / N) * (N - timeS + time0)
                    : Math.exp(timeS / N) * ((N - timeS + time0) / Math.exp(time0 / N) - (N - timeS + time1) / Math.exp(time1 / N));
            rate = rate - T / N;
            expected += expectedThisPiece;
        }
        // prob having coalesced = 1 - prob not having coalesced
        double norm = 1 - Math.exp(rate);
        expected = expected / norm;
        expected = expected + timeS;
        return expected;
    }

    private double getCumulativeTransitionPobability(double r, double timeS, double timeT, TransitionType type) {
        double cumulative;
        DoubleMatrix Omega;
        if (timeT < timeS) {
            Omega = computeTransitionPiecewiseUpToTimeT(r, timeT, type);
            if (type == TransitionType.CSC) {
                cumulative = Omega.get(0, 3);
            } else {
                cumulative = Omega.get(0, 2);
            }
        } else if (timeT == timeS) {
            Omega = computeTransitionPiecewiseUpToTimeT(r, timeS, type);
            if (type == TransitionType.CSC) {
                cumulative = Omega.get(0, 0) + Omega.get(0, 3);
            } else {
                cumulative = Omega.get(0, 0) + Omega.get(0, 2);
            }
        } else {
            Omega = computeTransitionPiecewiseUpToTimeT(r, timeS, type);
            double cumCoalFromStoT = cumulativeCoalesceFromStoT(timeS, timeT);
            if (type == TransitionType.CSC) {
                cumulative = Omega.get(0, 0) + cumCoalFromStoT * (Omega.get(0, 1) + Omega.get(0, 2)) + Omega.get(0, 3);
            } else {
                cumulative = Omega.get(0, 0) + cumCoalFromStoT * Omega.get(0, 1) + Omega.get(0, 2);
            }
        }
        return cumulative;
    }

    double cumulativeCoalesceFromStoT(double timeS, double timeT) {
        double Nt = getSizeInPiecewiseAtTimeT(timeT);
        double coal = coalesceFromStoT(timeS, timeT);
        double cumulative = 1 - Nt * coal;
        return cumulative;
    }

    private double getSizeInPiecewiseAtTimeT(double timeT) {
        double N = sizeVector[findIntervalForTime(timeT)];
        return N;
    }

    private double coalesceFromStoT(double timeS, double timeT) {
        if (timeT == Double.POSITIVE_INFINITY) {
            return 0.;
        }
        int indexFrom = findIntervalForTime(timeS);
        int indexTo = findIntervalForTime(timeT);
        double rate = 0;
        for (int i = indexFrom; i <= indexTo; i++) {
            rate += (Math.max(timeS, timeVector[i]) - Math.min(timeT, timeVector[i + 1])) / sizeVector[i];
        }
        double Nt = getSizeInPiecewiseAtTimeT(timeT);
        double p = 1 / Nt * Math.exp(rate);
        return p;
    }

    private DoubleMatrix computeTransitionPiecewiseUpToTimeT(double r, double time, TransitionType type) {
        int indexTo = findIntervalForTime(time);
        DoubleMatrix matrix = (type == TransitionType.CSC) ? org.jblas.DoubleMatrix.eye(4) : org.jblas.DoubleMatrix.eye(3);
        for (int i = 0; i <= indexTo - 1; i++) {
            DoubleMatrix M = getExponentiatedTransitionMatrix(sizeVector[i], r, timeVector[i + 1] - timeVector[i], type);
            matrix = matrix.mmul(M);
        }
        DoubleMatrix M = getExponentiatedTransitionMatrix(sizeVector[indexTo], r, time - timeVector[indexTo], type);
        matrix = matrix.mmul(M);
        return matrix;
    }

    private int findIntervalForTime(double time) {

        if (time == Double.POSITIVE_INFINITY) {
            return sizeVector.length - 1;
        }

        for (int i = 0; i < sizeVector.length; i++) {
            if (time >= timeVector[i] && time < timeVector[i + 1]) {
                return i;
            }
        }
        return -1;
    }

    private double notCoalesceFromStoT(double timeS, double timeT) {
        if (timeT == Double.POSITIVE_INFINITY) {
            return 0.;
        }
        int indexFrom = findIntervalForTime(timeS);
        int indexTo = findIntervalForTime(timeT);
        double rate = 0;
        for (int i = indexFrom; i <= indexTo; i++) {
            rate += (Math.max(timeS, timeVector[i]) - Math.min(timeT, timeVector[i + 1])) / sizeVector[i];
        }
        double p = Math.exp(rate);
        return p;
    }

    private DoubleMatrix computeTransitionPiecewiseFromTimeSToTimeT(double r, double timeS, double timeT, TransitionType type) {
        int indexFrom = findIntervalForTime(timeS);
        int indexTo = findIntervalForTime(timeT);
        DoubleMatrix matrix = (type == TransitionType.CSC) ? fourState_identity : threeState_identity;
        for (int i = indexFrom; i <= indexTo; i++) {
            DoubleMatrix M = getExponentiatedTransitionMatrix(sizeVector[i], r, (Math.min(timeT, timeVector[i + 1]) - Math.max(timeS, timeVector[i])), type);
            matrix = matrix.mmul(M);
        }
        return matrix;
    }

    double cumulativeCoalesceFromStoTsmart(double timeS, double timeT) {
        double cumulative = 1 - notCoalesceFromStoT(timeS, timeT);
        return cumulative;
    }

    private Pair<ArrayList<DoubleMatrix>, ArrayList<DoubleMatrix>> getOmegas(double r, TransitionType type) {
        ArrayList<DoubleMatrix> omegasAtBoundaries = new ArrayList<DoubleMatrix>(states + 1);
        ArrayList<DoubleMatrix> omegasAtExpectedTimes = new ArrayList<DoubleMatrix>(states);
        DoubleMatrix latestOmega = (type == TransitionType.CSC) ? fourState_identity : threeState_identity;
        omegasAtBoundaries.add(latestOmega.getRow(0));
        for (int i = 0; i < discretization.length - 1; i++) {
            double intervalStartTime = discretization[i];
            double intervalExpTime = expectedTimes[i];
            double intervalEndTime = discretization[i + 1];
            DoubleMatrix M = computeTransitionPiecewiseFromTimeSToTimeT(r, intervalStartTime, intervalExpTime, type);
            latestOmega = latestOmega.mmul(M);
            omegasAtExpectedTimes.add(latestOmega.getRow(0));
            if (intervalEndTime == Double.POSITIVE_INFINITY) {
                M = (type == TransitionType.CSC) ? fourState_timeInf : threeState_timeInf;
            } else {
                M = computeTransitionPiecewiseFromTimeSToTimeT(r, intervalExpTime, intervalEndTime, type);
            }
            latestOmega = latestOmega.mmul(M);
            omegasAtBoundaries.add(latestOmega.getRow(0));
        }
        return new Pair<ArrayList<DoubleMatrix>, ArrayList<DoubleMatrix>>(omegasAtBoundaries, omegasAtExpectedTimes);
    }

    private void computeCoalescentVectors() {
        probNotCoalesceBetweenExpectedTimes = new double[expectedTimes.length - 1];
        probNotCoalesceBetweenTimeIntervals = new double[expectedTimes.length];
        probCoalesceBetweenExpectedTimesAndUpperLimit = new double[expectedTimes.length];
        for (int i = 0; i < expectedTimes.length; i++) {
            double timeFrom = discretization[i];
            double timeTo = discretization[i + 1];
            double expTimeFrom = expectedTimes[i];
            double p;
            if (i < expectedTimes.length - 1) {
                double expTimeTo = expectedTimes[i + 1];
                p = notCoalesceFromStoT(expTimeFrom, expTimeTo);
                probNotCoalesceBetweenExpectedTimes[i] = p;
            }
            p = notCoalesceFromStoT(timeFrom, timeTo);
            probNotCoalesceBetweenTimeIntervals[i] = p;
            p = cumulativeCoalesceFromStoTsmart(expTimeFrom, timeTo);
            probCoalesceBetweenExpectedTimesAndUpperLimit[i] = p;
        }
        // do U and RR up to states - 2
        columnRatios = new DoubleMatrix(states - 1);
        for (int i = 1; i < states - 1; i++) {
            double thisCR = probNotCoalesceBetweenTimeIntervals[i] * (1 - probNotCoalesceBetweenTimeIntervals[i + 1]) / (1 - probNotCoalesceBetweenTimeIntervals[i]);
            if (thisCR == Double.NaN) {
                thisCR = 1.;
            }
            columnRatios.put(i, thisCR);
        }
    }

    public double[] getCoalDist() {
        double[] coalDist = new double[discretization.length - 1];
        double lastCoal = 0.;
        for (int i = 1; i < discretization.length; i++) {
            double from = discretization[i - 1];
            double to = discretization[i];
            double coal = cumulativeCoalesceFromStoT(0., to);
            coalDist[i - 1] = coal - lastCoal;
            lastCoal = coal;
        }
        return coalDist;
    }

//    static void unitTest() {
//        double[] sizeVector = {2000000., 200000., 20000., 20000., 20000., 20000., 20000., 20000.};
//        double[] timeVector = {0., 100., 200., 300., 400., 500., 600., Double.POSITIVE_INFINITY};
//        double[] discretization = {0., 50., 110., 305., 500., 633., Double.POSITIVE_INFINITY};
//        double recDist = 1.2e-8;
//        double mu = 1.65e-8;
//        ArrayList<Double> sizeVectorArray = new ArrayList<Double>();
//        for (double d : sizeVector) {
//            sizeVectorArray.add(d);
//        }
//        ArrayList<Double> timeVectorArray = new ArrayList<Double>();
//        for (double d : timeVector) {
//            timeVectorArray.add(d);
//        }
////        Utils.print(Transition.getTimeExponentialQuantiles(10, timeVectorArray, sizeVectorArray));
////        Utils.print(Transition.getTimeErlangQuantiles(10, timeVectorArray, sizeVectorArray));
//        Transition transition = new Transition(timeVector, sizeVector, discretization, TransitionType.CSC);
//        transition.computeCoalescentVectors();
////        Utils.print(transition.findIntervalForTime(0.));
////        Utils.print(transition.findIntervalForTime(101.1));
////        Utils.print(transition.findIntervalForTime(10000000));
////        Utils.print(transition.expectedTimeFromStoT(10, 201));
////        Utils.print(transition.cumulativeCoalesceFromStoTsmart(10, 201));
////        Utils.print(transition.notCoalesceFromStoT(10, 201));
////        Utils.print(transition.getExponentiatedTransitionMatrix(10000, 2.1e-5, 1201.1, TransitionType.SMC));
////        Utils.print(transition.getExponentiatedTransitionMatrix(10000, 2.1e-5, 1201.1, TransitionType.SMC1));
////        Utils.print(transition.getExponentiatedTransitionMatrix(10000, 2.1e-5, 1201.1, TransitionType.CSC));
////        Utils.print(Utils.doubleMatrixToString(transition.computeTransitionPiecewiseFromTimeSToTimeT(2.3e-5, 10, 201, Transition.TransitionType.CSC)));
////        Utils.print(transition.computeTransitionPiecewiseFromTimeSToTimeT(2.3e-5, 10, 201, Transition.TransitionType.CSC));
////        Utils.printMatrixOfDoubles(transition.probNotCoalesceBetweenExpectedTimes);
////        Utils.printMatrixOfDoubles(transition.probNotCoalesceBetweenTimeIntervals);
////        Utils.printMatrixOfDoubles(transition.probCoalesceBetweenExpectedTimesAndUpperLimit);
////        Utils.print(transition.columnRatios);
////        Utils.printMatrixOfDoubles(transition.timeVector);
////        Utils.printMatrixOfDoubles(transition.sizeFromVector);
////        Utils.print(transition.type);
////        Utils.printMatrixOfDoubles(transition.timeVectorPlusInfinity);
////        Utils.print(transition.states);
////        Utils.printMatrixOfDoubles(transition.expectedTimes);
////        Utils.printMatrixOfDoubles(transition.discretization);
////        Utils.print(transition.getCumulativeTransitionPobability(2.6e-6, 102.1, 1002.4, TransitionType.CSC));
////        Utils.print(transition.getTransitionFromStoInterval(2.6e-6, 1006.5, 102.1, 1002.4, Transition.TransitionType.CSC));
////        ArrayList<DoubleMatrix> res = transition.getLinearTimeDecodingQuantitiesGivenTransition(transition.getTransitionForDistance(2.4e-2));
////        Utils.print(res.get(0));
////        Utils.print(res.get(1));
////        Utils.print(res.get(2));
////        Utils.print(res.get(3));
////        Utils.print(transition.getSizeInPiecewiseAtTimeT(1002.4));
////        Utils.print(transition.coalesceFromStoT(102.1, 1002.4));
////        Utils.print(transition.computeTransitionPiecewiseUpToTimeT(2.1e-4, 1021.0, TransitionType.CSC));
////        Utils.print(transition.transtionMatrix(2.1e-2));
////        ArrayList<DoubleMatrix> res = transition.getLinearTimeDecodingQuantitiesAndMatrixGivenDistance(2.1e-2);
////        Utils.print(res.get(0));
////        Utils.print(res.get(1));
////        Utils.print(res.get(2));
////        Utils.print(res.get(3));
////        Utils.printMatrixOfDoubles(transition.getCoalDist());
////        DoubleMatrix T = transition.transtionMatrix(recDist);
////        double[] coalDist = transition.getCoalDist();
////        for (int i = 0; i < coalDist.length; i++) {
////            Utils.print(transition.discretization[i] + " " + transition.discretization[i + 1] + " " + coalDist[i]);
////        }
//    }

}
