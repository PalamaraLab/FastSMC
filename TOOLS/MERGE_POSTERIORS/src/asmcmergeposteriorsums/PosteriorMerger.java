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


package asmcmergeposteriorsums;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;

public class PosteriorMerger {

    public PosteriorMerger(String fileList, String outFileRoot, boolean normalize) throws IOException {
        this.fileListName = fileList;
        this.outFileRoot = outFileRoot;
        this.normalize = normalize;
        this.mergeMatricesFileList();
    }

    ArrayList<HashSet<Integer>> skipMapLinesNotInIntersection = new ArrayList<HashSet<Integer>>();
    HashMap<String, Integer> mapCount = new HashMap<String, Integer>();
    HashMap<String, String> SNPchr = new HashMap<String, String>();
    HashMap<String, String> SNPname = new HashMap<String, String>();
    HashMap<String, String> SNPgen = new HashMap<String, String>();
    HashMap<String, String> SNPphys = new HashMap<String, String>();

    public int getIntersectionMap(ArrayList<String> fileList) throws IOException {
        String line;
        for (String s : fileList) {
            BufferedReader br = Utils.openFile(s + ".map.gz");
            while ((line = br.readLine()) != null) {
                String[] sL = line.split("\\s+");
                String SNP = sL[0] + "\t" + sL[1] + "\t" + sL[3];
                SNPchr.put(SNP, sL[0]);
                SNPname.put(SNP, sL[1]);
                SNPgen.put(SNP, sL[2]);
                SNPphys.put(SNP, sL[3]);
                Integer SNPcnt = mapCount.get(SNP);
                if (SNPcnt == null) {
                    SNPcnt = 0;
                }
                SNPcnt++;
                if (!mapCount.containsKey(SNP)) {
                    mapArray.add(SNP);
                }
                mapCount.put(SNP, SNPcnt);
            }
        }
        int allBatches = 0;
        for (String s : mapCount.keySet()) {
            if (mapCount.get(s) == fileList.size()) {
                allBatches++;
                map.add(s);
            } else {
                mapArray.remove(s);
            }
        }
        return allBatches;
    }

    private void getRowColumns(String file) throws IOException {
        BufferedReader br = Utils.openFile(file);
        String line;
        int cnt = 0;
        int col = 0;
        while ((line = br.readLine()) != null) {
            col = line.split("\\s+").length;
            if (cnt > 0 && columns != col) {
                Utils.exit("Columns are not the same as in previous row: " + col);
            }
            columns = col;
            cnt++;
        }
        rows = cnt;
    }

    private void fillMatrixUsingRootList(String root, String string, float[][] matrix) throws IOException {
        String mapFile = root + ".map.gz";
        BufferedReader br = Utils.openFile(mapFile);
        ArrayList<String> thisMap = new ArrayList<String>();
        String line;
        while ((line = br.readLine()) != null) {
            String[] sL = line.split("\\s+");
            String SNP = sL[0] + "\t" + sL[1] + "\t" + sL[3];
            thisMap.add(SNP);
        }
        String name = root + "." + string + ".sumOverPairs.gz";
        br = Utils.openFile(name);
        int index = 0, cnt = 0;
        while ((line = br.readLine()) != null) {
            if (map.contains(thisMap.get(index))) {
                String[] splitStr = line.split("\\s+");
                for (int c = 0; c < splitStr.length; c++) {
                    matrix[cnt][c] += Float.parseFloat(splitStr[c]);
                }
                cnt++;
            }
            index++;
        }
        if (cnt != matrix.length) {
            Utils.exit("Rows mismatch: " + matrix.length + " --> " + cnt);
        }
    }

    public void mergeMatricesFileList() throws IOException {
        BufferedReader br = Utils.openFile(fileListName);
        String line;
        ArrayList<String> fileList = new ArrayList<String>();
        while ((line = br.readLine()) != null) {
            fileList.add(line);
        }
        rows = getIntersectionMap(fileList);
        columns = Utils.openFile(fileList.get(0) + ".00.sumOverPairs.gz").readLine().split("\\s+").length;
        Utils.print("Rows: " + rows + "\tcolumns: " + columns);
        sum00 = new float[rows][columns];
        sum01 = new float[rows][columns];
        sum11 = new float[rows][columns];
        sumAll = new float[rows][columns];
        for (int i = 0; i < fileList.size(); i++) {
            String root = fileList.get(i);
            this.fillMatrixUsingRootList(root, "00", sum00);
            this.fillMatrixUsingRootList(root, "01", sum01);
            this.fillMatrixUsingRootList(root, "11", sum11);
            for (int r = 0; r < rows; r++) {
                for (int c = 0; c < columns; c++) {
                    sumAll[r][c] = sum00[r][c] + sum01[r][c] + sum11[r][c];
                }
            }
        }
        if (normalize) {
            sum00 = Utils.normalizeMatrixRows(sum00);
            sum01 = Utils.normalizeMatrixRows(sum01);
            sum11 = Utils.normalizeMatrixRows(sum11);
            sumAll = Utils.normalizeMatrixRows(sumAll);
        }
    }

    String fileRoot;
    String fileListName;
    int jobs = -1;
    String outFileRoot;
    boolean normalize;

    TreeSet<String> map = new TreeSet<String>();
    ArrayList<String> mapArray = new ArrayList<String>();

    int rows, columns;
    float[][] sum00;
    float[][] sum01;
    float[][] sum11;
    float[][] sumAll;

    float[] expTime00;
    float[] expTime01;
    float[] expTime11;
    float[] expTimeAll;

    public PosteriorMerger(String fileName, int jobs, String outFileRoot, boolean normalize) throws IOException {
        this.fileRoot = fileName;
        this.jobs = jobs;
        this.outFileRoot = outFileRoot;
        this.normalize = normalize;
        this.mergeMatricesFileRoot();
    }

    public void mergeMatricesFileRoot() throws IOException {

        String name = fileRoot + "." + 1 + "-" + jobs + ".00.sumOverPairs.gz";
        getRowColumns(name);
        sum00 = new float[rows][columns];
        sum01 = new float[rows][columns];
        sum11 = new float[rows][columns];
        sumAll = new float[rows][columns];
        fillMatrixUsingJobIndeces("00", sum00);
        fillMatrixUsingJobIndeces("01", sum01);
        fillMatrixUsingJobIndeces("11", sum11);
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < columns; c++) {
                sumAll[r][c] = sum00[r][c] + sum01[r][c] + sum11[r][c];
            }
        }
        if (normalize) {
            sum00 = Utils.normalizeMatrixRows(sum00);
            sum01 = Utils.normalizeMatrixRows(sum01);
            sum11 = Utils.normalizeMatrixRows(sum11);
            sumAll = Utils.normalizeMatrixRows(sumAll);
        }
    }

    public void computeCoalescentTimes(String infoFile) throws IOException {
        BufferedReader br = Utils.openFile(infoFile);
        String line;
        int cnt = 0;
        float[] coalMeans = new float[columns];
        while ((line = br.readLine()) != null) {
            String[] splitStr = line.split("\\s+");
            coalMeans[cnt] = Float.parseFloat(splitStr[1]);
            cnt++;
        }
        if (cnt != columns) {
            Utils.exit("Wrong number of rows in info: " + cnt);
        }
        expTime00 = getCoalPerLine(sum00, coalMeans);
        expTime01 = getCoalPerLine(sum01, coalMeans);
        expTime11 = getCoalPerLine(sum11, coalMeans);
        expTimeAll = getCoalPerLine(sumAll, coalMeans);
    }

    private float[] getCoalPerLine(float[][] matrix, float[] coalMeans) {
        float[] ret = new float[rows];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                ret[i] += matrix[i][j] * coalMeans[j];
            }
        }
        return ret;
    }

    private void fillMatrixUsingJobIndeces(String string, float[][] matrix) throws IOException {
        for (int job = 1; job <= jobs; job++) {
            String name = fileRoot + "." + job + "-" + jobs + "." + string + ".sumOverPairs.gz";
            BufferedReader br = Utils.openFile(name);
            String line;
            int cnt = 0;
            while ((line = br.readLine()) != null) {
                String[] splitStr = line.split("\\s+");
                for (int c = 0; c < splitStr.length; c++) {
                    matrix[cnt][c] += Float.parseFloat(splitStr[c]);
                }
                cnt++;
            }
            if (cnt != rows) {
                Utils.exit("Rows mismatch: " + cnt);
            }
        }
    }

    public void writeMatricesToFile(String outFileRoot) throws IOException {
        BufferedWriter br = Utils.openGzipFileForWriting(outFileRoot + ".merged.00.sumOverPairs.gz");
        for (int i = 0; i < sum00.length; i++) {
            br.write(String.valueOf(sum00[i][0]));
            for (int j = 1; j < sum00[0].length; j++) {
                br.write("\t" + String.valueOf(sum00[i][j]));
            }
            br.write("\n");
        }
        br.close();
        br = Utils.openGzipFileForWriting(outFileRoot + ".merged.01.sumOverPairs.gz");
        for (int i = 0; i < sum01.length; i++) {
            br.write(String.valueOf(sum01[i][0]));
            for (int j = 1; j < sum01[0].length; j++) {
                br.write("\t" + String.valueOf(sum01[i][j]));
            }
            br.write("\n");
        }
        br.close();
        br = Utils.openGzipFileForWriting(outFileRoot + ".merged.11.sumOverPairs.gz");
        for (int i = 0; i < sum11.length; i++) {
            br.write(String.valueOf(sum11[i][0]));
            for (int j = 1; j < sum11[0].length; j++) {
                br.write("\t" + String.valueOf(sum11[i][j]));
            }
            br.write("\n");
        }
        br.close();
        br = Utils.openGzipFileForWriting(outFileRoot + ".merged.sumOverPairs.gz");
        for (int i = 0; i < sumAll.length; i++) {
            br.write(String.valueOf(sumAll[i][0]));
            for (int j = 1; j < sumAll[0].length; j++) {
                br.write("\t" + String.valueOf(sumAll[i][j]));
            }
            br.write("\n");
        }
        br.close();
    }

    public void writeExpCoalToFile(String outFileRoot) throws IOException {
        BufferedWriter br = Utils.openGzipFileForWriting(outFileRoot + ".merged.00.expCoalTime.gz");
        for (int i = 0; i < expTime00.length; i++) {
            br.write(String.valueOf(expTime00[i]) + "\n");
        }
        br.close();
        br = Utils.openGzipFileForWriting(outFileRoot + ".merged.01.expCoalTime.gz");
        for (int i = 0; i < expTime01.length; i++) {
            br.write(String.valueOf(expTime01[i]) + "\n");
        }
        br.close();
        br = Utils.openGzipFileForWriting(outFileRoot + ".merged.11.expCoalTime.gz");
        for (int i = 0; i < expTime11.length; i++) {
            br.write(String.valueOf(expTime11[i]) + "\n");
        }
        br.close();
        br = Utils.openGzipFileForWriting(outFileRoot + ".merged.expCoalTime.gz");
        for (int i = 0; i < expTimeAll.length; i++) {
            br.write(String.valueOf(expTimeAll[i]) + "\n");
        }
        br.close();
        if (!map.isEmpty()) {
            br = Utils.openGzipFileForWriting(outFileRoot + ".merged.map.gz");
            for (String s : mapArray) {
                String line = SNPchr.get(s) + "\t" + SNPname.get(s) + "\t" + SNPgen.get(s) + "\t" + SNPphys.get(s) + "\n";
                br.write(line);
            }
            br.close();
        }
    }

}
