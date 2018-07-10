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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Date;
import java.util.jar.JarFile;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import org.jblas.DoubleMatrix;

public class Utils {

    public static String doubleMatrixToString(DoubleMatrix m) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < m.columns; j++) {
                sb.append(m.get(i, j));
                sb.append("\t");
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    public static String doubleToRoundedString(double d, int decimals) {
        String format = "%." + decimals + "f";
        return String.format(format, d);
    }

    public static void printMatrixOfDoubles(double[][] m) {
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                System.out.format("\t%.10f", m[i][j]);
            }
            System.out.print("\n");
        }
        System.out.print("\n");
    }

    public static void printMatrixOfDoubles(double[] m) {
        for (int i = 0; i < m.length; i++) {
            System.out.format("\t%.10f", m[i]);
        }
        System.out.print("\n");
    }

    public static String doubleMatrixToString(double[] m) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < m.length; i++) {
            sb.append(m[i]);
            sb.append("\t");
        }
        sb.append("\n");
        return sb.toString();
    }

    public static String doubleMatrixToString(double[][] m) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                sb.append(m[i][j]);
                sb.append("\t");
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    public static boolean fileExists(String fileName) {
        if (fileName == null) {
            return false;
        }
        File f = new File(fileName);
        if (f.exists() && !f.isDirectory()) {
            return true;
        }
        return false;
    }

    public static void exit(String error) {
        System.err.println("ERROR:\t" + error);
        System.exit(1);
    }

    public static void print(Object message) {
        System.out.println(message.toString());
    }

    public static void print() {
        System.out.println();
    }

    public static void printf(Object message) {
        System.out.print(message.toString());
    }

    public static void warning(String message) {
        System.err.println("Warning:\t" + message);
    }

    public static BufferedReader openFile(String fileName) {
        try {
            Utils.print("Opening file " + fileName + " for reading.");
            if (isGzipped(fileName)) {
                return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fileName))));
            } else {
                return new BufferedReader(new FileReader(fileName));
            }
        } catch (IOException e) {
            Utils.exit("Could not open file " + fileName + ". " + e.toString());
            return null;
        }
    }

    public static BufferedWriter openGzipFileForWriting(String fileName) {
        try {
            Utils.print("Opening file " + fileName + " for writing.");
            BufferedWriter writer = null;
            GZIPOutputStream zip = new GZIPOutputStream(
                    new FileOutputStream(new File(fileName)));
            writer = new BufferedWriter(
                    new OutputStreamWriter(zip, "UTF-8"));
            return writer;
        } catch (IOException e) {
            Utils.exit("Could not open file " + fileName + ". " + e.toString());
            return null;
        }
    }

    public static BufferedWriter openFileForWriting(String fileName) {
        try {
            Utils.print("Opening file " + fileName + " for writing");
            return new BufferedWriter(new OutputStreamWriter(new FileOutputStream(fileName), "utf-8"));
        } catch (IOException e) {
            Utils.exit("Could not open file " + fileName + ". " + e.toString());
            return null;
        }
    }

    public static boolean isGzipped(String fileName) throws FileNotFoundException, IOException {
        File f = new File(fileName);
        InputStream is = new FileInputStream(f);
        byte[] signature = new byte[2];
        int nread = is.read(signature); //read the gzip signature
        return nread == 2 && signature[0] == (byte) 0x1f && signature[1] == (byte) 0x8b;
    }

    public static double getSumOfVector(double[] vector) {
        double sum = 0.f;
        for (int i = 0; i < vector.length; i++) {
            sum += vector[i];
        }
        return sum;
    }

    public static double getSumOfMatrixColumn(double[][] matrix, int colIndex) {
        double sum = 0.f;
        for (int i = 0; i < matrix.length; i++) {
            sum += matrix[i][colIndex];
        }
        return sum;
    }

    public static double[] elementWiseMultVectorScalar(double[] vector, double val) {
        double[] ret = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            ret[i] = vector[i] * val;
        }
        return ret;
    }

    public static double[] elementWiseMultVectorVector(double[] vector, double[] factors) {
        double[] ret = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            ret[i] = vector[i] * factors[i];
        }
        return ret;
    }

    public static double[][] elementWiseMultMatrixMatrix(double[][] matrix1, double[][] matrix2) {
        double[][] ret = new double[matrix1.length][matrix1[0].length];
        for (int i = 0; i < matrix1.length; i++) {
            for (int j = 0; j < matrix1[0].length; j++) {
                ret[i][j] = matrix1[i][j] * matrix2[i][j];
            }
        }
        return ret;
    }

    public static void sumToMatrix(double[][] matrix1, double[][] matrix2) {
        for (int i = 0; i < matrix1.length; i++) {
            for (int j = 0; j < matrix1[0].length; j++) {
                matrix1[i][j] += matrix2[i][j];
            }
        }
    }

    public static double[][] matrixMultiplySquareMatrices(double[][] matrix1, double[][] matrix2) {
        if (matrix1.length != matrix2.length || matrix1[0].length != matrix2[0].length) {
            Utils.exit("Matrix dimensions differ.");
        }
        int len = matrix1.length;
        double[][] res = new double[matrix1.length][matrix1[0].length];
        for (int row = 0; row < len; row++) {
            for (int column = 0; column < len; column++) {
                for (int k = 0; k < len; k++) {
                    res[row][column] += matrix1[row][k] * matrix2[k][column];
                }
            }
        }
        return res;
    }

    public static double[][] normalizeMatrixColumns(double[][] matrix) {
        double[][] ret = new double[matrix.length][matrix[0].length];
        for (int j = 0; j < matrix[0].length; j++) {
            double sum = 0.f;
            for (int i = 0; i < matrix.length; i++) {
                sum += matrix[i][j];
            }
            for (int i = 0; i < matrix.length; i++) {
                ret[i][j] = matrix[i][j] / sum;
            }
        }
        return ret;
    }

    public static void fillMatrixColumn(double[][] matrix, double[] vector, int pos) {
        for (int i = 0; i < vector.length; i++) {
            matrix[i][pos] = vector[i];
        }
    }

    private String getJarFolder() {
        String name = this.getClass().getName().replace('.', '/');
        String s = this.getClass().getResource("/" + name + ".class").toString();
        s = s.replace('/', File.separatorChar);
        s = s.substring(0, s.indexOf(".jar") + 4);
        s = s.substring(s.lastIndexOf(':') - 1);
        return s.substring(0, s.lastIndexOf(File.separatorChar) + 1);
    }

}

class Pair<K, V> {

    private K key;
    private V value;

    public Pair(K key, V value) {
        this.key = key;
        this.value = value;
    }

    public K getKey() {
        return this.key;
    }

    public V getValue() {
        return value;
    }
}
