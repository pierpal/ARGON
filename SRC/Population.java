package ARGON;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

/**
 *
 * @author Pier Palamara <ppalama@hsph.harvard.edu>
 */
public class Population {

    ArrayList<ArrayList<Long>> populationSizes = new ArrayList<ArrayList<Long>>();
    ArrayList<double[][]> samplingProbabilities = new ArrayList<double[][]>();
    ArrayList<ArrayList<Long>> populationIDOffsets = new ArrayList<ArrayList<Long>>();

    Population(String populationFile, boolean quiet) throws Exception {
        FileInputStream fstream = new FileInputStream(populationFile);
        DataInputStream in = new DataInputStream(fstream);
        BufferedReader br = new BufferedReader(new InputStreamReader(in));
        String strLine;
        int count = 0;
        long gen = 0;
        while ((strLine = br.readLine()) != null) {
            String[] tokens = strLine.split("\\s+");
            count++;
            ArrayList<double[]> migrationList;
            // if odd line, will contain pop sizes
            if (count % 2 == 1) {
                ArrayList<Long> thisGenSizes = new ArrayList<Long>();
                gen = Long.parseLong(tokens[0]);
                if (populationSizes.size() != gen) {
                    System.err.println("Detected missing generation in population file. Observed generation " + gen + ", expected " + populationSizes.size() + ".");
                    System.exit(0);
                }
                long totSize = 0L;
                ArrayList<Long> currOffsets = new ArrayList<Long>();
                currOffsets.add(totSize);
                for (int tokenIndex = 1; tokenIndex < tokens.length; tokenIndex++) {
                    long thisSize = Long.parseLong(tokens[tokenIndex]);
                    thisGenSizes.add(thisSize);
                    totSize += thisSize;
                    if (tokenIndex < tokens.length - 1) {
                        currOffsets.add(totSize);
                    }
                }
                populationIDOffsets.add(currOffsets);
                populationSizes.add(thisGenSizes);
                if (!quiet) {
                    System.out.println("Population sizes at generation " + gen + " are " + populationSizes.get(populationSizes.size() - 1));
                }
            } else {
                // if even line, will contain topology and migration
                migrationList = new ArrayList<double[]>();
                int maxFrom = 0;
                int maxTo = 0;
                for (int tokenIndex = 0; tokenIndex < tokens.length; tokenIndex++) {
                    String[] migrationParameters = tokens[tokenIndex].split("-|_", 3);
                    int fromPop = Integer.parseInt(migrationParameters[0]);
                    if (fromPop > maxFrom) {
                        maxFrom = fromPop;
                    }
                    int toPop = Integer.parseInt(migrationParameters[1]);
                    if (toPop > maxTo) {
                        maxTo = toPop;
                    }
                    if (migrationParameters.length == 3) {
                        double value = Double.parseDouble(migrationParameters[2]);
                        migrationList.add(new double[]{(double) fromPop, (double) toPop, value});
                    } else if (migrationParameters.length == 2) {
                        double value = 1.;
                        migrationList.add(new double[]{(double) fromPop, (double) toPop, value});
                    }
                }
                double[][] migrationMatrix = new double[maxFrom][maxTo];
                for (int i = 0; i < migrationList.size(); i++) {
                    int popFrom = (int) migrationList.get(i)[0];
                    int popTo = (int) migrationList.get(i)[1];
                    double mig = migrationList.get(i)[2];
                    migrationMatrix[popFrom - 1][popTo - 1] = mig;
                }
                boolean printWarning = false;
                for (int i = 0; i < maxFrom; i++) {
                    double sum = 0.;
                    for (int j = 0; j < maxTo; j++) {
                        sum += migrationMatrix[i][j];
                    }
                    if (sum != 1.) {
                        printWarning = true;
                        for (int j = 0; j < maxTo; j++) {
                            migrationMatrix[i][j] /= sum;
                        }
                    }
                }
                if (printWarning) {
                    System.err.println("Warning: some migration matrix rows did not sum to 1.0. Normalized.");
                }
                samplingProbabilities.add(migrationMatrix);
                if (!quiet) {
                    printMatrix(migrationMatrix);
                }
            }
        }
        if (count % 2 != 1) {
            System.err.println("Format error. Demographic history contains even number of lines. Must alternate sizes and split/migration lines.");
            System.exit(count);
        }
        if (populationSizes.get(populationSizes.size() - 1).size() == 1) {
            samplingProbabilities.add(new double[][]{{1.}});
        }
        if (!quiet) {
            System.out.print("Population sizes from generation " + (populationSizes.size()) + " on are " + populationSizes.get(populationSizes.size() - 1) + ". ");
            printMatrix(getSamplingProbabilities(populationSizes.size()));
            System.out.println();
        }
        if (populationSizes.get(populationSizes.size() - 1).size() != samplingProbabilities.get(samplingProbabilities.size() - 1).length) {
            System.err.println("The demographic model contains " + populationSizes.get(populationSizes.size() - 1).size()
                    + " populations at generation " + (populationSizes.size() - 1) + ", but only "
                    + samplingProbabilities.get(samplingProbabilities.size() - 1).length + " rows are present in the final migration matrix.");
            System.exit(0);
        }
        in.close();
    }

    Population(Long s) throws Exception {
        ArrayList<Long> sizes = new ArrayList<Long>();
        sizes.add(s);
        populationSizes.add(sizes);
        double[][] migrations = new double[1][1];
        migrations[0][0] = 1.;
        samplingProbabilities.add(migrations);
        ArrayList<Long> currOffsets = new ArrayList<Long>();
        currOffsets.add(0L);
        populationIDOffsets.add(currOffsets);
    }

    static void printMatrix(double[][] migrationMatrix) {
        StringBuilder migrationMatrixString = new StringBuilder("");
        for (int i = 0; i < migrationMatrix.length; i++) {
            if (i != 0) {
                migrationMatrixString.append("\n");
            }
            for (int j = 0; j < migrationMatrix[0].length; j++) {
                migrationMatrixString.append("\t");
                double rounded = Math.round(migrationMatrix[i][j] * 100000000000000.) / 100000000000000.;
                migrationMatrixString.append(Double.toString(rounded));
            }
        }
        System.out.println("Migration matrix:\n" + migrationMatrixString);
    }

    long getIDOffsetAt(int gen, int pop) {
        if (populationIDOffsets.size() <= gen) {
            return populationIDOffsets.get(populationSizes.size() - 1).get(pop);
        } else {
            return populationIDOffsets.get(gen).get(pop);
        }
    }

    long getSizeAt(int gen, int pop) {
        if (populationSizes.size() <= gen) {
            return populationSizes.get(populationSizes.size() - 1).get(pop);
        } else {
            return populationSizes.get(gen).get(pop);
        }
    }

    int getNumPopsAt(int gen) {
        return (populationSizes.size() <= gen) ? populationSizes.get(populationSizes.size() - 1).size() : populationSizes.get(gen).size();
    }

    double[][] getSamplingProbabilities(int gen) {
        return (samplingProbabilities.size() <= gen) ? samplingProbabilities.get(samplingProbabilities.size() - 1) : samplingProbabilities.get(gen);
    }

}
