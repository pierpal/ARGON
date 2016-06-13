package ARGON;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

/**
 *
 * @author Pier Palamara <ppalama@hsph.harvard.edu>
 */
public class Argon {

    private static final int MEG = (int) (Math.pow(1024, 2));
    private static final int bufferSize = 100 * MEG;

    final static String version = "0.1.160415";
    final static String date = "April 15, 2016";
    final static String authors = "Pier Palamara";
    final static String contact = "ppalama AT hsph DOT harvard " + "DOT edu";
    static final String versionString = "ARGON version " + version + " (" + date + "; " + authors + ", " + contact + ")";

    boolean quiet = false;

    int thresholdForPrintingBooleanRatherThanList;

    Random generator;

    final boolean shrinkSequenceData;
    final boolean writeOutFiles;
    final boolean printIBD;
    final boolean printMUT;

    final boolean debug = false;

    final boolean printAFS;
    final boolean writeAFS;
    final boolean usePopFile;
    final boolean useMapFile;
    final boolean specifiedMutRate;
    final boolean specifiedRecRate;
    final boolean doGeneConversion;
    final boolean gzipOutput;
    final boolean showDefaultsAndExamples;
    final boolean outputAlleleAge;
    final boolean seqOutFormat;
    final boolean hapsOutFormat;
    final boolean outputTrees;

    final float recRate;
    final float MAF;
    float mutationRate = 0.0000000165f;
    final float recRateFactor;
    final float IBDthresholdCM;
    final float geneConvRelativeRate;
    final float probGeneConversionGivenRec;
    final float GCPhysicalLengthRateForSampling;
    final float GCGeneticLengthRateForSampling;

    final TreeMap<Integer, Float> recRateMap;
    final TreeMap<Integer, Integer> physToGen;
    final TreeMap<Integer, Float> mutRateMap;
    final TreeMap<Integer, ArrayList<Integer>> recToPhysicalMap;

    final ArrayList<Integer> physicalPosMapArray = new ArrayList<Integer>();
    final ArrayList<Float> recRateArray = new ArrayList<Float>();
    final ArrayList<Float> mutRateArray = new ArrayList<Float>();
    final ArrayList<Integer> genPosArray = new ArrayList<Integer>();

    final HashMap<Integer, Integer> mapIndex = new HashMap<Integer, Integer>();

    final String popSizeFile;
    final String mapFile;
    final String fileBase;
    final String chromosome = "1";
    float[] AFS; // TODO have a matrix pop x pop
    int totMut = 0;

    final IBD IBDFactory;

    int numBlocksLeftToTrack;
    final int numPopsAtZero;
    final Long constantPopSize;
    final int sampleSizes[];
    final int genomeSizePhysical;
    int genomeSizeGenetic;
    final int minBlockSize;
    final int totModernSamples;
    final int meanGeneConversionTractLength;
    final int minimumGeneConversionTractLength;
    final int minimumGeneConversionTractLengthGenetic;

    final long seed;

    FileOutputStream MUToutput = null, IBDoutput = null, MAPoutput = null, AGEoutput = null, treeOutput = null;
    Writer MUTwriter = null, IBDwriter = null, MAPwriter = null, AGEwriter = null, treeWriter = null;

    // simulation, visit
    double ARG_runtime;
    Population demographicModel;
    ArrayList<Block> leaves;
    ArrayList<Block> roots;
    double Visit_runtime;
    int numMutations = 0;

    public Argon(
            boolean quietArg,
            boolean shrinkSequenceDataArg,
            boolean writeOutFilesArg,
            boolean printIBDArg,
            boolean printMUTArg,
            boolean printAFSArg,
            boolean writeAFSArg,
            boolean usePopFileArg,
            boolean useMapFileArg,
            boolean specifiedMutRateArg,
            boolean specifiedRecRateArg,
            boolean doGeneConversionArg,
            boolean gzipOutputArg,
            boolean showDefaultsAndExamplesArg,
            boolean outputAlleleAgeArg,
            boolean seqOutFormatArg,
            boolean hapsOutFormatArg,
            boolean outputTreesArg,
            float recRateArg,
            float MAFArg,
            float mutationRateArg,
            float recRateFactorArg,
            float IBDthresholdCMArg,
            float geneConvRelativeRateArg,
            float probGeneConversionGivenRecArg,
            String popSizeFileArg,
            String mapFileArg,
            String fileBaseArg,
            int numPopsArg,
            int sampleSizesArg[],
            int genomeSizePhysicalArg,
            int minBlockSizeArg,
            int totModernSamplesArg,
            int meanGeneConversionTractLengthArg,
            int minimumGeneConversionTractLengthArg,
            Long constantPopSizeArg,
            long seedArg
    ) throws Exception {

        this.quiet = quietArg;
        this.shrinkSequenceData = shrinkSequenceDataArg;
        this.writeOutFiles = writeOutFilesArg;
        this.printIBD = printIBDArg;
        this.printMUT = printMUTArg;
        this.printAFS = printAFSArg;
        this.writeAFS = writeAFSArg;
        this.usePopFile = usePopFileArg;
        this.useMapFile = useMapFileArg;
        this.specifiedMutRate = specifiedMutRateArg;
        this.specifiedRecRate = specifiedRecRateArg;
        this.doGeneConversion = doGeneConversionArg;
        this.gzipOutput = gzipOutputArg;
        this.showDefaultsAndExamples = showDefaultsAndExamplesArg;
        this.outputAlleleAge = outputAlleleAgeArg;
        this.seqOutFormat = seqOutFormatArg;
        this.hapsOutFormat = hapsOutFormatArg;
        this.outputTrees = outputTreesArg;
        this.recRate = recRateArg;
        this.MAF = MAFArg;
        this.mutationRate = mutationRateArg;
        this.recRateFactor = recRateFactorArg;
        this.IBDthresholdCM = IBDthresholdCMArg;
        this.geneConvRelativeRate = geneConvRelativeRateArg;
        this.probGeneConversionGivenRec = probGeneConversionGivenRecArg;
        this.popSizeFile = popSizeFileArg;
        this.mapFile = mapFileArg;
        this.fileBase = fileBaseArg;
        this.numPopsAtZero = numPopsArg;
        this.sampleSizes = sampleSizesArg;
        this.genomeSizePhysical = genomeSizePhysicalArg;
        this.minBlockSize = minBlockSizeArg;
        this.totModernSamples = totModernSamplesArg;
        this.meanGeneConversionTractLength = meanGeneConversionTractLengthArg;
        this.minimumGeneConversionTractLength = minimumGeneConversionTractLengthArg;
        this.constantPopSize = constantPopSizeArg;
        this.seed = seedArg;

        if (shrinkSequenceData && !seqOutFormat) {
            exit("\"-shrink\" can only be used with -seq-out.");
        }

        int totSamples = 0;
        for (int i = 0; i < sampleSizes.length; i++) {
            totSamples += sampleSizes[i];
            if (sampleSizes[i] <= 0) {
                exit("Sample sizes must be > 0.");
            }
        }
        DescendantsList.numSamples = totSamples;

        if (useMapFile) {
            recRateMap = new TreeMap<Integer, Float>();
            physToGen = new TreeMap<Integer, Integer>();
            mutRateMap = new TreeMap<Integer, Float>();
            recToPhysicalMap = new TreeMap<Integer, ArrayList<Integer>>();
            try {
                BufferedReader br = null;
                try {
                    br = new BufferedReader(new FileReader(mapFile));
                } catch (FileNotFoundException ex) {
                    exit("Could not open mutation rates file " + mapFile);
                }
                try {
                    String line = br.readLine();
                    int pos = 0;
                    while (line != null) {
                        String[] strSplit = line.split("\\s+");
                        if (strSplit.length != 4) {
                            exit("Error parsing line \"" + line + "\" in map file. Use format\n\tposition1\trecombinationRate1\tcumulativeRecombinationRate1\tmutationRate1\n\tposition1\trecombinationRate1\tcumulativeRecombinationRate1\tmutationRate1\n\t...");
                        }
                        pos = Integer.parseInt(strSplit[0]);
                        float recRateBin = Float.parseFloat(strSplit[1]);
                        int genPosBin = (int) Math.round(Float.parseFloat(strSplit[2]) * 1000000.f);
                        if (genPosBin == 0) {
                            genPosBin++;
                        }
                        float mutRateBin = Float.parseFloat(strSplit[3]);
                        recRateMap.put(pos, recRateBin);
                        physToGen.put(pos, genPosBin);
                        mutRateMap.put(pos, mutRateBin);
                        ArrayList<Integer> p;
                        if (!recToPhysicalMap.containsKey(genPosBin)) {
                            p = new ArrayList<Integer>();
                        } else {
                            p = recToPhysicalMap.get(genPosBin);
                        }
                        p.add(pos);
                        recToPhysicalMap.put(genPosBin, p);
                        line = br.readLine();
                    }
                    if (pos < genomeSizePhysical) {
                        exit("The map is shorter than the required genome size.");
                    }
                    br.close();
                    System.out.println("Parsed " + recRateMap.size() + " entries of the map. "
                            + "First position is " + recRateMap.firstKey() + ". Last position is " + recRateMap.lastKey());
                    int firstKey = recRateMap.firstKey();
                    if (firstKey != 1) {
                        if (firstKey < 1) {
                            exit("Map starts at value < 1.");
                        }
                        float firstRecPos = physToGen.get(firstKey);
                        float rRate = firstRecPos / (firstKey / 1000000.f); // in cM/Mb
                        float firstMutVal = mutRateMap.get(firstKey);
                        recRateMap.put(1, rRate);
                        mutRateMap.put(1, firstMutVal);
                        physToGen.put(1, 1);
                        ArrayList<Integer> p;
                        if (!recToPhysicalMap.containsKey(1)) {
                            p = new ArrayList<Integer>();
                        } else {
                            p = recToPhysicalMap.get(1);
                        }
                        p.add(1);
                        recToPhysicalMap.put(1, p);
                        System.out.println("Warning: recombination map does not start at 1. First entry inferred to be "
                                + "1\t" + rRate + "\t1\t" + firstMutVal);
                    }
                    Map.Entry<Integer, Integer> entry = physToGen.floorEntry(genomeSizePhysical);
                    int lastPhys = entry.getKey();
                    int lastGenPos = entry.getValue();
                    if (physToGen.lastEntry().getKey().equals(entry.getKey())) {
                        float lastRecRate = recRateMap.lastEntry().getValue();
                        System.out.println(lastRecRate);
                        genomeSizeGenetic = lastGenPos + (int) Math.round((genomeSizePhysical - lastPhys) * lastRecRate);
                        System.out.println("The specified map is shorter than the chosen physical length. Genome size in cM inferred to be " + (genomeSizeGenetic / 1000000.) + "\n");
                    } else if (lastPhys == genomeSizePhysical) {
                        genomeSizeGenetic = lastGenPos;
                    } else {
                        Map.Entry<Integer, Integer> entryTo = physToGen.ceilingEntry(genomeSizePhysical);
                        int nextPhys = entryTo.getKey();
                        int nextGenPos = entryTo.getValue();
                        float fractionPhysical = (genomeSizePhysical - lastPhys) / (nextPhys - lastPhys);
                        genomeSizeGenetic = lastGenPos + (int) Math.round(fractionPhysical * (nextGenPos - lastGenPos));
                        System.out.println("The specified map is longer than the chosen physical length. Genome size in cM inferred to be " + (genomeSizeGenetic / 1000000.) + "\n");
                    }
                } catch (IOException ex) {
                    exit("Could not read map file " + mapFile);
                }
            } catch (Exception e) {
                exit("Error parsing -map parameter.");
            }
        } else {
            if (genomeSizePhysical * recRate * 100000000.f > 2146000000.) {
                exit("Region genetic size cannot exceed " + 2146. + " cM");
            }
            genomeSizeGenetic = (int) Math.round(genomeSizePhysical * recRate * 100000000.f);
            recRateMap = null;
            physToGen = null;
            mutRateMap = null;
            recToPhysicalMap = null;
        }

        if (meanGeneConversionTractLength < minimumGeneConversionTractLength) {
            exit("The mean length of gene conversion tracts cannot be smaller than the minimum length.");
        }
        GCPhysicalLengthRateForSampling = 1.f / (meanGeneConversionTractLength - minimumGeneConversionTractLength);
        GCGeneticLengthRateForSampling = GCPhysicalLengthRateForSampling / recRateFactor;
        minimumGeneConversionTractLengthGenetic = (int) (minimumGeneConversionTractLength * recRateFactor);

        thresholdForPrintingBooleanRatherThanList = totSamples / (1 + (int) Math.ceil(Math.log10(totSamples)));
        System.out.println("Starting simulation. Parameters:");
        printParams(showDefaultsAndExamples);

        try {
            demographicModel = (usePopFile) ? new Population(popSizeFile, quiet) : new Population(constantPopSize);
            if (demographicModel.getNumPopsAt(0) != numPopsAtZero) {
                exit("User specified " + numPopsAtZero + " population(s), but " + demographicModel.populationSizes.get(0).size() + " are present in the demographic model at time 0.");
                System.exit(0);
            }
            for (int i = 0; i < demographicModel.getNumPopsAt(0); i++) {
                if (sampleSizes[i] > demographicModel.populationSizes.get(0).get(i)) {
                    exit("Population " + (i + 1) + " has size " + demographicModel.populationSizes.get(0).get(i)
                            + ", but user requested " + sampleSizes[i] + " samples.");
                    System.exit(0);
                }
            }
        } catch (Exception ex) {
            Logger.getLogger(Argon.class.getName()).log(Level.SEVERE, null, ex);
        }
        leaves = new ArrayList<Block>(totModernSamples);
        if (printIBD) {
            IBDFactory = new IBD(IBDthresholdCM, totModernSamples);
        } else {
            IBDFactory = null;
        }

    }

    public void simulate() {
        ArrayList<HashMap<Long, Individual>> populationsChildren = new ArrayList<HashMap<Long, Individual>>();
        for (int i = 0; i < demographicModel.getNumPopsAt(0); i++) {
            populationsChildren.add(new HashMap<Long, Individual>());
        }
        ArrayList<HashMap<Long, Individual>> populationsParents = new ArrayList<HashMap<Long, Individual>>();
        for (int i = 0; i < demographicModel.getNumPopsAt(1); i++) {
            populationsParents.add(new HashMap<Long, Individual>());
        }

        leaves = new ArrayList<Block>(totModernSamples);
        int sampleCnt = 0;
        for (int i = 0; i < numPopsAtZero; i++) {
            HashMap<Long, Individual> pop = populationsChildren.get(i);
            for (int j = 1; j <= sampleSizes[i]; j++) {
                Individual id = new Individual(j, genomeSizeGenetic);
                sampleCnt++;
                id.segments.get(0).ID = sampleCnt;
                id.segments.get(0).dList = new DescendantsList(sampleCnt);
                leaves.add(id.segments.get(0));
                pop.put((long) sampleCnt, id);
            }
            populationsChildren.set(i, pop);
        }

        int completed = 0;
        generator = new Random(seed);
        roots = new ArrayList<Block>();

        int gen = 0; // current generation
        long start = System.currentTimeMillis(), time = start;

        int nextRec;
        Individual nextIndividualChunk, parent;

        long toCoalesce = (totModernSamples - 1) * (long) genomeSizeGenetic;
        long coalescedBlocks = 0;

        long lastTime = -1l;

        System.out.println("Simulating...");
        while (coalescedBlocks != toCoalesce) {
            HashMap<Block, HashSet<Block>> childrenOfThisGeneration = new HashMap<Block, HashSet<Block>>();
            if (debug || !quiet && (long) ((System.currentTimeMillis() - start) / 1000.0) > lastTime) {
                lastTime = (long) Math.floor((System.currentTimeMillis() - start) / 1000.0);
                completed = genomeSizeGenetic - numBlocksLeftToTrack;
                DecimalFormat df = new DecimalFormat("#.##");
                int trackedInd = 0;
                StringBuilder trackString = new StringBuilder();
                for (int i = 0; i < populationsChildren.size(); i++) {
                    trackedInd += populationsChildren.get(i).size();
                    trackString.append(" ");
                    trackString.append(populationsChildren.get(i).size());
                }
//                if (verbose) {
//                    System.out.println("\n");
//                }
                System.out.print("\rGENERATION: " + gen + " tracked individuals: " + trackedInd + " [" + trackString + " ]"
                        + " elapsed time: " + (int) ((System.currentTimeMillis() - start) / 1000.0) + " seconds           ");
                time = System.currentTimeMillis();
            }

            // sample a parent for all current individuals we're following
            int numPops = demographicModel.getNumPopsAt(gen);
            double nextRecOrGCRate = (doGeneConversion) ? 1e-8f * (1.f + geneConvRelativeRate) : 1e-8f;
            for (int currPop = 0; currPop < numPops; currPop++) {
                HashMap<Long, Individual> individualsInCurrentPop = populationsChildren.get(currPop);
//                if (verbose) {
//                    System.out.println("Starting samples from population " + currPop + " of size " + individualsInCurrentPop.size());
//                }
                for (Object el : individualsInCurrentPop.keySet()) {
                    Individual ind = (Individual) individualsInCurrentPop.get(el);
                    LinkedList<Individual> parentQueue = new LinkedList<Individual>();
                    LinkedList<Integer> sampledPopQueue = new LinkedList<Integer>();
//                    if (verbose) {
//                        System.out.println();
//                    }
                    int GCend = -1;
                    while (ind != null) {
                        if (doGeneConversion) { //TODO: should avoid this and try to have unique piece of code for both gene conversion and no-gene conversion
//                            if (verbose) {
//                                System.out.println("Pop: " + currPop + " Ind: " + ind.ID + "; gen " + gen);
//                                ind.print();
//                            }
                            nextRec = 0;
                            if (ind.segments.get(0).start > GCend) {
                                // the previous tract was not gene conversion, or the gene conversion tract has already ended by the start of the current block
                                int distToNextRec = (int) Math.round(Tools.sampleExponential(nextRecOrGCRate, generator));
                                nextRec = ind.segments.get(0).start + Tools.roundToBlock(distToNextRec, minBlockSize);
//                                if (verbose) {
//                                    System.out.println("Not in GC tract. Sampled\t" + distToNextRec / 1000000. + ". Current start is: " + ind.segments.get(0).start / 1000000. + " (GCend is " + GCend / 1000000. + ")");
//                                }
                            } else {
                                // the previous tract was gene conversion, and the gene conversion tract has not ended by the start of the current block
                                nextRec = GCend;
//                                if (verbose) {
//                                    System.out.println("In GC tract. Will terminate at previously sampled position:\t" + nextRec);
//                                }
                            }
                            int parentTill = nextRec;
                            if (ind.segments.get(ind.segments.size() - 1).end >= parentTill && parentTill >= ind.segments.get(0).start) {
//                                if (verbose) {
//                                    System.out.println("Will recombine at " + parentTill);
//                                }
                                // recombination will happen, handle current individual and update current with remaining part, then iterate
                                ArrayList<Block> otherGuySegments = ind.splitAfterSomeBlocks(parentTill, ind.ID, gen);
                                nextIndividualChunk = new Individual(ind.ID, otherGuySegments);
                                nextIndividualChunk.ID = ind.ID;
                            } else {
                                nextIndividualChunk = null;
                            }

                            //handle current individual here
                            long parentID;
                            int sampledPopID;
                            HashMap<Long, Individual> sampledPop;
                            if (parentQueue.size() < 2) {
                                //sample a population
                                double[][] migMatrix = demographicModel.getSamplingProbabilities(gen);
                                double randomNumber = generator.nextDouble();
                                double totProb = 0.;
                                sampledPopID = 0;
//                                if (verbose) {
//                                    System.out.println("Sampling parent for population " + currPop + " Ind: " + ind.ID);
//                                }
                                // TODO: speed this up by storing non-zero prob vector in population class, iterate on that instead.
                                for (int i = 0; i < migMatrix[currPop].length; i++) {
//                                    if (verbose) {
//                                        System.out.println("\tsampling " + randomNumber + " " + (totProb + migMatrix[currPop][i]));
//                                    }
                                    if (totProb + migMatrix[currPop][i] >= randomNumber) {
                                        sampledPopID = i;
                                        break;
                                    }
                                    totProb += migMatrix[currPop][i];
                                }
                                sampledPop = populationsParents.get(sampledPopID);

                                //sample a parent
                                parentID = demographicModel.getIDOffsetAt(gen + 1, sampledPopID) + (long) Math.floor(generator.nextDouble() * (demographicModel.getSizeAt(gen + 1, sampledPopID))) + 1;
                                parent = sampledPop.get(parentID);
                            } else {
                                parent = parentQueue.poll();
                                sampledPopID = sampledPopQueue.poll();
                                sampledPop = populationsParents.get(sampledPopID);
                                parentID = parent.ID;
                            }
//                            if (verbose) {
//                                System.out.println("Sampled parent from pop: " + sampledPopID + " parentID " + parentID + " offset " + demographicModel.getIDOffsetAt(gen + 1, sampledPopID));
//                            }
                            //handle coalescence
                            if (parent == null) {
                                parent = new Individual(parentID, ind.segments);
                            } else {
                                ArrayList<Block> coalesced = parent.mergeArrayLists(ind.segments, gen + 1, parentID, totModernSamples, roots, childrenOfThisGeneration);
                                for (int i = 0; i < coalesced.size(); i++) {
                                    coalescedBlocks += coalesced.get(i).end - coalesced.get(i).start + 1;
                                }
                            }
                            if (!parent.segments.isEmpty()) {
                                sampledPop.put(parentID, parent);
                            } else {
                                sampledPop.remove(parentID);
                            } //done
                            ind = nextIndividualChunk;
                            parentQueue.add(parent);
                            sampledPopQueue.add(sampledPopID);
                            // if GCend == nextRec, I have just terminated a GC tract, will sample next GC/rec event. Else, may sample a GC
                            if (GCend != nextRec && parentTill < genomeSizeGenetic && generator.nextDouble() < probGeneConversionGivenRec) {
                                // enter gene conversion tract
                                GCend = Tools.roundToBlock(sampleGeneConversionTract(parentTill, debug), minBlockSize);;
//                                if (verbose) {
//                                    System.out.println("Gene conversion tract from " + parentTill + " to " + GCend);
//                                }
                            } else {
                                GCend = -1;
//                                if (verbose) {
//                                    if (GCend != nextRec) {
//                                        if (verbose) {
//                                            System.out.println("No gene conversion tract, standard recombination will happen (because of sampling or genome size).");
//                                        }
//                                    } else {
//                                        System.out.println("No gene conversion tract, standard recombination will happen (because a GC just happened).");
//                                    }
//                                }
                            }
                        } else {
                            if (debug) {
                                System.out.println("Pop: " + currPop + " Ind: " + ind.ID + "; gen " + gen);
                            }
                            if (debug) {
                                ind.print();
                            }
                            nextRec = 0;
                            while (nextRec == 0) {
                                nextRec = (int) Math.round(Tools.sampleExponential(1e-8f, generator));
                                nextRec = Tools.roundToBlock(nextRec, minBlockSize);
                            }
                            int parentTill = ind.segments.get(0).start + nextRec;
                            if (ind.segments.get(ind.segments.size() - 1).end >= parentTill) {
                                if (debug) {
                                    System.out.println("Will recombine at " + parentTill);
                                }
                                // recombination will happen, handle current individual and update current with remaining part, then iterate
                                ArrayList<Block> otherGuySegments = ind.splitAfterSomeBlocks(parentTill, ind.ID, gen);
                                nextIndividualChunk = new Individual(ind.ID, otherGuySegments);
                                nextIndividualChunk.ID = ind.ID;
                            } else {
                                nextIndividualChunk = null;
                            }

                            //handle current individual here
                            //sample a population
                            double[][] migMatrix = demographicModel.getSamplingProbabilities(gen);
                            double randomNumber = generator.nextDouble();
                            double totProb = 0.;
                            int sampledPopID = 0;
                            if (debug) {
                                System.out.println("Sampling parent for population " + currPop + " Ind: " + ind.ID);
                            }
                            for (int i = 0; i < migMatrix[currPop].length; i++) {
                                if (debug) {
                                    System.out.println("\tsampling " + randomNumber + " " + (totProb + migMatrix[currPop][i]));
                                }
                                if (totProb + migMatrix[currPop][i] >= randomNumber) {
                                    sampledPopID = i;
                                    break;
                                }
                                totProb += migMatrix[currPop][i];
                            }
                            HashMap<Long, Individual> sampledPop = populationsParents.get(sampledPopID);

                            //sample a parent
                            long parentID = demographicModel.getIDOffsetAt(gen + 1, sampledPopID) + (long) Math.floor(generator.nextDouble() * (demographicModel.getSizeAt(gen + 1, sampledPopID))) + 1;
                            parent = sampledPop.get(parentID);
                            //handle coalescence
                            if (parent == null) {
                                parent = new Individual(parentID, ind.segments);
                            } else {
                                ArrayList<Block> coalesced = parent.mergeArrayLists(ind.segments, gen + 1, parentID, totModernSamples, roots, childrenOfThisGeneration);
                                for (int i = 0; i < coalesced.size(); i++) {
                                    coalescedBlocks += coalesced.get(i).end - coalesced.get(i).start + 1;
                                }
                            }
                            if (debug) {
                                System.out.println("\tsampled from pop: " + sampledPopID + " parentID " + parentID + " offset " + demographicModel.getIDOffsetAt(gen + 1, sampledPopID) + "\n");
                            }
                            if (parent.segments.size() != 0) {
                                sampledPop.put(parentID, parent);
                            } else {
                                sampledPop.remove(parentID);
                            }
                            ind = nextIndividualChunk;
                        }
                    }
                }
            } // done for each pop
            // assign parents for next generation
            populationsChildren = populationsParents;
            gen++;
            populationsParents = new ArrayList<HashMap<Long, Individual>>();
            for (int i = 0; i < demographicModel.getNumPopsAt(gen + 1); i++) {
                populationsParents.add(new HashMap<Long, Individual>());
            }
        }

        Collections.sort(roots);

        long timeB = System.currentTimeMillis();
        // insert breakpoints into maps
        if (useMapFile) {
            System.out.print("\nAdding breakpoints into map ... ");
            for (Block b : roots) {
                addGenPosToMaps(b.start);
            }
            addGenPosToMaps(genomeSizeGenetic + 1);
            int cnt = 0;
            for (Integer genPos : recToPhysicalMap.keySet()) {
                ArrayList<Integer> aPhys = recToPhysicalMap.get(genPos);
                Collections.sort(aPhys);
                for (Integer physPos : aPhys) {
                    physicalPosMapArray.add(physPos);
                    mutRateArray.add(mutRateMap.get(physPos));
                    genPosArray.add(genPos);
                    recRateArray.add(recRateMap.get(physPos));
                    mapIndex.put(genPos, cnt);
                    cnt++;
                }
            }
            recRateMap.clear();
            mutRateMap.clear();
            physToGen.clear();
            System.out.println((System.currentTimeMillis() - timeB) / 1000.0 + " seconds.");
        }

        System.out.print("\n" + roots.size() + " blocks\t");

        if (writeOutFiles) {
            String name = (gzipOutput) ? fileBase + ".map.gz" : fileBase + ".map";
            try {
                MAPoutput = new FileOutputStream(name);
            } catch (FileNotFoundException e) {
                exit("file " + name + " could not be created.");
            }
            try {
                MAPwriter = (gzipOutput) ? new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(MAPoutput), "UTF-8"), bufferSize) : new BufferedWriter(new OutputStreamWriter(MAPoutput), bufferSize);
            } catch (IOException e) {
                exit("could not open output stream for " + name);
            }
            if (outputAlleleAge) {
                name = (gzipOutput) ? fileBase + ".age.gz" : fileBase + ".age";
                try {
                    AGEoutput = new FileOutputStream(name);
                } catch (FileNotFoundException e) {
                    exit("file " + name + " could not be created.");
                }
                try {
                    AGEwriter = (gzipOutput) ? new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(AGEoutput), "UTF-8"), bufferSize) : new BufferedWriter(new OutputStreamWriter(AGEoutput), bufferSize);
                } catch (IOException e) {
                    exit("could not open output stream for " + name);
                }
            }
            if (outputTrees) {
                name = (gzipOutput) ? fileBase + ".trees.gz" : fileBase + ".trees";
                try {
                    treeOutput = new FileOutputStream(name);
                } catch (FileNotFoundException e) {
                    exit("file " + name + " could not be created.");
                }
                try {
                    treeWriter = (gzipOutput) ? new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(treeOutput), "UTF-8"), bufferSize) : new BufferedWriter(new OutputStreamWriter(treeOutput), bufferSize);
                } catch (IOException e) {
                    exit("could not open output stream for " + name);
                }
            }
            if (printMUT) {
                if (seqOutFormat) {
                    name = (gzipOutput) ? fileBase + ".mut.gz" : fileBase + ".mut";
                    try {
                        MUToutput = new FileOutputStream(name);
                    } catch (FileNotFoundException e) {
                        exit("file " + name + " could not be created.");
                    }
                    try {
                        MUTwriter = (gzipOutput) ? new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(MUToutput), "UTF-8"), bufferSize) : new BufferedWriter(new OutputStreamWriter(MUToutput), bufferSize);
                    } catch (IOException e) {
                        exit("could not open output stream for " + name);
                    }
                } else if (hapsOutFormat) {
                    name = (gzipOutput) ? fileBase + ".hap.gz" : fileBase + ".hap";
                    try {
                        MUToutput = new FileOutputStream(name);
                    } catch (FileNotFoundException e) {
                        exit("file " + name + " could not be created.");
                    }
                    try {
                        MUTwriter = (gzipOutput) ? new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(MUToutput), "UTF-8"), bufferSize) : new BufferedWriter(new OutputStreamWriter(MUToutput), bufferSize);
                    } catch (IOException e) {
                        exit("could not open output stream for " + name);
                    }
                } else {
                    name = (gzipOutput) ? fileBase + ".vcf.gz" : fileBase + ".vcf";
                    try {
                        MUToutput = new FileOutputStream(name);
                    } catch (FileNotFoundException e) {
                        exit("file " + name + " could not be created.");
                    }
                    try {
                        MUTwriter = (gzipOutput) ? new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(MUToutput), "UTF-8"), bufferSize) : new BufferedWriter(new OutputStreamWriter(MUToutput), bufferSize);
                    } catch (IOException e) {
                        exit("could not open output stream for " + name);
                    }
                }
            }
            if (printIBD && writeOutFiles) {
                try {
                    IBDoutput = (gzipOutput) ? new FileOutputStream(fileBase + ".ibd.gz") : new FileOutputStream(fileBase + ".ibd");
                } catch (FileNotFoundException e) {
                    exit("file " + fileBase + ".mut.gz could not be created.");
                }
                try {
                    IBDwriter = (gzipOutput) ? new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(IBDoutput), "UTF-8"), bufferSize) : new BufferedWriter(new OutputStreamWriter(IBDoutput), bufferSize);
                } catch (IOException e) {
                    exit("could not open output stream for IBD output");
                }
            }
            try {
                MAPwriter.write("BREAKPOINTS");
                if (!useMapFile) {
                    float recRateFactor = recRate * 100000000.0f;
                    for (Block root : roots) {
                        MAPwriter.write("\t" + ((int) Math.round(root.start / recRateFactor)));
                    }
                } else {
                    for (int i : physicalPosMapArray) {
                        int limit = recToPhysicalMap.get(genomeSizeGenetic).get(recToPhysicalMap.get(genomeSizeGenetic).size() - 1);
                        if (i <= limit) {
                            MAPwriter.write("\t" + i);
                        }
                    }
                }
                MAPwriter.write("\n");
            } catch (IOException ex) {
                Logger.getLogger(Argon.class.getName()).log(Level.SEVERE, null, ex);
            }
        } else {
            if (printMUT) {
                System.out.print("\nBREAKPOINTS");
                if (!useMapFile) {
                    float recRateFactor = recRate * 100000000.0f;
                    for (Block root : roots) {
                        System.out.print("\t" + ((int) Math.round(root.start / recRateFactor)));
                    }
                } else {
                    for (int i : physicalPosMapArray) {
                        System.out.print("\t" + i);
                    }
                }
                System.out.println();
            }
        }

        ARG_runtime = (System.currentTimeMillis() - start) / 1000.0;

        System.out.println("ARG runtime: " + ARG_runtime);

    }

    public void visit() throws Exception {
        long start = System.currentTimeMillis(), time = start;

        // for IBD printing
        int IBD_from = 0;
        int IBD_to = 0;
        int IBD_fromGen = 0;
        int IBD_toGen = 0;
        ArrayList<DescendantsList> IBD_sets = null;

        if (printMUT || printIBD || writeAFS || printAFS) { //TODO: efficient visit of ARG
            Comparator<Block> BGC = new BlockGenComparator();
            DescendantsList.numSamples = totModernSamples;
            // use to traverse bottom up
            PriorityQueue<Block> blockQueue = new PriorityQueue<Block>(leaves.size(), BGC);
            // contains list of descendants
            HashMap<Block, DescendantsList> currentBlockMap = new HashMap<Block, DescendantsList>();
            HashSet<Block> activeBlocks = new HashSet<Block>(); // create list of those from previous gen to be removed
            for (int i = 1; i <= leaves.size(); i++) {
                currentBlockMap.put(leaves.get(i - 1), new DescendantsList((int) leaves.get(i - 1).ID));
            }

            HashMap<Block, Integer> dListRequests = new HashMap<Block, Integer>();
            HashMap<Block, Integer> dListUsage = new HashMap<Block, Integer>();
            int cnt = 0, oldPercentage = 0;

            if (writeOutFiles) {
                System.out.println();
            }
            if (!(seqOutFormat || hapsOutFormat)) {
                printVCFheader();
            }
            if (hapsOutFormat) {
                printSamplesFile();
            }

            if (outputTrees) {
                for (Block currentRoot : roots) {
                    StringBuilder treeString = new StringBuilder(currentRoot.gen + "\t" + currentRoot.start + "\t" + currentRoot.end + "\t");
                    writeNewickTree(treeString, currentRoot, 0, currentRoot.start, currentRoot.end);
                    if (!writeOutFiles) {
                        System.out.println("Tree\t" + treeString + ";");
                    } else {
                        treeWriter.write(treeString + ";\n");
                    }
                }
            }

            for (Block currentRoot : roots) {
                if (writeOutFiles || (!printMUT && (writeAFS || printAFS))) {
//                    System.out.println();
                    cnt++;
                    int percentage = 100 * cnt / roots.size();
                    if (!quiet && (percentage > oldPercentage || percentage == 0)) {
                        oldPercentage = percentage;
                        System.out.print("\rVisiting ARG of " + totModernSamples + " samples, ~" + genomeSizePhysical / 1000000 + " Mb to compute output\t" + percentage + "%");
                    }
                }

                LinkedList<Block> topDownList = new LinkedList<Block>();
                topDownList.add(currentRoot);
                blockQueue.add(currentRoot);
                activeBlocks.add(currentRoot);
                while (true) {
                    Block current = topDownList.poll();
                    for (Block child : current.children) {
                        if (!activeBlocks.contains(child)) {
                            blockQueue.add(child);
                            topDownList.add(child);
                            activeBlocks.add(child);
                        } else {
                        }
                    }
                    if (topDownList.isEmpty()) {
                        break;
                    }
                }
                HashSet<Block> doneThisGen = new HashSet();
                while (!blockQueue.isEmpty()) {
                    Block currentBlock = blockQueue.poll();
                    //do block
                    DescendantsList myList = currentBlockMap.get(currentBlock);
                    if (myList == null) {
                        myList = new DescendantsList();
                    }
                    for (Block child : currentBlock.children) {
                        DescendantsList childsList = currentBlockMap.get(child);
                        if (childsList == null) {
                            childsList = new DescendantsList();
                        }
                        myList.add(childsList);
                        Integer requests; // get number of requests on child
                        if (dListRequests.containsKey(child)) {
                            requests = dListRequests.get(child);
                        } else {
                            requests = 0;
                        }
                        requests++;
                        dListRequests.put(child, requests);
                        if (child.parents.size() == requests) { //if requests to child are done from each parent, safe to remove
                            doneThisGen.add(child);
                        }
                    }
                    currentBlockMap.put(currentBlock, myList);
                    if (printMUT || writeAFS || printAFS) {
                        for (Block parentBlock : currentBlock.parents) {
                            int size = myList.getSize();
                            float freq = (float) size / totModernSamples;
                            int distanceToFather = parentBlock.gen - currentBlock.gen;
                            int fromPhys, toPhys;
                            if (!useMapFile) {
                                fromPhys = (int) Math.round(parentBlock.start / recRateFactor);
                                toPhys = (int) Math.round(parentBlock.end / recRateFactor);
                                outputMut(fromPhys, toPhys, mutationRate, currentBlock.gen, parentBlock.gen, distanceToFather, size, freq, myList);
                            } else {
                                int recMutblockLastIndex = mapIndex.get(parentBlock.end + 1);
                                for (int recMutblockIndex = mapIndex.get(parentBlock.start); recMutblockIndex < recMutblockLastIndex; recMutblockIndex++) {
                                    fromPhys = physicalPosMapArray.get(recMutblockIndex);
                                    toPhys = physicalPosMapArray.get(recMutblockIndex + 1) - 1;
                                    float mut = mutRateArray.get(recMutblockIndex);
                                    outputMut(fromPhys, toPhys, mut, currentBlock.gen, parentBlock.gen, distanceToFather, size, freq, myList);
                                }
                            }
                        }
                    }
                }
                doneThisGen.add(currentRoot);
                for (Block doneBlock : doneThisGen) {
                    if (doneBlock.gen != 0) {
                        if (printIBD) {
                            IBD_sets = new ArrayList<DescendantsList>();
                            int fromPhys, toPhys;
                            if (!useMapFile) {
                                fromPhys = (int) Math.round(doneBlock.start / recRateFactor);
                                toPhys = (int) Math.round(doneBlock.end / recRateFactor);
                            } else {
                                int recMutblockIndexFrom = mapIndex.get(doneBlock.start);
                                fromPhys = physicalPosMapArray.get(recMutblockIndexFrom);
                                int recMutblockIndexTo = mapIndex.get(doneBlock.end + 1);
                                toPhys = physicalPosMapArray.get(recMutblockIndexTo) - 1;
                            }
                            IBD_from = fromPhys;
                            IBD_to = toPhys;
                            IBD_fromGen = doneBlock.start;
                            IBD_toGen = doneBlock.end;
                        }
                        for (Block c : doneBlock.children) {
                            DescendantsList dL = currentBlockMap.get(c);
                            Integer uses; // get number of requests on child
                            if (dListUsage.containsKey(c)) {
                                uses = dListUsage.get(c);
                            } else {
                                uses = 0;
                            }
                            uses++;
                            dListUsage.put(c, uses);
                            if (c.parents.size() == uses) { //if requests to child are done from each parent, safe to remove
                                dListRequests.remove(c);
                                dListUsage.remove(c);
                                currentBlockMap.remove(c);
                            }
                            if (printIBD) {
                                IBD_sets.add(dL);
                            }
                        }
                        if (printIBD) {
                            int IBD_gen = doneBlock.gen;
                            long IBD_ID = doneBlock.ID;
                            try {
                                StringBuilder IBDstring = IBDFactory.addIBDset(IBD_from, IBD_to, IBD_fromGen, IBD_toGen, IBD_gen, IBD_ID, IBD_sets);
                                if (IBDstring.length() != 0) {
                                    if (writeOutFiles) {
                                        IBDwriter.write(IBDstring.toString());
                                    } else {
                                        System.out.print("IBD\t" + IBDstring.toString());
                                    }
                                }
                            } catch (Exception e) {
                                System.err.println("Error computing IBD segments" + e.toString());
                                throw e;
                            }
                        }
                    }
                    activeBlocks.remove(doneBlock);
                }
            }
        }

        if (printIBD) {
            try {
                StringBuilder IBDstring = IBDFactory.finalize(genomeSizeGenetic);
                if (writeOutFiles) {
                    IBDwriter.write(IBDstring.toString());
                } else {
                    System.out.print("IBD\t" + IBDstring.toString());
                }
            } catch (Exception e) {
                exit("Error computing IBD segments: " + e.toString());
            }
        }
        if (writeOutFiles) {
            try {
                MAPwriter.flush();
                MAPwriter.close();
                if (outputAlleleAge) {
                    AGEwriter.flush();
                    AGEwriter.close();
                }
                if (printMUT) {
                    MUTwriter.flush();
                    MUTwriter.close();
                }
                if (outputTrees) {
                    treeWriter.flush();
                    treeWriter.close();
                }
                if (printIBD) {
                    IBDwriter.flush();
                    IBDwriter.close();
                }
            } catch (IOException ex) {
                Logger.getLogger(Argon.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        if (writeAFS) {
            try {
                PrintWriter writer = new PrintWriter(fileBase + ".afs", "UTF-8");
                for (int i = 1; i < AFS.length; i++) {
                    writer.println(i + "\t" + AFS[i] + "\t" + AFS[i] / totMut);
                }
                writer.close();
            } catch (UnsupportedEncodingException e) {
                exit("file " + fileBase + ".afs could not be created.");
            } catch (FileNotFoundException ex) {
                exit("file " + fileBase + ".afs could not be created.");
            }
        }
        if (printAFS) {
            System.out.println();
            for (int i = 1; i < AFS.length; i++) {
                System.out.println("AFS\t" + i + "\t" + AFS[i] + "\t" + AFS[i] / totMut);
            }
        }

        if (writeAFS || printAFS) {
            AFS = new float[totModernSamples];
            if (!(printMUT || printIBD)) {
                // not interested in sequence, just compute AFS
                DescendantsList.onlyStoreNumber = true;
            }
        } else {
            AFS = null;
        }

        Visit_runtime = (System.currentTimeMillis() - start) / 1000.0;
        System.out.println("\nARG visit runtime: " + Visit_runtime);

    }

    int sampleGeneConversionTract(int fromGenPos, boolean verbose) {
        int nextRec;
        if (useMapFile) {
            // find current physical position
            int phys = genToPhys(fromGenPos);
            // find next physical position
            int nextPhys = phys + minimumGeneConversionTractLength + (int) Math.round(Tools.sampleExponential(GCPhysicalLengthRateForSampling, generator));
            // convert it to genetic position
            Map.Entry<Integer, Integer> entry = physToGen.floorEntry(nextPhys);
            int physFloor = entry.getKey();
            int genFloor = entry.getValue();
            entry = physToGen.ceilingEntry(nextPhys);
            int physCeil = entry.getKey();
            int genCeil = entry.getValue();
            float recBlockPhysLen = (physCeil - physFloor);
            float recBlockGenLen = (genCeil - genFloor);
            // find current location in genetic map coordinates
            if (recBlockPhysLen != 0f) { // TODO: may not need this, or this may not be enough to avoid problems
                float frac = (nextPhys - physFloor) / recBlockPhysLen;
                nextRec = Tools.roundToBlock((int) Math.round(genFloor + frac * recBlockGenLen), minBlockSize);
            } else {
                nextRec = Tools.roundToBlock(fromGenPos + genCeil, minBlockSize);;
            }
//            if (verbose) {
//                System.out.println("Sampling GC tract using genetic map. Current start: " + fromGenPos + " to " + nextRec + " genLen " + (nextRec - fromGenPos) + " physLen " + (nextPhys - phys));
//            }
        } else {
            int recDist = minimumGeneConversionTractLengthGenetic + (int) Math.round(Tools.sampleExponential(GCGeneticLengthRateForSampling, generator));
//            if (verbose) {
//                System.out.println("Sampling GC tract not using map. Sampled\t" + recDist);
//            }
            nextRec = Tools.roundToBlock(fromGenPos + recDist, minBlockSize);
        }
        return nextRec;
    }

    void printParams(boolean showDefaultsAndExamples) {
        // TODO: avoid hard-coding default values.
        System.out.println("\nParameter values:");
        System.out.println("\t-N\t\tPopulation size is " + ((popSizeFile != "") ? popSizeFile : constantPopSize) + ".");
        if (showDefaultsAndExamples) {
            System.out.println("\t\t\t(Default = 1000; example: \"-N 10000\")");
        }
        System.out.println("\t-pop\t\tThere will be " + numPopsAtZero + " population(s).");
        for (int i = 0; i < numPopsAtZero; i++) {
            System.out.println("\t\t\tPopulation " + (i + 1) + " will have " + sampleSizes[i] + " samples.");
        }
        if (showDefaultsAndExamples) {
            System.out.println("\t\t\t(Default = 1 pop, 1000 samples; example: \"-pop 1 1000\" or \"-pop 2 1000 2000\")");
        }
        System.out.println("\t-len\t\tMinimum block size is " + minBlockSize + " microMorgans.");
        if (showDefaultsAndExamples) {
            System.out.println("\t\t\t(Default = 1 microMorgan; example: \"-len 1\" or \"-len 10000\")");
        }
        System.out.println("\t-size\t\tWill simulate a chromosome of length "
                + (genomeSizeGenetic / 1000000.0) + " cM and " + (genomeSizePhysical / 1000000.0) + " Mb.");
        if (showDefaultsAndExamples) {
            System.out.println("\t\t\t(Default = 10 cM, 10 Mb; example: \"-size 10\")");
        }
        if (!useMapFile) {
            System.out.println("\t-rec\t\tRecombination rate per base pair is " + recRate + ".");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = 1.0E-8; example: \"-rec 1E-8\")");
            }
            System.out.println("\t-mut\t\tMutation rate per base pair is " + mutationRate + ".");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = 1.65E-8; example: \"-mut 1.65E-8\")");
            }
            System.out.println("\t-map\t\tMap file not specified.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = no map; example: \"-map map.txt\")");
            }
        } else {
            System.out.println("\t-rec\t\tRecombination rate per base pair specified in map file.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = 1.0E-8; example: \"-rec 1E-8\")");
            }
            System.out.println("\t-mut\t\tMutation rate per base pair specified in map file.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = 1.65E-8; example: \"-mut 1.65E-8\")");
            }
            System.out.println("\t-map\t\tUsing map file " + mapFile + ".");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = no map; example: \"-map map.txt\")");
            }
        }
        System.out.println("\t-seed\t\tRandom seed is " + seed + ".");
        if (showDefaultsAndExamples) {
            System.out.println("\t\t\t(Default = random; example: \"-seed 1234\")");
        }
        if (doGeneConversion) {
            System.out.println("\t-GC\t\tNon-crossover gene conversion (NCOGC) is active. Non-crossover/crossover ratio is " + geneConvRelativeRate + ".");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = NCOGC disabled; example: \"-GC 1.0\")");
            }
            System.out.println("\t-meanGC\t\tNCOGC tracts will have a mean length of " + meanGeneConversionTractLength + " bp.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = 300; example: -meanGC 300\")");
            }
            System.out.println("\t-minGC\t\tNCOGC tracts will have a minimum length of " + minimumGeneConversionTractLength + " bp.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = 0; example: -minGC 0\")");
            }
        } else {
            System.out.println("\t-GC\t\tNon-crossover gene conversion (NCOGC) is not active.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = non-crossover GC disabled; example: \"-GC 1.0\")");
            }
            System.out.println("\t-meanGC\t\tMean NCOGC tract length is unused (NCOGC is not active).");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = 300; example: \"-meanGC 300\")");
            }
            System.out.println("\t-minGC\t\tMinimum NCOGC tract length is unused (NCOGC is not active).");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = 0; example: \"-minGC 0\")");
            }
        }

        System.out.println("\nOutput options:");
        if (writeOutFiles) {
            System.out.println("\t-out\t\tWill write output in files " + fileBase + ".*");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = output to files ARGON.*; example: \"-out ARGON\")");
            }
            System.out.println("\t-screen\t\tNot active (will write to file instead).");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = output to files; example: \"-screen\")");
            }
        } else {
            System.out.println("\t-out\t\tOutput file will be ignored because -screen was activated.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = output to files ARGON.*; example: \"-out ARGON\")");
            }
            System.out.println("\t-screen\t\tIs active, will write output on screen.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = output to files; example: \"-screen\")");
            }
        }
        if (printMUT) {
            String format = (!seqOutFormat && !hapsOutFormat) ? "VCF" : (seqOutFormat) ? "seq" : "haps/samples";
            System.out.println("\t-seq\t\tWill print sequence (" + format + " format). Minimum allele frequency " + MAF + ".");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = print sequence, MAF = 0.0; example: \"-seq 0.0\")");
            }
        } else {
            System.out.println("\t-seq\t\tWill not print sequence.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = print sequence, MAF = 0.0; example: \"-seq 0.0\")");
            }
        }
        if (seqOutFormat) {
            System.out.println("\t-seq-out\tWill use seq file format for output.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = use VCF file format; example: \"-seq-out\")");
            }
        } else {
            System.out.println("\t-seq-out\tWill not use seq file format for output.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = use VCF file format; example: \"-seq-out\")");
            }
        }
        if (hapsOutFormat) {
            System.out.println("\t-haps-out\tWill use haps/samples file format for output.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = use VCF file format; example: \"-haps-out\")");
            }
        } else {
            System.out.println("\t-haps-out\tWill not use haps/samples file format for output.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = use VCF file format; example: \"-haps-out\")");
            }
        }
        if (printIBD) {
            System.out.println("\t-IBD\t\tWill output IBD segments longer than " + IBDthresholdCM + " cM.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = do not print IBD; example: \"-IBD 1.0\")");
            }
        } else {
            System.out.println("\t-IBD\t\tWill not output IBD.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = do not print IBD; example: \"-IBD 1.0\")");
            }
        }
//        if (printAFS) {
//            System.out.println("\t-printAFS\twill print the allele frequency spectrum.");
//        } else {
//            System.out.println("\t-printAFS\twill not print the allele frequency spectrum.");
//        }
//        if (writeAFS) {
//            System.out.println("\t-writeAFS\twill save the allele frequency spectrum to file.");
//        } else {
//            System.out.println("\t-writeAFS\twill not save the allele frequency spectrum to file.");
//        }
        if (shrinkSequenceData) {
            System.out.println("\t-shrink\t\tWill output sequence as sample list insteaf of binary list when " + thresholdForPrintingBooleanRatherThanList + " samples carry the derived allele.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = do not shrink sequence; example: \"-shrink\")");
            }
        } else {
            System.out.println("\t-shrink\t\tNot active. Will always output sequence as binary list.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = do not shrink sequence; example: \"-shrink\")");
            }
        }
        if (quiet) {
            System.out.println("\t-quiet\t\tActive. Will not print progress.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = print progress; example: \"-quiet\")");
            }
        } else {
            System.out.println("\t-quiet\t\tNot active. Will print progress.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = print progress; example: \"-quiet\")");
            }
        }
        if (outputAlleleAge) {
            System.out.println("\t-age\t\tWill output age of alleles.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = do not output age of alleles; example: \"-age\")");
            }
        } else {
            System.out.println("\t-age\t\tWill not output age of alleles.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = do not output age of alleles; example: \"-age\")");
            }
        }
        if (gzipOutput) {
            System.out.println("\t-gz\t\tWill compress output.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = do not compress output; example: \"-gz\")");
            }
        } else {
            System.out.println("\t-gz\t\tWill not compress output.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = do not compress output; example: \"-gz\")");
            }
        }
        if (outputTrees) {
            System.out.println("\t-trees\t\tWill output Newick trees.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = do not output trees; example: \"-trees\")");
            }
        } else {
            System.out.println("\t-trees\t\tWill not output Newick trees.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = do not output trees; example: \"-trees\")");
            }
        }
        if (showDefaultsAndExamples) {
            System.out.println("\t-help\t\tPrinting parameter defaults and examples.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = do not print parameter defaults and examples; example: \"-help\")");
            }
        } else {
            System.out.println("\t-help\t\tNot printing parameter defaults and examples.");
            if (showDefaultsAndExamples) {
                System.out.println("\t\t\t(Default = do not print parameter defaults and examples; example: \"-help\")");
            }
        }
        System.out.println();
    }

    public int genToPhys(int genPos) {
        if (!recToPhysicalMap.containsKey(genPos)) {
            Map.Entry<Integer, ArrayList<Integer>> previousEntry = recToPhysicalMap.floorEntry(genPos);
            Map.Entry<Integer, ArrayList<Integer>> nextEntry = recToPhysicalMap.ceilingEntry(genPos);
            double previousRecPos = previousEntry.getKey();
            ArrayList<Integer> previousPhysPos = previousEntry.getValue();
            int pPhys = previousPhysPos.get(0);
            double nextRecPos = nextEntry.getKey();
            ArrayList<Integer> nextPhysPos = nextEntry.getValue();
            int nPhys = nextPhysPos.get(0);
            double frac = (((double) genPos) - previousRecPos) / (nextRecPos - previousRecPos);
            int newPhysPos = pPhys + (int) Math.round(frac * (nPhys - pPhys));
            return newPhysPos;
        } else {
            return recToPhysicalMap.get(genPos).get(0);
        }
    }

    public void addGenPosToMaps(int genPos) {
        if (!recToPhysicalMap.containsKey(genPos)) {
            Map.Entry<Integer, ArrayList<Integer>> previousEntry = recToPhysicalMap.floorEntry(genPos);
            Map.Entry<Integer, ArrayList<Integer>> nextEntry = recToPhysicalMap.ceilingEntry(genPos);
            double previousRecPos = previousEntry.getKey();
            ArrayList<Integer> previousPhysPos = previousEntry.getValue();
            int pPhys = previousPhysPos.get(0);
            float previousRecRate = recRateMap.get(pPhys);
            float previousMutRate = mutRateMap.get(pPhys);
            double nextRecPos = nextEntry.getKey();
            ArrayList<Integer> nextPhysPos = nextEntry.getValue();
            int nPhys = nextPhysPos.get(0);
            double frac = (((double) genPos) - previousRecPos) / (nextRecPos - previousRecPos);
            int newPhysPos = pPhys + (int) Math.round(frac * (nPhys - pPhys));
            ArrayList<Integer> p = new ArrayList<Integer>();
            p.add(newPhysPos);
            recToPhysicalMap.put(genPos, p);
            recRateMap.put(newPhysPos, previousRecRate);
            mutRateMap.put(newPhysPos, previousMutRate);
        }
    }

    public static void writeNewickTree(StringBuilder out, Block node, int distFromParent, int fromPos, int toPos) {
        if (node.children.isEmpty()) {
            out.append(node.gen).append("_").append(node.ID).append(":").append(distFromParent);
        } else {
            out.append("(");
            for (int i = 0; i < node.children.size() - 1; i++) {
                int distToChild = node.gen - node.children.get(i).gen;
                writeNewickTree(out, node.children.get(i), distToChild, fromPos, toPos);
                out.append(",");
            }
            int distToChild = node.gen - node.children.get(node.children.size() - 1).gen;
            writeNewickTree(out, node.children.get(node.children.size() - 1), distToChild, fromPos, toPos);
            out.append(")");
            if (distFromParent > 0) {
                out.append(node.gen).append("_").append(node.ID).append(":").append(distFromParent);
            } else {
                out.append(node.gen).append("_").append(node.ID);
            }
        }
    }

    public void outputMut(int fromPhys, int toPhys, float mut, int myGen, int parentGen, int distanceToFather, int size, float freq, DescendantsList myList) throws IOException {
        int physLen = toPhys - fromPhys + 1;
        double lambda = mut * distanceToFather * physLen;
        int mutations = Tools.samplePoisson(lambda, generator);
        if (writeAFS || printAFS) {
            AFS[size]++;
            totMut += mutations;
        }
        if (mutations > 0 && freq >= MAF) {
            StringBuffer list;
            if (seqOutFormat) { // this format is obsolete, should remove it
                if (shrinkSequenceData && myList.getSize() < thresholdForPrintingBooleanRatherThanList) {
                    list = myList.getSetList();
                } else {
                    list = myList.getBinaryList();
                }
                String str = fromPhys + "\t" + toPhys + "\t" + mutations + "\t" + freq + "\t" + list + "\n";
                if (writeOutFiles && printMUT) {
                    try {
                        MUTwriter.write(str);
                    } catch (IOException ex) {
                        Logger.getLogger(Argon.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    if (outputAlleleAge) {
                        int mutAge = myGen + generator.nextInt(distanceToFather) + 1;
                        AGEwriter.write(chromosome + "\t" + fromPhys + "\t" + toPhys + "\t" + mutAge + "\t" + myGen + "\t" + parentGen + "\n");
                    }
                } else if (printMUT) {
                    System.out.print("MUT\t" + str);
                    if (outputAlleleAge) {
                        int mutAge = myGen + generator.nextInt(distanceToFather) + 1;
                        System.out.println("AGE\t" + chromosome + "\t" + fromPhys + "\t" + toPhys + "\t" + mutAge + "\t" + myGen + "\t" + parentGen);
                    }
                }
            } else if (hapsOutFormat) {
                for (int k = 0; k < mutations; k++) {
                    numMutations++;
                    list = myList.getBinaryListWithSpaces();
                    int pos = fromPhys + generator.nextInt(toPhys - fromPhys + 1);
                    String SNPname = "SNP_" + pos + "_" + numMutations;
                    String str = chromosome + ":" + pos + "_1_2 " + SNPname + " " + pos + " 1 2 " + list + "\n";
                    if (writeOutFiles && printMUT) {
                        try {
                            MUTwriter.write(str);
                        } catch (IOException ex) {
                            Logger.getLogger(Argon.class.getName()).log(Level.SEVERE, null, ex);
                        }
                        if (outputAlleleAge) {
                            int mutAge = myGen + generator.nextInt(distanceToFather) + 1;
                            AGEwriter.write(chromosome + "\t" + SNPname + "\t" + pos + "\t" + fromPhys + "\t" + toPhys + "\t" + mutAge + "\t" + myGen + "\t" + parentGen + "\n");
                        }
                    } else if (printMUT) {
                        if (outputAlleleAge) {
                            System.out.print("MUTa\t" + str);
                            int mutAge = myGen + generator.nextInt(distanceToFather) + 1;
                            System.out.print("AGE\t" + chromosome + "\t" + SNPname + "\t" + pos + "\t" + fromPhys + "\t" + toPhys + "\t" + mutAge + "\t" + myGen + "\t" + parentGen);
                        }
                    }
                }
            } else {
                for (int k = 0; k < mutations; k++) {
                    numMutations++;
                    list = myList.getBinaryList();
                    int pos = fromPhys + generator.nextInt(toPhys - fromPhys + 1);
                    StringBuilder outString = new StringBuilder();
                    outString.append(chromosome).append("\t").append(pos).append("\tSNP_").append(numMutations).append("\t1\t2\t.\t.\tPR\tGT");
                    for (int i = 0; i < list.length(); i += 2) {
                        outString.append("\t").append(list.charAt(i)).append("|").append(list.charAt(i + 1));
                    }
                    outString.append("\n");
                    if (writeOutFiles && printMUT) {
                        try {
                            MUTwriter.write(outString.toString());
                        } catch (IOException ex) {
                            Logger.getLogger(Argon.class.getName()).log(Level.SEVERE, null, ex);
                        }
                        if (outputAlleleAge) {
                            String SNPname = "SNP_" + pos + "_" + numMutations;
                            int mutAge = myGen + generator.nextInt(distanceToFather) + 1;
                            AGEwriter.write(chromosome + "\t" + SNPname + "\t" + pos + "\t" + fromPhys + "\t" + toPhys + "\t" + mutAge + "\t" + myGen + "\t" + parentGen + "\n");
                        }
                    } else if (printMUT) {
                        System.out.print("MUT\t" + outString.toString());
                        if (outputAlleleAge) {
                            String SNPname = "SNP_" + pos + "_" + numMutations;
                            int mutAge = myGen + generator.nextInt(distanceToFather) + 1;
                            System.out.println("AGE\t" + chromosome + "\t" + SNPname + "\t" + pos + "\t" + fromPhys + "\t" + toPhys + "\t" + mutAge + "\t" + myGen + "\t" + parentGen);
                        }
                    }
                }
            }
        }
    }

    public void printVCFheader() {
        SimpleDateFormat sdate = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        Date date = new Date();
        String datestr = "##fileDate=[" + sdate.format(date) + "]";
        String versionstr = "##source=ARGONv" + version;
        if (writeOutFiles && printMUT) {
            try {
                MUTwriter.write("##fileformat=VCFv4.2\n");
                MUTwriter.write(datestr + "\n");
                MUTwriter.write(versionstr + "\n");
                MUTwriter.write("##contig=<ID=" + chromosome + ",length=" + genomeSizePhysical + ">\n");
                MUTwriter.write("##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Simulated data, reference alleles are ancestral\"\n");
                MUTwriter.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
                MUTwriter.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
                for (int i = 0; i < sampleSizes.length; i++) {
                    if (sampleSizes[i] % 2 != 0) {
                        exit("If the VCF file format is used, the number of samples should be even for all populations.");
                    }
                    int cnt = 0;
                    int popInd = i + 1;
                    for (int j = 0; j < sampleSizes[i]; j += 2) {
                        cnt++;
                        MUTwriter.write("\t" + popInd + "_" + cnt);
                    }
                }
                MUTwriter.write("\n");

            } catch (IOException ex) {
                Logger.getLogger(Argon.class
                        .getName()).log(Level.SEVERE, null, ex);
            }
        } else if (printMUT) {
            System.out.print("MUT ##fileformat=VCFv4.2\n");
            System.out.print("MUT " + datestr + "\n");
            System.out.print("MUT " + versionstr + "\n");
            System.out.print("MUT ##contig=<ID=" + chromosome + ",length=" + genomeSizePhysical + ">\n");
            System.out.print("MUT ##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Simulated data, reference alleles are ancestral\"\n");
            System.out.print("MUT ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            System.out.print("MUT #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
            for (int i = 0; i < sampleSizes.length; i++) {
                if (sampleSizes[i] % 2 != 0) {
                    exit("If the VCF file format is used, the number of samples should be even for all populations.");
                }
                int cnt = 0;
                int popInd = i + 1;
                for (int j = 0; j < sampleSizes[i]; j += 2) {
                    cnt++;
                    System.out.print("\t" + popInd + "_" + cnt);
                }
            }
            System.out.print("\n");
        }
    }

    public void printSamplesFile() {
        if (writeOutFiles && printMUT) {
            String name = fileBase + ".samples";
            FileOutputStream samplesOutput = null;
            OutputStreamWriter samplesWriter = null;
            try {
                samplesOutput = new FileOutputStream(name);
            } catch (FileNotFoundException e) {
                exit("file " + name + " could not be created.");
            }
            try {
                samplesWriter = new OutputStreamWriter(samplesOutput);
            } catch (Exception e) {
                exit("could not open output stream for " + name);
            }
            try {
                samplesWriter.write("ID_1 ID_2 missing\n");
                samplesWriter.write("0 0 0\n");
                for (int i = 0; i < sampleSizes.length; i++) {
                    int cnt = 0;
                    int popInd = i + 1;
                    for (int j = 0; j < sampleSizes[i]; j += 2) {
                        cnt++;
                        String ID1 = popInd + "_" + cnt;
                        String ID2 = popInd + "_" + cnt;
                        samplesWriter.write(ID1 + " " + ID2 + " 0\n");
                    }
                }
            } catch (IOException ex) {
                Logger.getLogger(Argon.class.getName()).log(Level.SEVERE, null, ex);
            }
            try {
                samplesWriter.flush();
                samplesWriter.close();
            } catch (IOException ex) {
                Logger.getLogger(Argon.class.getName()).log(Level.SEVERE, null, ex);
            }
        } else if (printMUT) {
            System.out.print("MUT ID_1 ID_2 missing\n");
            System.out.print("MUT 0 0 0\n");
            for (int i = 0; i < sampleSizes.length; i++) {
                int cnt = 0;
                int popInd = i + 1;
                for (int j = 0; j < sampleSizes[i]; j += 2) {
                    cnt++;
                    String ID1 = popInd + "_" + cnt;
                    String ID2 = popInd + "_" + cnt;
                    System.out.print("MUT " + ID1 + " " + ID2 + " 0\n");
                }
            }
        }
    }

    public void exit(String error) {
        printParams(true);
        System.err.println("ERROR:\t" + error + "\n");
        System.exit(1);
    }

}
