/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ARGON;

/**
 *
 * @author Pier Palamara <ppalama@hsph.harvard.edu>
 */
public class ArgonMain {

    public static void main(String[] args) throws Exception {

        boolean quiet = false;

        boolean shrinkSequenceData = false;
        boolean writeOutFiles = true;
        boolean printIBD = false;
        boolean printMUT = true;

        boolean printAFS = false;
        boolean writeAFS = false;
        boolean usePopFile = false;
        boolean useMapFile = false;
        boolean specifiedMutRate = false;
        boolean specifiedRecRate = false;
        boolean doGeneConversion = false;
        boolean gzipOutput = false;
        boolean popFlagUsed = false;
        boolean outputAlleleAge = false;
        boolean seqOutFormat = false;
        boolean hapsOutFormat = false;
        boolean outputTrees = false;

        float recRate = 0.00000001f;
        float MAF = 0.0f;
        float mutationRate = 0.0000000165f;
        float recRateFactor = recRate * 100000000.0f;
        float IBDthresholdCM = 1.f;
        float geneConvRelativeRate = 1.f;
        float probGeneConversionGivenRec = geneConvRelativeRate / (1.f + geneConvRelativeRate);

        String popSizeFile = "";
        String mapFile = "";
        String fileBase = "ARGON";

        int numPops = 1;
        Long constantPopSize = 1000l;
        int sampleSizes[] = null;
        int genomeSizePhysical = 10000000;
        int minBlockSize = 1;
        int totModernSamples = 0;
        int meanGeneConversionTractLength = 300;
        int minimumGeneConversionTractLength = 0;

        long seed = System.currentTimeMillis();

        boolean showDefaultsAndExamples = false;

        System.out.println();
        System.out.println(Argon.versionString + "\n");

        //parse arguments
        int argIndex = 0;
        try {
            while (argIndex < args.length) {
                String arg = args[argIndex++];
//            System.out.println("arg " + argIndex + "\t" + arg);
                if (arg.equals("-N")) {
                    try {
                        constantPopSize = Long.parseLong(args[argIndex]);
//                        System.out.println("Parsed constant population size " + constantPopSize);
                    } catch (NumberFormatException e) {
                        usePopFile = true;
                        popSizeFile = args[argIndex];
//                        System.out.println("Parsed population size file " + popSizeFile);
                    }
                    argIndex++;
                } else if (arg.equals("-pop")) {
                    try {
                        numPops = Integer.parseInt(args[argIndex++]);
//                    System.out.println("There will be " + numPops + " population(s).");
                        if (popFlagUsed) {
                            System.out.println(popFlagUsed);
                            exit("-pop parameter specified twice.");
                            System.out.println(popFlagUsed);
                        }
                        sampleSizes = new int[numPops];
                        for (int i = 0; i < numPops; i++) {
                            sampleSizes[i] = Integer.parseInt(args[argIndex++]);
                            totModernSamples += sampleSizes[i];
//                        System.out.println("\tPop " + (i + 1) + " will have " + sampleSizes[i] + " samples.");
                        }
                        popFlagUsed = true;
                    } catch (Exception e) {
                        exit("Error parsing -pop parameters.");
                    }
                } else if (arg.equals("-age")) {
                    outputAlleleAge = true;
                } else if (arg.equals("-seq-out")) {
                    seqOutFormat = true;
                } else if (arg.equals("-trees")) {
                    outputTrees = true;
                } else if (arg.equals("-haps-out")) {
                    hapsOutFormat = true;
                } else if (arg.equals("-size")) {
                    try {
                        double val = Double.parseDouble(args[argIndex++]);
                        if (val > 2146.) {
                            exit("Region size cannot exceed " + 2146.);
                        }
                        genomeSizePhysical = (int) (val * 1000000);
                    } catch (NumberFormatException e) {
                        exit("Error parsing -size parameter.");
                    }
                } else if (arg.equals("-len")) {
                    try {
                        minBlockSize = Integer.parseInt(args[argIndex++]);
                    } catch (NumberFormatException e) {
                        exit("Error parsing -len parameter.");
                    }
                    int excessBlocks = genomeSizePhysical % minBlockSize;
                    if (excessBlocks > 0) {
                        genomeSizePhysical = Tools.roundToBlock(genomeSizePhysical, minBlockSize);
                        System.out.println("Warning: adjusting sequence length to " + genomeSizePhysical + " to match minimum block size requirements.\n");
                    }
                } else if (arg.equals("-rec")) {
                    try {
                        recRate = Float.parseFloat(args[argIndex]);
                        recRateFactor = recRate * 100000000.0f;
                        specifiedRecRate = true;
                    } catch (NumberFormatException e) {
                        exit("Error parsing -rec parameter.");
                    }
                    argIndex++;
                } else if (arg.equals("-seed")) {
                    try {
                        seed = Long.parseLong(args[argIndex++]);
                    } catch (NumberFormatException e) {
                        exit("Error parsing -seed parameter.");
                    }
                } else if (arg.equals("-mut")) {
                    try {
                        mutationRate = Float.parseFloat(args[argIndex]);
                        specifiedMutRate = true;
                    } catch (NumberFormatException e) {
                        exit("Error parsing -mut parameter.");
                    }
                    argIndex++;
                } else if (arg.equals("-map")) {
                    useMapFile = true;
                    mapFile = args[argIndex];
                    argIndex++;
                } else if (arg.equals("-seq")) {
                    printMUT = true;
                    try {
                        MAF = Float.parseFloat(args[argIndex++]);
                        if (MAF < 0 || MAF > 1.0) {
                            exit("MAF value must to be between 0 and 1.");
                        }
                    } catch (NumberFormatException e) {
                        exit("Error parsing -maf parameter.");
                    }
                } else if (arg.equals("-out")) {
                    writeOutFiles = true;
                    fileBase = args[argIndex++];
                } else if (arg.equals("-screen")) {
                    writeOutFiles = false;
                } else if (arg.equals("-help")) {
                    showDefaultsAndExamples = true;
                } else if (arg.equals("-IBD")) {
                    printIBD = true;
                    IBDthresholdCM = Float.parseFloat(args[argIndex++]);
                } else if (arg.equals("-GC")) {
                    doGeneConversion = true;
                    geneConvRelativeRate = Float.parseFloat(args[argIndex++]);
                    probGeneConversionGivenRec = geneConvRelativeRate / (1.f + geneConvRelativeRate);
                } else if (arg.equals("-meanGC")) {
                    meanGeneConversionTractLength = Integer.parseInt(args[argIndex++]);
                    if (meanGeneConversionTractLength < 0) {
                        exit("mean gene conversion length cannot be less than 0.");
                    }
                } else if (arg.equals("-minGC")) {
                    minimumGeneConversionTractLength = Integer.parseInt(args[argIndex++]);;
                    if (minimumGeneConversionTractLength < 0) {
                        minimumGeneConversionTractLength = 0;
                    }
                } else if (arg.equals("-printAFS")) {
                    printAFS = true;
                } else if (arg.equals("-gz")) {
                    gzipOutput = true;
                } else if (arg.equals("-quiet")) {
                    quiet = true;
                } else if (arg.equals("-saveAFS")) {
                    writeAFS = true;
                } else if (arg.equals("-shrink")) {
                    shrinkSequenceData = true;
                } else {
                    exit("Error reading unknown input argument " + arg + ".");
                }
            }
        } catch (Exception e) {
            exit("Error parsing input arguments.");
        }

        if ((specifiedMutRate || specifiedRecRate) && useMapFile) {
            exit("if a map is specified, the -rec and -mut flags cannot be used");
        }

        if (doGeneConversion && minBlockSize > 1) {
            exit("cannot use block size approximation (-len > 1) if gene conversion is simulated.");
        }

        if (sampleSizes == null) {
            sampleSizes = new int[]{1000};
            totModernSamples = sampleSizes[0];
            numPops = 1;
        }

        Argon argon = new Argon(
                quiet,
                shrinkSequenceData,
                writeOutFiles,
                printIBD,
                printMUT,
                printAFS,
                writeAFS,
                usePopFile,
                useMapFile,
                specifiedMutRate,
                specifiedRecRate,
                doGeneConversion,
                gzipOutput,
                showDefaultsAndExamples,
                outputAlleleAge,
                seqOutFormat,
                hapsOutFormat,
                outputTrees,
                recRate,
                MAF,
                mutationRate,
                recRateFactor,
                IBDthresholdCM,
                geneConvRelativeRate,
                probGeneConversionGivenRec,
                popSizeFile,
                mapFile,
                fileBase,
                numPops,
                sampleSizes,
                genomeSizePhysical,
                minBlockSize,
                totModernSamples,
                meanGeneConversionTractLength,
                minimumGeneConversionTractLength,
                constantPopSize,
                seed
        );

        argon.simulate();
        argon.visit();

        System.out.println("Total runtime: " + (argon.ARG_runtime + argon.Visit_runtime));
        System.out.println();
    }

    public static void exit(String error) {
        printHelp();
        System.err.println("ERROR:\t" + error + "\n");
        System.exit(1);
    }

    public static void printHelp() {
        
        
        // TODO: avoid hard-coding default values. Should handle args using ARGParse library or similar.
        System.out.println("\nParameter values:");

        System.out.println("\t-N\t\tPopulation size.");
        System.out.println("\t\t\t(Default = 1000; example: \"-N 10000\")");

        System.out.println("\t-pop\t\tNumber of samples.");
        System.out.println("\t\t\t(Default = 1 pop, 1000 samples; example: \"-pop 1 1000\" or \"-pop 2 1000 2000\")");

        System.out.println("\t-len\t\tMinimum block size (MicroMorgans).");
        System.out.println("\t\t\t(Default = 1 mM; example: \"-len 1\" or \"-len 10000\")");

        System.out.println("\t-size\t\tChromosome length (Mb).");
        System.out.println("\t\t\t(Default = 10 cM, 10 Mb; example: \"-size 10\")");

        System.out.println("\t-rec\t\tRecombination rate per base pair.");
        System.out.println("\t\t\t(Default = 1.0E-8; example: \"-rec 1E-8\")");

        System.out.println("\t-mut\t\tMutation rate per base pair.");
        System.out.println("\t\t\t(Default = 1.65E-8; example: \"-mut 1.65E-8\")");

        System.out.println("\t-map\t\tMap file.");
        System.out.println("\t\t\t(Default = no map; example: \"-map map.txt\")");

        System.out.println("\t-seed\t\tRandom seed.");
        System.out.println("\t\t\t(Default = system time; example: \"-seed 1234\")");

        System.out.println("\t-GC\t\tActivate non-crossover gene conversion (NCOGC) (specify non-crossover/crossover ratio).");
        System.out.println("\t\t\t(Default = NCOGC disabled; example: \"-GC 1.0\")");

        System.out.println("\t-meanGC\t\tMean length of NCOGC tracts (specify length in bp).");
        System.out.println("\t\t\t(Default = 300; example: -meanGC 300\")");

        System.out.println("\t-minGC\t\tMinimum length of NCOGC tracts (specify length in bp).");
        System.out.println("\t\t\t(Default = 0; example: -minGC 0\")");


        
        System.out.println("\nOutput options:");

        System.out.println("\t-out\t\tRoot for file output.");
        System.out.println("\t\t\t(Default = output to files ARGON.*; example: \"-out ARGON\")");

        System.out.println("\t-screen\t\tOutput to screen.");
        System.out.println("\t\t\t(Default = not active, output to files; example: \"-screen\")");

        System.out.println("\t-seq\t\tPrint sequence (specify minimum allele frequency).");
        System.out.println("\t\t\t(Default = print sequence, MAF = 0.0; example: \"-seq 0.0\")");

        System.out.println("\t-seq-out\tUse seq file format for output.");
        System.out.println("\t\t\t(Default = do not seq file format for output (use VCF); example: \"-seq-out\")");

        System.out.println("\t-haps-out\tUse haps/samples file format for output.");
        System.out.println("\t\t\t(Default = do not haps/samples file format for output (use VCF); example: \"-haps-out\")");

        System.out.println("\t-IBD\t\tOutput IBD segments (specify minimum length threshold in cM).");
        System.out.println("\t\t\t(Default = do not print IBD; example: \"-IBD 1.0\")");

        System.out.println("\t-shrink\t\tActivate to output sequence as sample list insteaf of binary list when this requires less output (only when -seq format is used).");
        System.out.println("\t\t\t(Default = do not shrink sequence; example: \"-shrink\")");

        System.out.println("\t-quiet\t\tDo not print progress.");
        System.out.println("\t\t\t(Default = print progress; example: \"-quiet\")");

        System.out.println("\t-age\t\tOutput age of alleles.");
        System.out.println("\t\t\t(Default = do not output age of alleles; example: \"-age\")");

        System.out.println("\t-gz\t\tCompress output.");
        System.out.println("\t\t\t(Default = do not compress output; example: \"-gz\")");

        System.out.println("\t-trees\t\tOutput Newick trees to *.trees(.gz) file.");
        System.out.println("\t\t\t(Default = do not output trees; example: \"-trees\")");

        System.out.println("\t-help\t\tPrint parameter defaults and examples.");
        System.out.println("\t\t\t(Default = do not print parameter defaults and examples; example: \"-help\")");

        System.out.println();
    }

}
