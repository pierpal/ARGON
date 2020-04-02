package ARGON;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.TreeSet;

/**
 *
 * @author Pier Palamara <ppalama@hsph.harvard.edu>
 */
public class Individual {

    static final boolean verboseCoalesce = false;
    static final boolean verboseRecomb = false;
    static final boolean stopTrackingBlockVerbose = false;

    static int tempCounterBlocksDone = 0;

    ArrayList<Block> segments = new ArrayList<Block>();
    long ID;

    Individual(Individual id) {
        this.ID = id.ID;
        this.segments = id.segments;
    }

    Individual(long ID, ArrayList<Block> segments) {
        this.ID = ID;
        this.segments = segments;
    }

    Individual(int ID) {
        this.ID = ID;
    }

    Individual(long ID, int numBlocks) {
        this.ID = ID;
        Block b = new Block(ID, 0, 1, numBlocks, 1);
        segments.add(b);
    }

    void print() {
        System.out.print(ID);
        Block b = null;
        for (int i = 0; i < segments.size(); i++) {
            b = segments.get(i);
            System.out.print(" " + b.start + "-" + b.end);
        }
        System.out.println();
    }

    static void print(ArrayList<Block> segments) {
        Block b = null;
        for (int i = 0; i < segments.size(); i++) {
            b = segments.get(i);
            System.out.print(" " + b.getRangeString());
        }
        System.out.println();
    }

    public ArrayList<Block> splitAfterSomeBlocks(int targetBlock, long ID, int gen, boolean inGC) {
//        if (verboseRecomb) {
//            System.out.println("************ START RECOMB at gen " + gen + " targetPos: " + targetBlock);
//            System.out.println("GC is " + inGC);
//            System.out.print("start this ind ");
//            print(segments);
//        }
        ArrayList<Block> newMe = new ArrayList<Block>();
        ArrayList<Block> otherInd = new ArrayList<Block>();
        int pos = 0;
        for (; pos < segments.size(); pos++) {
            if (segments.get(pos).start >= targetBlock) {
                break; // from here on it's next ind
            } else if (segments.get(pos).start < targetBlock && segments.get(pos).end >= targetBlock) {
                Block recombinant = segments.get(pos);
                // When we find a block containing the target, we split it up by
                //    start  - target-1
                //    target - end
                // If inGC is true, we are terminating a GC block, so the lefmost fragment is labeled GC
                Block chunkInNewThisInd = new Block(ID, gen, recombinant.start, targetBlock - 1, recombinant.numberOfDescendants, inGC);
                Block chunkInOtherInd = new Block(ID, gen, targetBlock, recombinant.end, recombinant.numberOfDescendants);
                newMe.add(chunkInNewThisInd);
                otherInd.add(chunkInOtherInd);
                chunkInNewThisInd.representedChildren = recombinant.representedChildren;
                chunkInOtherInd.representedChildren = recombinant.representedChildren;
//                if (verboseRecomb) {
//                    System.out.println(chunkInNewThisInd.getRangeString() + " now represents " + chunkInNewThisInd.representedChildren.getRangeString());
//                    System.out.println(chunkInOtherInd.getRangeString() + " now represents " + chunkInOtherInd.representedChildren.getRangeString());
//                }
                pos++;
                break;
            } else {
                newMe.add(segments.get(pos));
            }
        }
        for (; pos < segments.size(); pos++) {
            otherInd.add(segments.get(pos));
        }
        segments = newMe;
//        if (verboseRecomb) {
//            System.out.print("end this ind ");
//            print(segments);
//            System.out.print("end other ind ");
//            print(otherInd);
//            System.out.println("************ END RECOMB");
//        }
        return otherInd;
    }

    public ArrayList<Block> mergeArrayLists(ArrayList<Block> otherBlocks, int gen, long ID, int samples, ArrayList<Block> roots, HashMap<Block, HashSet<Block>> childrenOfThisGeneration) {
//        if (verboseCoalesce) {
//            System.out.println("************* START COALESCENCE ************");
//            System.out.print("segments ");
//            print(segments);
//            System.out.print("other ");
//            print(otherBlocks);
//        }
        ArrayList<Block> destination = new ArrayList<Block>(), mySegments = segments;
        ArrayList<Block> overlaps = new ArrayList<Block>();
        ArrayList<Block> temp;
        TreeSet<Block> created = new TreeSet<Block>(); // will contain all new blocks created in destination, which will be parents of all fragmented blocks in sources
        PriorityQueue<Block> fragmentedInSources = new PriorityQueue<Block>(); // <-- all blocks in sources that happen to be fragmented in the process, so that they will be linked to a parent in the created set. Using these structures because a block could be fragmented multiple times, and can do the parent linking at the end. FIXME: Could be done more efficiently and broken blocks could be merged sometimes.
        int count = 0;
        while (segments.size() > 0 || otherBlocks.size() > 0) {
            count++;
//            if (verboseCoalesce) {
//                System.out.print("it " + count + " segments ");
//                print(segments);
//                System.out.print("it " + count + " other ");
//                print(otherBlocks);
//                System.out.print("it " + count + " destination ");
//                print(destination);
//            }
            if (otherBlocks.size() == 0) {
                temp = mySegments;
            } else if (mySegments.size() == 0) {
                temp = otherBlocks;
            } else if (mySegments.get(0).start < otherBlocks.get(0).start) {
                temp = mySegments;
            } else {
                temp = otherBlocks;
            }
            if (destination.size() == 0) { // first iteration
                destination.add(temp.get(0));
                temp.remove(0);
            } else { // not first
                Block testedBlockInSource = temp.get(0);
                Block testedBlockInDestination = destination.get(destination.size() - 1);
                int descendantsSum = testedBlockInDestination.numberOfDescendants + testedBlockInSource.numberOfDescendants;
                Block coalesceRange = testedBlockInSource.overlap(testedBlockInDestination);
                if (coalesceRange != null) {
                    coalesceRange.numberOfDescendants = descendantsSum;
                    coalesceRange.ID = ID;
                    coalesceRange.gen = gen;
//                    if (verboseCoalesce) {
//                        System.out.println("overlap! " + coalesceRange.getRangeString() + " " + testedBlockInDestination.getRangeString());
//                    }
                    fragmentedInSources.add(testedBlockInSource);
//                    if (verboseCoalesce) {
//                        System.out.println("F_added " + testedBlockInSource.getRangeString());
//                    }
                    int end = Math.max(testedBlockInSource.end, testedBlockInDestination.end);
                    if (created.contains(testedBlockInDestination)) {
                        created.remove(testedBlockInDestination); // in the list of created blocks, the last is the one we're fragmenting
//                        if (verboseCoalesce) {
//                            System.out.println("C_removed " + testedBlockInDestination.getRangeString());
//                        }
                    } else {
                        fragmentedInSources.add(testedBlockInDestination);
//                        if (verboseCoalesce) {
//                            System.out.println("F_added " + testedBlockInDestination.getRangeString());
//                        }
                    }
                    coalesceRange.isCoalescent = true; // this is a coalescent block and will show up in marginal genealogies
                    created.add(coalesceRange); // add this block to the list of created blocks
//                    if (verboseCoalesce) {
//                        System.out.println("C1_added " + coalesceRange.getRangeString());
//                    }
                    destination.remove(destination.size() - 1);
                    if (testedBlockInDestination.start < coalesceRange.start) {
                        Block headBlock = new Block(ID, gen, testedBlockInDestination.start, coalesceRange.start - 1, testedBlockInDestination.numberOfDescendants);//, testedBlockInDestination.isGC);
                        headBlock.representedChildren = testedBlockInDestination.representedChildren;
                        created.add(headBlock); // add this block to the list of created blocks
//                        if (verboseCoalesce) {
//                            System.out.println("C2_added " + headBlock.getRangeString() + " " + gen);
//                        }
                        destination.add(headBlock);
                    }
                    coalesceRange.numberOfDescendants = testedBlockInDestination.numberOfDescendants + testedBlockInSource.numberOfDescendants;
                    if (coalesceRange.numberOfDescendants == samples) {
                        roots.add(coalesceRange);
                    } else {
                        destination.add(coalesceRange);
                    }
                    if (coalesceRange.end != end) {
                        int desc;
                        Block representing = null;
                        if (coalesceRange.end < testedBlockInDestination.end) {
                            desc = testedBlockInDestination.numberOfDescendants;
                            representing = testedBlockInDestination.representedChildren;
                        } else {
                            desc = testedBlockInSource.numberOfDescendants;
                            representing = testedBlockInSource.representedChildren;
                        }
                        Block tailBlock = new Block(ID, gen, Math.min(coalesceRange.end, testedBlockInDestination.end) + 1, end, desc);
                        tailBlock.representedChildren = representing;
                        created.add(tailBlock); // add this block to the list of created blocks
//                        if (verboseCoalesce) {
//                            System.out.println("C3_added " + tailBlock.getRangeString());
//                        }
                        created.add(tailBlock); // add this block to the list of created blocks
                        destination.add(tailBlock); //TODO descendants
                    }
                    temp.remove(0);
                    overlaps.add(coalesceRange);
                } else {
                    destination.add(temp.get(0));
                    temp.remove(0);
                }
            }
//            if (verboseCoalesce) {
//                System.out.println("size " + segments.size() + " " + otherBlocks.size());
//            }
        }
//        if (verboseCoalesce) {
//            System.out.print("segments ");
//            print(segments);
//            System.out.print("other ");
//            print(otherBlocks);
//            System.out.print("destination ");
//            print(destination);
//            System.out.print("overlap ");
//            print(overlaps);
//            System.out.println("************* END COALESCENCE  ************");
//        }

        int lastPosCheckedInCreated = 0;
        ArrayList<Block> orderedCreated = new ArrayList<Block>(created);
//        if (verboseCoalesce) {
//            for (Block cr : orderedCreated) {
//                System.out.println("created " + cr.getRangeString());
//            }
//            System.out.println("");
//        }
        while (fragmentedInSources.size() > 0) {
            Block fragment = fragmentedInSources.poll();
            HashSet<Block> childrenThatIshouldCorrect = childrenOfThisGeneration.get(fragment);
//            if (verboseCoalesce) {
//                System.out.println("fragment broken in source: " + fragment.getRangeString());
//            }
            while (fragment.start > orderedCreated.get(lastPosCheckedInCreated).start) {
                lastPosCheckedInCreated++;
            }
            int currentPos = lastPosCheckedInCreated;
            while (currentPos < orderedCreated.size() && fragment.overlap(orderedCreated.get(currentPos)) != null) {
                Block currBlock = orderedCreated.get(currentPos);
//                if (verboseCoalesce) {
//                    System.out.println("\toverlaps with created block: " + currBlock.getRangeString());
//                }
                if (fragment.gen == gen) {
                    if (fragment.isCoalescent) {
                        // if this is a piece of a coalescent block at the current generation, it must be a coalescent block too.
                        currBlock.isCoalescent = true;
                        currBlock.representedChildren = currBlock;
                    }
//                    if (verboseCoalesce) {
//                        System.out.println("\tand comes from this generation. I should link this to recent fragments");
//                    }
                    HashSet<Block> childrenOfThisBlock = childrenOfThisGeneration.get(currBlock);
                    if (childrenOfThisBlock == null) {
                        childrenOfThisBlock = new HashSet<Block>();
                    }
                    for (Block currentChild : childrenThatIshouldCorrect) {
                        childrenOfThisBlock.add(currentChild);
                        currentChild.removeParent(fragment);
                        boolean blockDone = currentChild.addParent(currBlock);
//                        }
//                        if (verboseCoalesce) {
//                            System.out.println("\tfor " + currentChild.getRangeString() + " removed " + fragment.getRangeString() + " added " + currBlock.getRangeString());
//                        }
                    }
                    childrenOfThisGeneration.put(currBlock, childrenOfThisBlock);
                } else {
//                    if (verboseCoalesce) {
//                        System.out.println("\t comes from a lower generation.");
//                    }
                    if (currBlock.isCoalescent) {
//                        if (verboseCoalesce) {
//                            System.out.println("\t... it is a standard coalescence.");
//                        }
                        boolean blockDone = fragment.representedChildren.addParent(currBlock);
//                        if (verboseCoalesce) {
//                            System.out.println("\t\tadded " + currBlock.getRangeString() + " to parents of " + fragment.representedChildren.getRangeString());
//                        }
                        HashSet<Block> childrenOfThisBlock = childrenOfThisGeneration.get(currBlock);
                        if (childrenOfThisBlock == null) {
                            childrenOfThisBlock = new HashSet<Block>();
                        }
                        childrenOfThisBlock.add(fragment.representedChildren);
//                        if (verboseCoalesce) {
//                            System.out.println("\t\tadded " + fragment.representedChildren.getRangeString() + " to children of " + currBlock.getRangeString());
//                        }
                        childrenOfThisGeneration.put(currBlock, childrenOfThisBlock);
                    } else {
//                        if (verboseCoalesce) {
//                            System.out.println("\t\tThis is a leftover. Will update representing block and gen.");
//                        }
                        currBlock.representedChildren = fragment.representedChildren;
                        currBlock.gen = fragment.representedChildren.gen;
                    }
                }
                currentPos++;
            }
        }

        segments = destination;
        return overlaps;

    }
}
