package ARGON;

import java.util.ArrayList;
import java.util.Comparator;

/**
 *
 * @author Pier Palamara <ppalama@hsph.harvard.edu>
 */

public class Block implements Comparable {

    public int start;
    public int end;
    int numberOfDescendants;
    public boolean isCoalescent = false;
    public Block representedChildren = this;
    ArrayList<Block> parents = new ArrayList<Block>();
    ArrayList<Block> children = new ArrayList<Block>();
    public long ID;
    public int gen;
    public int haveFoundAParent = 0;
//    public int ID;
//    public static int count = 0; //remove
//    public static int killed = 0; //remove
    DescendantsList dList;

//    @Override
//    protected void finalize() { //remove
//        killed++; //remove
//    } //remove
    
    Block(long ID, int gen, int start, int end, int descendants) {
//        count++; //remove
//        this.ID = count; //remove
        this.ID = ID;
        this.gen = gen;
        this.start = start;
        this.end = end;
        this.numberOfDescendants = descendants;
    }

    public boolean addParent(Block parent) {
        parents.add(parent);
        parent.children.add(this);
        haveFoundAParent += parent.end - parent.start + 1;
        return haveFoundAParent == (this.end - this.start + 1);
    }

    public void removeParent(Block parent) {
        parents.remove(parent);
        parent.children.remove(this);
        haveFoundAParent -= parent.end - parent.start + 1;
    }

    public Block overlap(Block b) {
//        System.out.println(this.start + "-" + this.end + " " + b.start + "-" + b.end + " " + !((this.start > b.end + 1) || (b.start > this.end + 1)));
        if (!((this.start > b.end) || (b.start > this.end))) {
            return new Block(0, 0, Math.max(this.start, b.start), Math.min(b.end, this.end), 0);
//            coalesceRange.start = Math.max(this.start, b.start);
//            coalesceRange.end = Math.min(b.end, this.end);
        } else {
            return null;
        }
    }

    public void reUse(int ID, int gen, int start, int end, int descendants, boolean isCoalescent, Block whoIRepresent, ArrayList<Block> parents, ArrayList<Block> children) {
        this.start = start;
        this.end = end;
        this.numberOfDescendants = descendants;
        this.isCoalescent = isCoalescent;
        this.representedChildren = whoIRepresent;
        this.parents = parents;
        this.children = children;
    }

    public int compareTo(Object b) {
        return this.start - ((Block) b).start;
    }

    public String getRangeString() {
        return "gen" + gen + ".ID" + ID + ".R" + start + "-" + end + "_D" + numberOfDescendants + "_C" + ((isCoalescent) ? 1 : 0);
    }
    
}

class BlockGenComparator implements Comparator<Block> {

    @Override
    public int compare(Block b1, Block b2) {
        return b1.gen - b2.gen;
    }
}
