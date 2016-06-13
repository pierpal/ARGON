package ARGON;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.TreeMap;

/**
 *
 * @author Pier Palamara <ppalama@hsph.harvard.edu>
 */
public class IBD {

    float threshold;
    int samples;
    String[][] ancestors;
    Integer[][] ancStart;
    TreeMap<Integer, Integer> map = new TreeMap<Integer, Integer>();

    public IBD(float threshold, int samples) {
        this.threshold = threshold;
        this.samples = samples;
        ancestors = new String[samples][samples];
        ancStart = new Integer[samples][samples];
        for (int i = 0; i < samples; i++) {
            for (int j = 0; j < samples; j++) {
                ancestors[i][j] = "";
                ancStart[i][j] = 1;
            }
        }
    }

    public StringBuilder addIBDset(int from, int to, int fromGen, int toGen, int gen, long ID, ArrayList<DescendantsList> sets) throws Exception {
        StringBuilder IBDstring = new StringBuilder();
        map.put(fromGen, from);
        map.put(toGen, to);
        DescendantsList currentSetList = sets.get(0);
        HashSet<Integer> currentSet = currentSetList.getHashSetClone();
        for (int i = 1; i < sets.size(); i++) {
            DescendantsList nextSetList = sets.get(i);
            HashSet<Integer> nextSet = nextSetList.getHashSetClone();
            for (Integer j : currentSet) {
                for (Integer k : nextSet) {
                    int id1, id2;
                    if (k > j) {
                        id1 = j - 1;
                        id2 = k - 1;
                    } else {
                        id1 = k - 1;
                        id2 = j - 1;
                    }
                    String anc = ancestors[id1][id2];
                    Integer start = ancStart[id1][id2];
                    if (anc.compareTo("") == 0) {
                        ancestors[id1][id2] = gen + "\t" + ID;
                        ancStart[id1][id2] = fromGen;
                    } else if (anc.compareTo(gen + "\t" + ID) != 0) {
                        int len = fromGen - start;
                        ancestors[id1][id2] = gen + "\t" + ID;
                        ancStart[id1][id2] = fromGen;
                        if (len / 1000000.0f >= threshold) {
                            IBDstring.append(new StringBuilder((id1 + 1) + "\t" + (id2 + 1) + "\t" + map.get(start) + "\t" + map.get(fromGen - 1) + "\t" + start + "\t" + (fromGen - 1) + "\t" + len / 1000000.0 + "\t" + anc + "\n"));
                        }
                    }
                }
            }
            currentSet.addAll(nextSet);
        }
        return IBDstring;
    }

    public StringBuilder finalize(int toGen) {
        StringBuilder IBDstring = new StringBuilder();
        for (int i = 0; i < samples; i++) {
            for (int j = i + 1; j < samples; j++) {
                String anc = ancestors[i][j];
                Integer start = ancStart[i][j];
                int len = toGen - start + 1;
                if (len / 1000000.0f >= threshold) {
                    IBDstring.append(new StringBuilder((i + 1) + "\t" + (j + 1) + "\t" + map.get(start) + "\t" + map.get(toGen) + "\t" + start + "\t" + toGen + "\t" + len / 1000000.0 + "\t" + anc + "\n"));
                }
            }
        }
        return IBDstring;
    }

}
