package ARGON;

import java.util.BitSet;
import java.util.HashSet;

/**
 *
 * @author Pier Palamara <ppalama@hsph.harvard.edu>
 */
public final class DescendantsList {

    HashSet<Integer> list = new HashSet<Integer>();
    BitSet set;

    public static int numSamples = 0;
    public static boolean onlyStoreNumber = false;

//    static int created = 0;
//    static int killed = 0;
//    @Override
//    protected void finalize() { //remove
//        killed++; //remove
//    }
    public DescendantsList() {
        if (onlyStoreNumber) {
            // NOT IMPLEMENTED
            list.add(0);
        } else {

        }
//        created++;
    }

    public DescendantsList(int ID) {
//        created++;
        if (onlyStoreNumber) {
            // NOT IMPLEMENTED
            list.add(1);
        } else {
            list.add(ID);
        }
//        System.out.println("created " + getBinaryList() + " " + numSamples);
    }

    public DescendantsList(DescendantsList l1, DescendantsList l2) {
//        created++;
        if (onlyStoreNumber) {
            list.clear();
            list.add(l1.getSize() + l2.getSize());
        } else {
            list.addAll(l1.list);
            list.addAll(l2.list);
        }
    }

    public DescendantsList(DescendantsList l) {
//        created++;
        if (onlyStoreNumber) {
            int s = this.getSize();
            list.clear();
            list.add(l.getSize() + s);
        } else {
            list.addAll(l.list);
        }
    }

    public void add(DescendantsList l) {
        if (onlyStoreNumber) {
            int s = this.getSize();
            list.clear();
            list.add(l.getSize() + s);
        } else {
            if (l.list == null || list == null) {
                BitSet other;
                if (l.list == null) {
                    other = l.set;
                } else {
                    other = l.transformToBitSet();
                }
                if (list != null) {
                    set = transformToBitSet();
                    list = null;
                }
                set.or(other);
            } else {
                this.list.addAll(l.list);
                if (numSamples < 8 * this.list.size()) {
                    set = transformToBitSet();
                    list = null;
                }
            }
        }
    }

    public StringBuffer getBinaryListWithSpaces() {
        if (onlyStoreNumber) {
            // NOT IMPLEMENTED
            return null;
        } else {
            if (list != null) {
                StringBuffer listBuff = new StringBuffer();
                for (int i = 1; i <= numSamples; i++) {
                    listBuff.append((list.contains(i) ? 1 : 0));
                    listBuff.append(" ");
                }
                if (listBuff.length() > 0) {
                    listBuff.deleteCharAt(listBuff.length() - 1);
                }
                return listBuff;
            } else {
                StringBuffer listBuff = new StringBuffer();
                for (int i = 0; i < numSamples; i++) {
                    listBuff.append((set.get(i)) ? 1 : 0);
                    listBuff.append(" ");
                }
                if (listBuff.length() > 0) {
                    listBuff.deleteCharAt(listBuff.length() - 1);
                }
                return listBuff;
            }
        }
    }

    public BitSet transformToBitSet() {
        if (onlyStoreNumber) {
            // NOT IMPLEMENTED
            return null;
        } else {
            BitSet newSet = new BitSet(numSamples);
            for (Integer i : list) {
                newSet.set(i - 1);
            }
            return newSet;
        }
    }

    public StringBuffer getBinaryList() {
        if (onlyStoreNumber) {
            // NOT IMPLEMENTED
            return null;
        } else {
            if (list != null) {
                StringBuffer listBuff = new StringBuffer();
                for (int i = 1; i <= numSamples; i++) {
                    listBuff.append((list.contains(i) ? 1 : 0));
                }
                return listBuff;
            } else {
                StringBuffer listBuff = new StringBuffer();
                for (int i = 0; i < numSamples; i++) {
                    listBuff.append((set.get(i)) ? 1 : 0);
                }
                return listBuff;
            }
        }
    }

    public StringBuffer getSetList() {
        if (onlyStoreNumber) {
            // NOT IMPLEMENTED
            return null;
        } else {
            if (list != null) {
                StringBuffer listBuff = new StringBuffer();
                for (Integer i : list) {
                    listBuff.append(i);
                    listBuff.append(" ");
                }
                return listBuff;
            } else {
                StringBuffer listBuff = new StringBuffer();
                for (int i = 0; i < numSamples; i++) {
                    if (set.get(i)) {
                        // samples start at 1, so add 1
                        listBuff.append(i + 1);
                        listBuff.append(" ");
                    }
                }
                return listBuff;
            }
        }
    }

    public HashSet<Integer> getHashSetClone() {
        if (onlyStoreNumber) {
            // NOT IMPLEMENTED
            return null;
        } else {
            HashSet<Integer> l = new HashSet<Integer>();
            if (list != null) {
                l.addAll(list);
                return l;
            } else {
                for (int i = 0; i < numSamples; i++) {
                    if (set.get(i)) {
                        // samples start at 1, so add 1
                        l.add(i + 1);
                    }
                }
                return l;
            }
        }
    }

    public int getSize() {
        if (onlyStoreNumber) {
            int ret = 0;
            for (int i : list) {
                ret = i;
            }
            return ret;
        } else if (list != null) {
            return list.size();
        } else {
            return set.cardinality();
        }
    }

}
