import org.jgrapht.UndirectedGraph;
import org.jgrapht.graph.DefaultEdge;

import java.util.Random;

class S {
    Integer[] v;
    int[] dvSum;
    Random rand;

    public S(UndirectedGraph<Integer, DefaultEdge> g, int s1, Random rand) {
        v = new Integer[s1];
        dvSum = new int[s1];
        this.rand = rand;

        Integer[] allV = g.vertexSet().toArray(new Integer[0]);
        for (int i = 0; i < (int) s1; i++) {
            v[i] = allV[rand.nextInt(allV.length)];
        }
        dvSum[0] = g.edgesOf(v[0]).size();
        for (int i = 1; i < (int) s1; i++) {
            dvSum[i] = dvSum[i - 1] + g.edgesOf(v[i]).size();
        }
    }

    public Integer getRandVertex() {
        int low = 0, high = dvSum.length - 1;
        int pos = rand.nextInt(dvSum[high]);

        while (low + 1 < high) {
            int mid = low + (high - low) / 2;
            if (pos < dvSum[mid]) high = mid;
            else if (pos > dvSum[mid]) low = mid;
            else return v[mid];
        }
        return v[low];
    }

    public int getDegreeSum() {
        return dvSum[dvSum.length - 1];
    }
}
