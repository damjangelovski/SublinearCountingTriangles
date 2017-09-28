import org.jgrapht.*;
//import edu.uci.ics.jung.*;
import org.jgrapht.generate.RandomGraphGenerator;
import org.jgrapht.graph.*;

import java.util.*;

public class TriangleCounter {

    private static final double c = 1;// for number of iterations, not discussed
    private static final double c1 = 1; //for s1, TODO: See Theorem 4 and proof, pg. 8/9 counting triangles
    private static final double c2 = 1; //for s2, sufficiently large, TODO: See pg. 14 counting triangles
    private static final double ch = 287; //TODO: See Claim 2 proof, pg. 8 counting triangles

    private static final double cd = 100;  //for beta; cd>1; TODO: See pg.5 "Approximating Average Parameters of Graphs"

    private static Random rand;
    private static long startTime;

    private static UndirectedGraph<Integer, DefaultEdge> g;
    private static int n;

    public static void main(String[] args) {


        rand = new Random();
        startTime = System.currentTimeMillis();

        g = generateGraph();
        System.out.println("graph generated in: " + (System.currentTimeMillis() - startTime) + " milliseconds");
        n = g.vertexSet().size();

        startTime = System.currentTimeMillis();
        System.out.println("Average degree is: " + getAvgDeg());
        System.out.println("Linear average degree finished in: " + (System.currentTimeMillis() - startTime) + " milliseconds");
        startTime = System.currentTimeMillis();
        System.out.println("Approximate average degree with e=0.5 is: " + averageDegreeApproximation(0.5));
        System.out.println("Sub-linear average degree finished in: " + (System.currentTimeMillis() - startTime) + " milliseconds");

        /*startTime = System.currentTimeMillis();
        System.out.println("Triangle estimate with e=1 is: " + estimate(1));
        System.out.println("Sub-linear triangle count finished in: " + (System.currentTimeMillis() - startTime) + " milliseconds");
        startTime = System.currentTimeMillis();
        System.out.println("True No. of triangles is: "+getNumTriangles());
        System.out.println("superLinear triangle count finished in: " + (System.currentTimeMillis() - startTime) + " milliseconds");
*/
    }

    //not sub-linear
    public static double getAvgDeg() {
        return g.edgeSet().size() * 2.0 / g.vertexSet().size() * 1.0;
    }

    //not sub-linear
    public static int getNumTriangles() {
        int count = 0;
        for (Integer v : g.vertexSet()) {
            for (DefaultEdge e : g.edgesOf(v)) {
                Integer v1 = g.getEdgeSource(e);
                if (v.equals(v1)) {
                    v1 = g.getEdgeTarget(e);
                }
                for (DefaultEdge e2 : g.edgesOf(v)) {
                    if (e2 == e) {
                        continue;
                    }
                    Integer v2 = g.getEdgeSource(e2);
                    if (v.equals(v2)) {
                        v2 = g.getEdgeTarget(e2);
                    }
                    if (g.containsEdge(v2, v1)) {
                        count++;
                    }
                }
            }
        }
        return count / 6;
    }

    // algorithm defined by Goldreich and Ron
    public static double averageDegreeApproximation(double e) {
        double l = 1; //lower bound on average degree. With no previous information, put 1
        double beta=e/cd;
        double k = Math.sqrt(n / l) * Math.pow(e, -4.5) * Math.pow(Math.log(n), 2) * Math.log(1 / e);
        //double k = 500;
        int t = (int) Math.ceil(Math.log(n) / Math.log1p(beta)) + 1;

        Set<Integer> s = new HashSet<Integer>((int) k);
        Map<Integer, Set<Integer>> b = new HashMap<Integer, Set<Integer>>();
        Integer[] allV = g.vertexSet().toArray(new Integer[0]);

        System.out.println("size of sample is: "+k);
        System.out.println("beta is: "+beta);
        System.out.println("t is: "+t);
        System.out.println("Max i is: "+Math.ceil(Math.log(n) / Math.log1p(beta)));
        System.out.println("Bucket percentage needs to be bigger than "+Math.sqrt(e * l / (6. * n)) / t);
        for (int j = 1; j <= k; j++) {
            int v = allV[rand.nextInt(n)];
            s.add(v);
            int dv = g.edgesOf(v).size();
            int i = (int) Math.ceil(Math.log(dv) / Math.log1p(beta));
            if (b.get(i) == null) {
                b.put(i, new HashSet<Integer>());
            }
            b.get(i).add(v);
        }
        int numNodes=0;
        for (Map.Entry<Integer, Set<Integer>> entry : b.entrySet()) {
            int i=entry.getKey();
            numNodes+=entry.getValue().size();
            System.out.println("for i="+i+" degree is between "+Math.ceil(Math.pow(1+beta,i-1))+" and "
                    +Math.floor(Math.pow(1+beta,i)) +" Bi size is: "+entry.getValue().size());
            if (!(entry.getValue().size() / (1.0 * s.size()) >= Math.sqrt(e * l / (6. * n)) / t)) {
                b.remove(entry.getKey());
                System.out.println("is deleted");
            }
        }
        System.out.println("number of nodes is "+numNodes);
        Set<Integer> L=b.keySet();
        double sum=0;
        for (Map.Entry<Integer, Set<Integer>> entry : b.entrySet()) {
            double ai=0;
            double bi=0;
            for (Integer v : entry.getValue()) {
                DefaultEdge[] edgesV = g.edgesOf(v).toArray(new DefaultEdge[0]);
                DefaultEdge randE = edgesV[rand.nextInt(edgesV.length)];
                Integer u = getOtherEndpoint(randE, v);
                int du = g.edgesOf(v).size();
                Integer j=(int) Math.ceil(Math.log(du) / Math.log1p(beta));
                if(!L.contains(j)){
                    ai++;
                }
                else {
                    bi++;
                }
            }
            System.out.println("for i: "+entry.getKey()+" ai:"+ai+" bi:"+bi);
            ai/=entry.getValue().size();
            sum+=(1+ai)*entry.getValue().size()*Math.pow(1+beta,entry.getKey());
        }
        //numNodes is the real size of the sample, as duplicate sample nodes are not counted
        //return sum/k;
        return sum/numNodes;
    }

    //counting triangles, page 15
    public static double estimate(double e) {
        double e1 = e / (3. * ch);
        double d = getAvgDeg(); //TODO: change to sub-linear algorithm.
        double m = n * d / 2.;
        double t = n * n * n;
        while (t >= 1) {
            double x = n * n * n;
            for (int i = 1; i < c * Math.log(Math.log(n)) / e; i++) {
                double xi = estimateWithAdvice(m, t, e1);
                //System.out.println("estimation with t="+t+ "xi="+xi);
                if (x > xi) {
                    x = xi;
                }
            }
            if (x >= t) {
                return x;
            }
            t = t / 2.;
        }
        return -1;
    }

    //counting triangles, page 12
    private static double estimateWithAdvice(double m, double t, double e) {
        double s1 = (c1 * Math.pow(e, -3) / Math.log(n / e)) * Math.log(n / Math.pow(t, 1 / 3.)); // log base is n/e
        double s11 = c1 * Math.pow(e, -3) * Math.log(n / e) / Math.log(n / Math.pow(t, 1 / 3.)); // log base is n/t^(1/3)
        double s12 = c1 * Math.pow(e, -3) * Math.log(n / e) * (n / Math.pow(t, 1 / 3.)); // log base is e
        double s2 = c2 * Math.pow(e, -4) * Math.pow(Math.log(n), 2) * (Math.pow(m, 3 / 2.) / t);
        double ch=binomial((int)Math.pow(t,1/3.)*12,3)/(t);
        System.out.println("estimateWithAdvice started t="+t+" s1="+s1+" s2="+s2+" s11="+s11+" s12="+s12+" new ch="+ch);


        S s = new S(g, 1+(int) s1, rand); //TODO: change second argument with the correct s1 formula
        double y = 0;

        for (int i = 1; i <= s2; i++) {
            Integer v = s.getRandVertex();
            if (isHeavy(v, m, e, t)) {
                continue;
            }
            DefaultEdge[] edgesV = g.edgesOf(v).toArray(new DefaultEdge[0]);
            DefaultEdge randE = edgesV[rand.nextInt(edgesV.length)];
            Integer u = getSmallerEndpoint(randE);
            Integer x = getOtherEndpoint(randE, v);
            DefaultEdge[] edgesU = g.edgesOf(u).toArray(new DefaultEdge[0]);
            Set<Integer> bannedVertexes = new HashSet<Integer>();
            bannedVertexes.add(v);
            bannedVertexes.add(x);

            int r;
            if (edgesU.length <= Math.sqrt(m)) {
                if (rand.nextDouble() < edgesU.length / Math.sqrt(m)) {
                    r = 1;
                } else {
                    //r = 0;
                    continue;
                }
            } else {
                r = (int) (edgesU.length / Math.sqrt(m));
            }

            int isXLight = isHeavy(x, m, e, t) ? 0 : 1;
            double z = 0;
            for (int j = 1; j <= r; j++) {
                DefaultEdge randE2 = edgesU[rand.nextInt(edgesU.length)];
                while (bannedVertexes.contains(getOtherEndpoint(randE2, u))) {
                    randE2 = edgesU[rand.nextInt(edgesU.length)];
                }
                Integer w = getOtherEndpoint(randE2, u);
                bannedVertexes.add(w);
                if (g.containsEdge(v, w) && g.containsEdge(x, w) && x.equals(getSmallerEndpoint(x, w))) {
                    int isWLight = isHeavy(w, m, e, t) ? 0 : 1;
                    z += Math.max(edgesU.length, Math.sqrt(m)) / (1. + isXLight + isWLight);
                }

            }
            y += z / r;
        }
        return (n / (double) ((int) s1 * (int) s2)) * s.getDegreeSum() * y;
    }

    //counting triangles, page 9
    private static boolean isHeavy(Integer v, double m, double e, double t) {
        DefaultEdge[] edgesV = g.edgesOf(v).toArray(new DefaultEdge[0]);
        if (edgesV.length > 2 * m / Math.pow(e * t, 1 / 3.)) {
            return true;
        }
        Set<Double> setX = new HashSet<Double>();
        for (int i = 1; i <= 10 * Math.log(n); i++) {
            double Xi = 0;
            for (int j = 1; j <= 20 * Math.pow(m, 3 / 2.) / (t * e * e); j++) {
                double Yj = 0;
                DefaultEdge randE = edgesV[rand.nextInt(edgesV.length)];
                Integer u = getSmallerEndpoint(randE);
                Integer x = getOtherEndpoint(randE, v);
                DefaultEdge[] edgesU = g.edgesOf(u).toArray(new DefaultEdge[0]);
                Set<Integer> bannedVertexes = new HashSet<Integer>();
                bannedVertexes.add(v);
                bannedVertexes.add(x);
                for (int k = 1; k <= edgesU.length / Math.sqrt(m); k++) {
                    DefaultEdge randE2 = edgesU[rand.nextInt(edgesU.length)];
                    while (bannedVertexes.contains(getOtherEndpoint(randE2, u))) {
                        randE2 = edgesU[rand.nextInt(edgesU.length)];
                    }
                    Integer w = getOtherEndpoint(randE2, u);
                    bannedVertexes.add(w);
                    if (g.containsEdge(v, w) && g.containsEdge(x, w) && x.equals(getSmallerEndpoint(x, w))) {
                        Yj += edgesU.length;
                    }
                }
                Yj = Yj / (edgesU.length / Math.sqrt(m));
                Xi += Yj;
            }
            Xi = Xi * edgesV.length / (20 * Math.pow(m, 3 / 2.) / (t * e * e));
            setX.add(Xi);
        }
        Double[] arrX = setX.toArray(new Double[0]);
        Arrays.sort(arrX);
        double median = median(arrX);
        if (median > Math.pow(t, 2 / 3.) / Math.pow(e, 1 / 3.)) {
            return true;
        }
        return false;
    }

    //counting triangles, page 7
    private static Integer getSmallerEndpoint(DefaultEdge e) {
        Integer v1 = g.getEdgeSource(e);
        Integer v2 = g.getEdgeTarget(e);
        return getSmallerEndpoint(v1, v2);
    }

    private static Integer getSmallerEndpoint(Integer v1, Integer v2) {
        int dv1 = g.edgesOf(v1).size();
        int dv2 = g.edgesOf(v2).size();

        if (dv1 == dv2) {
            if (v1 < v2) {
                return v1;
            } else {
                return v2;
            }
        }
        if (dv1 < dv2) {
            return v1;
        } else {
            return v2;
        }
    }

    private static Integer getOtherEndpoint(DefaultEdge e, Integer v) {
        int v2 = g.getEdgeSource(e);
        if (v2 != v) {
            return v2;
        } else {
            return g.getEdgeTarget(e);
        }
    }

    private static UndirectedGraph<Integer, DefaultEdge> createGraph() {
        UndirectedGraph<Integer, DefaultEdge> g =
                new SimpleGraph<Integer, DefaultEdge>(DefaultEdge.class);

        // add the vertices
        g.addVertex(1);
        g.addVertex(2);
        g.addVertex(3);
        g.addVertex(4);

        // add edges to create a circuit
        g.addEdge(1, 2);
        g.addEdge(2, 3);
        g.addEdge(3, 4);
        g.addEdge(4, 1);
        g.addEdge(3, 1);

        return g;
    }

    private static UndirectedGraph<Integer, DefaultEdge> generateGraph() {
        RandomGraphGenerator<Integer, DefaultEdge> randomGenerator
                = new RandomGraphGenerator<Integer, DefaultEdge>(1000, 100000);
        UndirectedGraph<Integer, DefaultEdge> randomGraph = new SimpleGraph<Integer, DefaultEdge>(DefaultEdge.class);
        VertexFactory<Integer> factory = new IntVertexFactory();
        randomGenerator.generateGraph(randomGraph, factory, null);

        return randomGraph;
    }

    public static double median(Double[] m) {
        int middle = m.length / 2;
        if (m.length % 2 == 1) {
            return m[middle];
        } else {
            return (m[middle - 1] + m[middle]) / 2.0;
        }
    }

    private static long binomial(int n, int k)
    {
        if (k>n-k)
            k=n-k;

        long b=1;
        for (int i=1, m=n; i<=k; i++, m--)
            b=b*m/i;
        return b;
    }
}
