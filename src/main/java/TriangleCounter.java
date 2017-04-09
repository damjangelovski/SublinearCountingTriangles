import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import org.apache.commons.cli.*;
import org.ejml.simple.SimpleMatrix;
import org.jgrapht.*;
//import edu.uci.ics.jung.*;
import org.jgrapht.generate.RandomGraphGenerator;
import org.jgrapht.graph.*;
import com.google.common.math.*;

import java.util.*;

public class TriangleCounter {

    private static final double c = 1;// for number of iterations, not discussed
    private static final double c1 = 2000; //for s1, TODO: See Theorem 4 and proof, pg. 8/9 counting triangles
    private static final double c2 = 2000; //for s2, sufficiently large, TODO: See pg. 14 counting triangles
    private static final double ch = 100; //TODO: See Claim 2 proof, pg. 8 counting triangles

    private static final double cd = 1;  //for beta; cd>1; TODO: See pg.5 "Approximating Average Parameters of Graphs"

    private static Random rand;
    private static long startTime;

    private static UndirectedGraph<Integer, DefaultEdge> g;
    private static int n;

    private static int s_set,s_iters,s_heavy;

    public static void main(String [] args){
        Options options = new Options();

        Option n = new Option("n", "n verts", true, "n verts ");
        n.setRequired(true);
        options.addOption(n);

        Option e = new Option("e", "edges count", true, "no edges");
        e.setRequired(true);
        options.addOption(e);

        Option  corr_perc=
                new Option("corr_p",
                        "corr_percent",
                        true, "corr_percent");
        corr_perc.setRequired(true);
        options.addOption(corr_perc);

        Option  err_perc=
                new Option("err_p",
                        "err_p",
                        true, "err_p");
        err_perc.setRequired(true);
        options.addOption(err_perc);


        Option  sset=
                new Option("s_st",
                        "s_st",
                        true, "s_st");
        sset.setRequired(true);
        options.addOption(sset);

        Option  siter=
                new Option("s_iter",
                        "s_iter",
                        true, "s_iter");
        siter.setRequired(true);
        options.addOption(siter);

        Option  shev=
                new Option("s_heavy",
                        "s_heavy",
                        true, "s_heavy");
        shev.setRequired(true);
        options.addOption(shev);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException pe) {
            System.out.println(pe.getMessage());
            formatter.printHelp("utility-name", options);

            System.exit(1);
            return;
        }

        int nV=Integer.parseInt(cmd.getOptionValue("n")),
                nE=Integer.parseInt(cmd.getOptionValue("e"));

        int s_st=Integer.parseInt(cmd.getOptionValue("s_st")),
                s_iter=Integer.parseInt(cmd.getOptionValue("s_iter")),
                s_heav=Integer.parseInt(cmd.getOptionValue("s_heavy"));


        s_set=s_st;
        s_iters=s_iter;
        s_heavy=s_heav;

        double corr_p=Double.parseDouble(cmd.getOptionValue("corr_p")),
                err_p=Double.parseDouble(cmd.getOptionValue("err_p"));


        mainRuning(nV,nE,corr_p,err_p,false);

    }
    public static void mainRuning(int nV,int nE,double corr_p,double err_p,boolean doPrint) {


        rand = new Random();
        startTime = System.currentTimeMillis();

//        g = generateGraph(31000,351202);


//        g = generateGraph(750,217340);
        g = generateGraph(nV,nE);
        System.out.println("graph generated in: " + (System.currentTimeMillis() - startTime) + " milliseconds");
        n = g.vertexSet().size();

//        startTime = System.currentTimeMillis();
//        System.out.println("Average degree is: " + getAvgDeg());
//        System.out.println("Linear average degree finished in: " + (System.currentTimeMillis() - startTime) + " milliseconds");
//        startTime = System.currentTimeMillis();
//        System.out.println("Approximate average degree with e=0.5 is: " + averageDegreeApproximation(0.5));
//        System.out.println("Sub-linear average degree finished in: " + (System.currentTimeMillis() - startTime) + " milliseconds");


//        startTime = System.currentTimeMillis();
//        System.out.println("Triangle estimate with e=1 is: " + estimate(0.1));
//        System.out.println("Sub-linear triangle count finished in: " + (System.currentTimeMillis() - startTime) + " milliseconds");



        startTime = System.currentTimeMillis();
        int tri_c=getTri(doPrint);

        System.out.println("Real tri count: "+tri_c);
        System.out.println("Real Time Complexity: " + (System.currentTimeMillis() - startTime) + " milliseconds");

        startTime = System.currentTimeMillis();
        double d= getAvgDeg();
        double m = n * d / 2;
        int tri_c_corrupt= (int) (tri_c+(tri_c*makeRand()*corr_p));
        System.out.println("Tri corrupt_num is: " +tri_c_corrupt);

        startTime=System.currentTimeMillis();
        System.out.println("Estimate tri count: " +
                estimateWithAdvice(m,tri_c_corrupt,err_p,doPrint));
        System.out.println("Our Time Complexity: " + (System.currentTimeMillis() - startTime) + " milliseconds");




    }
    public static double makeRand(){
        return (rand.nextDouble()*2)-1;
    }
    //not sub-linear
    public static double getAvgDeg() {
        return g.edgeSet().size() * 2.0 / g.vertexSet().size() * 1.0;
    }
    public static int getTri(boolean doP){
        SimpleMatrix adj=new SimpleMatrix(n,n);
        int eCn=0;
        for (DefaultEdge de :
                g.edgeSet()) {
           int v=g.getEdgeSource(de),u=g.getEdgeTarget(de);
           adj.set(v,u,1);
           adj.set(u,v,1);
        if (doP && eCn%100000==12){
            System.out.println("Ecn is: "+eCn);
        }
        eCn++;
//            break;

        }


        SimpleMatrix trip=adj.mult(adj).mult(adj);


        return (int)trip.trace()/6;
    }
    //not sub-linear
    public static int getNumTriangles() {
        int count = 0;
        for (Integer v : g.vertexSet()) {
            if (v%100==11){
                System.out.println("Im in: "+v);
            }
            for (DefaultEdge e : g.edgesOf(v)) {
                Integer v1 = g.getEdgeSource(e);
                if (v == v1) {
                    v1 = g.getEdgeTarget(e);
                }
                for (DefaultEdge e2 : g.edgesOf(v)) {
                    if (e2 == e) {
                        continue;
                    }
                    Integer v2 = g.getEdgeSource(e2);
                    if (v == v2) {
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
        int t = (int) Math.ceil(Math.log(n) / Math.log1p(1 + beta)) + 1;

        Set<Integer> s = new HashSet<Integer>((int) k);
        Map<Integer, Set<Integer>> b = new HashMap<Integer, Set<Integer>>();
        Integer[] allV = g.vertexSet().toArray(new Integer[0]);

        for (int j = 1; j <= k; j++) {
            int v = allV[rand.nextInt(n)];
            s.add(v);
            int dv = g.edgesOf(v).size();
            int i = (int) Math.ceil(Math.log(dv) / Math.log1p(1 + beta));
            if (b.get(i) == null) {
                b.put(i, new HashSet<Integer>());
            }
            b.get(i).add(v);
        }
        for (Map.Entry<Integer, Set<Integer>> entry : b.entrySet()) {
            if (!(entry.getValue().size() / (1.0 * s.size()) >= Math.sqrt(e * l / (6 * n)) / t)) {
                b.remove(entry.getKey());
            }
        }
        Set<Integer> L=b.keySet();
        double sum=0;
        for (Map.Entry<Integer, Set<Integer>> entry : b.entrySet()) {
            double ai=0;
            for (Integer v : entry.getValue()) {
                DefaultEdge[] edgesV = g.edgesOf(v).toArray(new DefaultEdge[0]);
                DefaultEdge randE = edgesV[rand.nextInt(edgesV.length)];
                Integer u = getOtherEndpoint(randE, v);
                int du = g.edgesOf(v).size();
                Integer j=(int) Math.ceil(Math.log(du) / Math.log1p(1 + beta));
                if(!L.contains(j)){
                    ai++;
                }
            }
            ai/=entry.getValue().size();
            sum+=(1+ai)*entry.getValue().size()*Math.pow(1+beta,entry.getKey());
        }
        return sum/k;
    }

    //counting triangles, page 15
    public static double estimate(double e) {
        double e1 = e;// / (3 * ch);
        double d = getAvgDeg(); //TODO: change to sub-linear algorithm.
        double m = n * d / 2;
        double t_chengel = n * n * n;
        int t_ch_counter=1;
        while (t_chengel >= 1) {

            double t_iner=n*n*n;
            for (int iner=0;iner<t_ch_counter;iner++) {



                double x = n * n * n;
                for (int i = 1; i <= c * Math.log(Math.log(n)) / e; i++) {

// UNF                   t_iner=((rand.nextDouble()-0.5)/5)

                    double xi = estimateWithAdvice(m, t_iner, e1,true);
                    //System.out.println("estimation with t="+t+ "xi="+xi);
                    if (x > xi) {
                        x = xi;
                    }
                }
                if (x >= t_iner) {
                    return x;
                }
                t_iner/=2.0;
            }

            t_chengel = t_chengel / 2.0;
            t_ch_counter+=1;
        }
        return -1;
    }

    //counting triangles, page 12
    private static double estimateWithAdvice(double m, double t, double e,boolean doP) {
        double s1 = c1 * Math.pow(e, -3) * Math.log(n / e) * (n / Math.pow(t, 1.0 / 3));
        double s2 = c2 * Math.pow(e, -4) * Math.pow(Math.log(n), 2) * (Math.pow(m, 3.0 / 2) / t);

        s1=s_set;
        s2=s_iters;
//        ArrayList<Integer> counter=new ArrayList<Integer>();
        Multiset<Integer> counter = HashMultiset.create();
        //System.out.println("estimateWithAdvice started t="+t+" s1="+s1+" s2="+s2);
        S s = new S(g, (int) s1, rand);
        double y = 0;

        for (int i = 1; i <= s2; i++) {

            if (doP && i%500==11){
                System.out.println("ind is: "+i+" of "+s2);
            }


            Integer v = s.getRandVertex();
            if (isHeavy_Kiko(v, m, e, t)) {
                continue;
            }
            DefaultEdge[] edgesV = g.edgesOf(v).toArray(new DefaultEdge[0]);
            DefaultEdge randE = edgesV[rand.nextInt(edgesV.length)];
            Integer u = getSmallerEndpoint(randE);
            Integer x = getOtherEndpoint(randE, v);
            DefaultEdge[] edgesU = g.edgesOf(u).toArray(new DefaultEdge[0]);
            Set<Integer> bannedVertexes = new HashSet<Integer>();
//            bannedVertexes.add(v);
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

            int isXLight = isHeavy_Kiko(x, m, e, t) ? 0 : 1;

            double z = 0;
            for (int j = 1; j <= r; j++) {
                DefaultEdge randE2 = edgesU[rand.nextInt(edgesU.length)];
                while (bannedVertexes.contains(getOtherEndpoint(randE2, u))) {
                    randE2 = edgesU[rand.nextInt(edgesU.length)];
                }
                Integer w = getOtherEndpoint(randE2, u);
//                bannedVertexes.add(w);
                int errC=unwrap(g.containsEdge(v, w), g.containsEdge(x, w) ,x == getSmallerEndpoint(x, w));
                counter.add( errC);



                if (g.containsEdge(v, w) && g.containsEdge(x, w) && x == getSmallerEndpoint(x, w)) {
                    int isWLight = isHeavy_Kiko(w, m, e, t) ? 0 : 1;

                    z += Math.max(edgesU.length, Math.sqrt(m)) / (1 + isXLight + isWLight);
                }

            }
            y += z / r;
        }


        Map<Integer, Double> myPercentiles =
                Quantiles.percentiles().indexes(10,20,40,70,95, 90, 99).compute(counter);

        if (doP) {
            System.out.println(myPercentiles);

            for (int indexC = 0; indexC < 8; indexC++) {
                System.out.print("Counter of : " + indexC + " is: " + counter.count(indexC) + " ; ");


            }
            System.out.println();
        }


        return (n / (double) ((int) s1 * (int) s2)) * s.getDegreeSum() * y;
    }
    private static int unwrap(boolean b1,boolean b2, boolean b3){
        return Byte.parseByte(String.format("%d%d%d",makeInt(b1),makeInt(b2),makeInt(b3)),2);

    }
    private static int makeInt(boolean b){
        return (b) ? 1 : 0;
    }

    //counting triangles, page 9
    private static boolean isHeavy_Kiko(Integer v, double m, double e, double t) {


//        double s=(c*Math.log(n/e)/Math.pow(e,3))*n/Math.pow(t,0.33333333);

        double s_2=4.0 * Math.pow(m, 3.0 / 2) / (t * e *100);
         s_2=s_heavy;

        DefaultEdge[] edgesV = g.edgesOf(v).toArray(new DefaultEdge[0]);
        if (edgesV.length > 2 * m / Math.pow(e * t, 1.0 / 3)) {
            return true;
        }
        double[] listX = new double[(int) Math.ceil( 10 * Math.log(n)) ];
        for (int i = 1; i <= 10 * Math.log(n); i++) {
            double Xi = 0;
            for (int j = 1; j <= s_2; j++) {
                double Yj = 0;
                DefaultEdge randE = edgesV[rand.nextInt(edgesV.length)];
                Integer u = getSmallerEndpoint(randE);
                Integer x = getOtherEndpoint(randE, v);
                DefaultEdge[] edgesU = g.edgesOf(u).toArray(new DefaultEdge[0]);
                Set<Integer> bannedVertexes = new HashSet<Integer>();
//                bannedVertexes.add(v);
                bannedVertexes.add(x);
                for (int k = 1; k <= edgesU.length / Math.sqrt(m); k++) {
                    DefaultEdge randE2 = edgesU[rand.nextInt(edgesU.length)];
                    while (bannedVertexes.contains(getOtherEndpoint(randE2, u))) {
                        randE2 = edgesU[rand.nextInt(edgesU.length)];
                    }
                    Integer w = getOtherEndpoint(randE2, u);
//                    bannedVertexes.add(w);
                    if (g.containsEdge(v, w) && g.containsEdge(x, w) && x == getSmallerEndpoint(x, w)) {
                        Yj += edgesU.length;
                    }
                }
                Yj = Yj / (edgesU.length);
                Xi += Yj;
            }
            Xi = Xi * edgesV.length / s_2;
            listX[i]=Xi;
        }


        double median = median(listX);
        if (median > Math.pow(t, 2.0 / 3) / Math.pow(e, 1.0 / 3)) {
            return true;
        }
        return false;
    }

    //counting triangles, page 9
//    private static boolean isHeavy(Integer v, double m, double e, double t) {
//        DefaultEdge[] edgesV = g.edgesOf(v).toArray(new DefaultEdge[0]);
//        if (edgesV.length > 2 * m / Math.pow(e * t, 1 / 3)) {
//            return true;
//        }
//        Set<Double> setX = new HashSet<Double>();
//        for (int i = 1; i <= 10 * Math.log(n); i++) {
//            double Xi = 0;
//            for (int j = 1; j <= 20 * Math.pow(m, 3 / 2) / (t * e * e); j++) {
//                double Yj = 0;
//                DefaultEdge randE = edgesV[rand.nextInt(edgesV.length)];
//                Integer u = getSmallerEndpoint(randE);
//                Integer x = getOtherEndpoint(randE, v);
//                DefaultEdge[] edgesU = g.edgesOf(u).toArray(new DefaultEdge[0]);
//                Set<Integer> bannedVertexes = new HashSet<Integer>();
//                bannedVertexes.add(v);
//                bannedVertexes.add(x);
//                for (int k = 1; k <= edgesU.length / Math.sqrt(m); k++) {
//                    DefaultEdge randE2 = edgesU[rand.nextInt(edgesU.length)];
//                    while (bannedVertexes.contains(getOtherEndpoint(randE2, u))) {
//                        randE2 = edgesU[rand.nextInt(edgesU.length)];
//                    }
//                    Integer w = getOtherEndpoint(randE2, u);
//                    bannedVertexes.add(w);
//                    if (g.containsEdge(v, w) && g.containsEdge(x, w) && x == getSmallerEndpoint(x, w)) {
//                        Yj += edgesU.length;
//                    }
//                }
//                Yj = Yj / (edgesU.length / Math.sqrt(m));
//                Xi += Yj;
//            }
//            Xi = Xi * edgesV.length / (20 * Math.pow(m, 3 / 2) / (t * e * e));
//            setX.add(Xi);
//        }
//        Double[] arrX = setX.toArray(new Double[0]);
//        Arrays.sort(arrX);
//        double median = median(arrX);
//        if (median > Math.pow(t, 2 / 3) / Math.pow(e, 1 / 3)) {
//            return true;
//        }
//        return false;
//    }
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

    private static UndirectedGraph<Integer, DefaultEdge> generateGraph(int v,int e) {
        RandomGraphGenerator<Integer, DefaultEdge> randomGenerator
                = new RandomGraphGenerator<Integer, DefaultEdge>(v,e);
        UndirectedGraph<Integer, DefaultEdge> randomGraph = new SimpleGraph<Integer, DefaultEdge>(DefaultEdge.class);
        VertexFactory<Integer> factory = new IntVertexFactory();
        randomGenerator.generateGraph(randomGraph, factory, null);

        return randomGraph;
    }

    public static double median(double[] m) {
        Arrays.sort(m);
        int middle = m.length / 2;
        if (m.length % 2 == 1) {
            return m[middle];
        } else {
            return (m[middle - 1] + m[middle]) / 2.0;
        }
    }
}
