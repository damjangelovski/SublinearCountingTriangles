import org.jgrapht.VertexFactory;

class IntVertexFactory implements VertexFactory {

    Integer count;

    public IntVertexFactory() {
        count = 0;
    }


    public Integer createVertex() {
        return count++;
    }
}
