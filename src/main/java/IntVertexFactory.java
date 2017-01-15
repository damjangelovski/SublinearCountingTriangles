import org.jgrapht.VertexFactory;

class IntVertexFactory implements VertexFactory {

    Integer count;

    public IntVertexFactory() {
        count = 0;
    }

    @Override
    public Integer createVertex() {
        return count++;
    }
}
