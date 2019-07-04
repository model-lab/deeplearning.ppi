package cn.modellab.deeplearning.ppi;


import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.indexing.NDArrayIndex;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.util.ModelSerializer;
import org.nd4j.evaluation.classification.Evaluation;
import org.nd4j.evaluation.classification.ROC;
import org.nd4j.linalg.io.ClassPathResource;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

interface Predicate<T> {
    boolean test(T t);
}

public class testModel {
    static INDArray boolOp(INDArray arr1, Predicate<Boolean []> predicate) {
        INDArray result = Nd4j.create(arr1.columns(),1);

        //
        for (int i = 0; i < arr1.length(); i++) {
            boolean answer = predicate.test(new Boolean[]{arr1.getDouble(i) == 1.0});
            result.putScalar(i, answer ? 1.0 : 0.0);
        }
        return result;
    }

    protected static final Logger log = LoggerFactory.getLogger(testModel.class);
    public static void main(String[] args) throws Exception {
        int modelNum=32;
        INDArray readCSV = Nd4j.readNumpy(new ClassPathResource("datasets/Colland04_test_dataset.csv").getFile().getPath(),",");
        int nn = readCSV.rows();
        INDArray Data = readCSV.get(NDArrayIndex.all(),NDArrayIndex.interval(1, readCSV.columns()));
        INDArray Label = readCSV.get(NDArrayIndex.all(),NDArrayIndex.point(0));
        INDArray equal0 = boolOp(Label.eq(0), new Predicate<Boolean[]>() {
            @Override
            public boolean test(Boolean[] booleans) {
                return booleans[0];
            }
        });
        INDArray testDataLabel = Nd4j.concat(1, equal0, Label.reshape(Label.columns(),1));
        INDArray testData = Data.reshape('c', new int[]{nn, 6, 1, 343});
        Nd4j.getRandom().setSeed(12345);
        INDArray v = Nd4j.zeros(nn,2);
        for(int i=0;i<modelNum;i++){
            MultiLayerNetwork modelTest=ModelSerializer.restoreMultiLayerNetwork(new ClassPathResource("models/model_"+i+".net").getFile(),true);
            INDArray v1 = modelTest.output(testData);
            v= v1.add(v);
            log.info("Evaluation {} finish...",i);
        }
        INDArray vv=v.div(modelNum);
        Evaluation eval = new Evaluation(2);
        eval.eval(testDataLabel,vv);
        log.info("Evaluation results...");
        log.info(eval.stats());

        ROC evalROC = new ROC();
        evalROC.eval(testDataLabel,vv);
        log.info(evalROC.stats());
    }
}

