package cn.modellab.deeplearning.ppi;

import org.datavec.api.records.reader.RecordReader;
import org.datavec.api.records.reader.impl.csv.CSVRecordReader;
import org.datavec.api.split.FileSplit;
import org.deeplearning4j.datasets.datavec.RecordReaderDataSetIterator;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.inputs.InputType;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.LSTM;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.deeplearning4j.optimize.listeners.ScoreIterationListener;
import org.deeplearning4j.util.ModelSerializer;
import org.nd4j.evaluation.classification.Evaluation;
import org.nd4j.linalg.activations.Activation;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.io.ClassPathResource;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.nd4j.linalg.learning.config.Adam;
import org.nd4j.linalg.lossfunctions.LossFunctions;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.util.Random;

public class trainModel {
    protected static final Logger log = LoggerFactory.getLogger(trainModel.class);
    protected static int channels = 6;
    protected static int batchSize = 64;

    protected static long seed = 42;
    protected static int epochs = 150;
    protected static int labelIndex = 0;
    protected static int numClasses = 2;
    protected static int nrow = 1;
    protected static int ncol = 343;
    protected static Random randNumGen = new Random(seed);

    public static void main(String[] args) throws Exception {
        int numLinesToSkip = 0;
        char delimiter = ',';

        DataSetPreProcessor preProcessor = new DataSetPreProcessor() {
            //@Override
            public void preProcess(org.nd4j.linalg.dataset.api.DataSet toPreProcess) {
                INDArray d = toPreProcess.getFeatures();
                toPreProcess.setFeatures(d.reshape('c', new int[]{d.rows(), channels, nrow, ncol}));
                toPreProcess.shuffle();
            }
        };

        File train_csvFile = new ClassPathResource("datasets/train_dataset_0.csv").getFile();
        File test_csvFile = new ClassPathResource("datasets/Colland04_test_dataset.csv").getFile();

        RecordReader recordReaderTrain = new CSVRecordReader(numLinesToSkip, delimiter);
        RecordReader recordReaderTest = new CSVRecordReader(numLinesToSkip, delimiter);
        recordReaderTrain.initialize(new FileSplit(train_csvFile));
        recordReaderTest.initialize(new FileSplit(test_csvFile));
        DataSetIterator iteratorTrain = new RecordReaderDataSetIterator(recordReaderTrain, batchSize, labelIndex, numClasses);
        iteratorTrain.setPreProcessor(preProcessor);
        DataSetIterator iteratorTest = new RecordReaderDataSetIterator(recordReaderTest, batchSize, labelIndex, numClasses);
        iteratorTest.setPreProcessor(preProcessor);
        Adam.builder();

        MultiLayerConfiguration conf = new NeuralNetConfiguration.Builder()
                .seed(seed)
                .activation(Activation.RELU)
                .weightInit(WeightInit.XAVIER)
                .updater(new Adam(1.4792765037709115E-5,0.9,0.99,1E-8))
                .l2(1.978894473010271E-5)
                .list()
                .layer(0, new DenseLayer.Builder().nIn(2058).nOut(909).dropOut(0.9)
                        .build())
                .layer(1, new DenseLayer.Builder().nOut(396).dropOut(0.9)
                        .build())
                .layer(2, new DenseLayer.Builder().nOut(396).dropOut(0.9)
                        .build())
                .layer(3, new DenseLayer.Builder().nOut(427).dropOut(0.9)
                        .build())
                .layer(4, new DenseLayer.Builder().nOut(226)
                        .build())
                .layer(5, new LSTM.Builder().nOut(226)
                        .build())
                .layer(6, new OutputLayer.Builder(LossFunctions.LossFunction.NEGATIVELOGLIKELIHOOD)
                        .activation(Activation.SOFTMAX)
                        .nOut(2).build())
                .setInputType(InputType.convolutional(nrow, ncol, channels))
                .build();

        MultiLayerNetwork net = new MultiLayerNetwork(conf);
        net.init();
        net.setListeners(new ScoreIterationListener(100));
        log.debug("Total num of params: {}", net.numParams());

        double maxAcc = 0;
        int maxEpoch = 0;
        for (int n = 0; n < epochs; n++) {
            net.fit(iteratorTrain);
            log.info(String.format("Epoch %d finished training", n + 1));
            log.info("Evaluate model....");
            Evaluation eval = net.evaluate(iteratorTest);

            log.info(eval.stats());
            if (eval.accuracy() > maxAcc) {
                maxEpoch = n;
                ModelSerializer.writeModel(net, "/Users/cx/Desktop/model/model.net", true);
                epochs = n + 30;
            }
            log.info("Max accuracy of model....");
            log.info(maxEpoch + "\t" + maxAcc);
        }
    }
}
