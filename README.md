# deeplearning.ppi

>***In this repository, you may likely see test project - those content is only for special examples.***

## About

Protein-protein interactions (PPIs) are central to most biological processes. Although efforts have been devoted to the development of a methodology for predicting PPIs and protein interaction networks, the application of most existing methods is limited because they need information about protein homology or the interaction marks of the protein partners. In the present work, we propose a method for PPI prediction using only the property of protein sequences. This method was developed based on a deep learning framework combined with a K-means-based conjoint triad feature for describing protein sequences. The prediction ability of our approach is better than that of other sequence-based PPI prediction methods.

## Dependencies
Note that this repository contains all examples for test deeplearning.ppi project. It will download about 0.8g of dependencies from maven central when you are first starting out. That being said, this makes it easier to get started without worrying about what to download. 

This project is released under an Apache 2.0 license. By contributing code to this repository, you agree to make your contribution available under an Apache 2.0 license.

## Install and Run (On IntelliJ)

1. Type on the command line:

```
$ git clone https://github.com/model-lab/deeplearning.ppi.git
$ cd deeplearning.ppi
$ mvn clean install
```

2. Open IntelliJ and select import project. Then select the "deeplearning.ppi" home directory.
3. Select import projects from external models to ensure that Maven is selected.
4. Continue with the wizard options. Select the SDK that starts with the JDK. (you may need to click on the plus sign to see the relevant options... ) then click finish. Wait a moment, let IntelliJ download all dependencies.
5. A progress bar appears at the bottom right. Select a example from the tree directory on the left.
6. Right-click the file and run.

## Known issues with JavaFX

If you are running on JDK 1.7 or inferior, the maven-enforcer plugin will require you to set the variable JAVAFX_HOME before building. That variable should point to a directory containing jfxrt.jar, a file that is part of the JavaFX 2.0 distrbution.

Please set it to an instance of JavaFX that matches the JDK with which you are trying to use this project. Usually, the Sun JDK comes with JavaFX. However, OpenJDK does not and you may have to install OpenJFX, a free distribution of JavaFX.

Beware that your editor (e.g. IntelliJ) may not be using the JDK that is your system default (and that you may ancounter on the command line).
