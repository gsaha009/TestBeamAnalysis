# btdqmanalysis

Root based analysis framework for analyzing 2S module data collected during beam test. All the 
source code is in src directory. The main piece of code is BTAnalyzer.h and BTAnalyzer.cc. It 
reads relevant branches from the dqm input trees and does some basic plots. The main.cc is the 
program which drives this code. btdriver.py is the  python interface to the BTAnalyzer code.


#To checkout the code

originally cloned from here :  
```
https://github.com/sroychow/btdqmanalysis.git 
```

#Compile

```
cd TestBeamAnalysis
make
```

At this point a bin directory will be created with the executable and the libraries.

#execution 
#from bash
```
./bin/baseline.exe <inputfileName> <#of modules to be read from the tree> <outputfileName>
```
#python 
```
python btdriver.py <inputfileName> <#of modules to be read from the tree>
```
For this case, the output file name will be inputFilename_histos.root. 

##Notes
Please try to produce the output histo files with the names inputFilename_histos.root so that multi run harvesting of plots is easier.
compiler : gcc4.9 or higher.
ROOT: 6.16 or higher.
