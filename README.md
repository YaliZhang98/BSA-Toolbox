# Sequence Analysis Toolbox

## Introduction
In this program, Flask, a popular framework in python is used to build a web application. This application is a toolbox that can get a uploaded fasta or fastq file or a pure nucleotide sequence submitted by users, and then performs different analysis and functions on nucleotide sequences.  

There are five different functions in this toolbox:  
1) GC content calculation. In this function, GC content of each sequence in the uploaded file will be calculated separately. In additon, the total GC content of all nucleotide sequences in file will be calculated as well.  
2) Sequence length detection. In this function, the maximum and minimum length of sequences in uploaded file will be found, and the specific sequence with the longest and shortest length also will be pointed at the same time. In addition, the average length of all sequences in file will be calculated.  
3) Reversed complement of nucleotide sequences. In this function, reversed complementary sequence of each nucleotide sequence in the uploaded file will be found respectively. All reversed complements will be output in a fasta format file.  
4) Translation of nucleotide sequence. In this function, all nucleotide sequence in the uploaded file will be translated to corresponding amino acid sequence separately. All amino acid sequences will be output in a fasta format file.  
5) GC content plot. In this function, all sequences in the upload file will be merged into one sequence, and the GC content changes will be ploted over this one sequence. The window size and step will be defined by user.  

## Files included
In this application, there are nultiple files:  
1) A python script that can activate web application.  
2) A readme file that contains the introduction and usage method of this application.  
3) A templates folder that contains all html scripts.  
4) A static folder that will contains the latest output results file on the web application. Existing files are reference output of the example files in example_file folder. The files in this folder will be overwritten after the user starts application.  
5) A example_file that contains one fasta file and one fastq file that can be used to do the test on this web application.  

## Usage

### Download or clone application to your own directory
```shell
# clone application through terminal
git clone https://github.com/YaliZhang98/Binp29_final_project.git
```

### Installation of packages used in application
Here will use conda to install some packages. If you have not installed conda,  follow commands below to install it.  

#### Conda installation
```shell
wget https://repo.continuum.io.miniconda/Miniconda3-latest-Linux-x86_64.sh -o ~miniconda.sh # download conda (conda 4.9.2)
ls # whether there is a file named miniconda*.sh
bash Miniconda3-latest-Linux-x86_64.sh # Run the bash program
cd miniconda3/
cd etc
cd profile.d/
pwd # find the route of conda installed
export PATH=route:$PATH # add the route into your ~/.bashrc folder so that we can use it anywhere.
source ~/.bashrc

conda -V # it is success if there is a version of conda
```
#### Packages installation
Here will install Flask used in python for building web applications. The  numpy required for calculation in python. And the matplotlib needed for plotting in python.  
```shell
# install the Flask package in a virtual environment using conda
conda create -n flaskenv python=3.6 flask

# activate virtual flaskenv environment
conda activate flaskenv 

# check whether flaskenv package was installed successfully
flask --version
#Python 3.6.13
#Flask 1.1.2
#Werkzeug 1.0.1

# install package will be used in python script
pip install numpy
pip install matplotlib
```

### Start web application
```shell
python sequence_analysis_toolbox.py
# Press "Ctrl + click" to open the web of application
# Press "Ctrl + C" to quit
```

## Example
There are one fasta file and one fastq file in example_file folder can be used to test the program.  
1) After starting the web application, choose one of functions in toolbox. Then manually enter a pure nucleotide sequence or upload a fasta or fastq file.  
2) Press "Reset" botton can clear the data you entered before. Press "Submit" botton can submit your data and then jump to the result page.  
3) On result web page, the results will be printed on the web page directly. There is also a botton above the output results called "download file" or "download figure" that enable user download results in a file and save it on the local computer.  




