# -*- coding: utf-8 -*-
"""
 Title: sequence_analysis_toolbox.py
 Date: 2021-03-15
 Author: Yali Zhang 
 
 Description: 
     In this script, Flask, a popular framework in python is used to build a web 
     application. This application is a toolbox that can get a uploaded fasta or 
     fastq file or a pure nucleotide sequence submitted by users, and then performs 
     different analysis and functions on nucleotide sequences. 

    There are five different functions in this toolbox:
    
    1) GC content calculation. In this function, GC content of each sequence in the 
    uploaded file will be calculated separately. In additon, the total GC content 
    of all nucleotide sequences in file will be calculated as well.
    
    2) Sequence length detection. In this function, the maximum and minimum length 
    of sequences in uploaded file will be found, and the specific sequence with 
    the longest and shortest length also will be pointed at the same time. In 
    addition, the average length of all sequences in file will be calculated.
    
    3) Reversed complement of nucleotide sequences. In this function, reversed 
    complementary sequence of each nucleotide sequence in the uploaded file will be 
    found respectively. All reversed complements will be output in a fasta format file.
    
    4) Translation of nucleotide sequence. In this function, all nucleotide sequence 
    in the uploaded file will be translated to corresponding amino acid sequence 
    separately. All amino acid sequences will be output in a fasta format file.
    
    5) GC content plot. In this function, all sequences in the upload file will be 
    merged into one sequence, and the GC content changes will be ploted over this one 
    sequence. The window size and step will be defined by user.
     
 Usage:
     python sequence_analysis_toolbox.py
     Ctrl+click to open the web page.
     Ctrl+c to stop the program
"""


from flask import Flask,request,render_template
from numpy import mean

app = Flask(__name__)


# Preprocessing functions of file ---------------------------------------------


# limitation of upload file format.
fasta_format = ['fasta','fna','fa']
fastq_format = ['fastq','fq']


# function that can extract the file format part, if it meets the requirements, return TRUE. 
# If it is not the specified file format, return FALSE.

def file_fasta(filename): # for fasta file
    return '.' in filename and filename.rsplit('.',1)[1] in fasta_format
 
def file_fastq(filename): # for fastq file
    return '.' in filename and filename.rsplit('.',1)[1] in fastq_format


# function that can extract sequence ID and nucleotide of each sequence in fasta or fastq file.

def sequence_fa(fasta):
    seq_dic = {}
    seq = ''
    seq_id = ''
    for line in fasta:
        line = line.decode('utf-8') # change the code of upload file 
        line = line.strip()
        if line.startswith(">"): # find id line
            if seq: # skip first sequence
                seq_dic[seq_id] = seq # add sequence and its id in dictionary
            seq_id = line # set new sequence id
            seq = '' # renew a new sequence
        else:
            line = line.upper()
            seq += line # parse sequence in one line
    seq_dic[seq_id] = seq # add last sequence and its id into dictionary
    return(seq_dic)

def sequence_fq(fastq):
    fq_seq_dic = {}
    fq_seq = ''
    fq_id = ''
    count = 0 # flag of number of line
    for line in fastq:
        line = line.decode('utf-8') # change the code of upload file 
        line = line.strip()
        count += 1 # count of each line
        
        if count%4 == 1: # find the id line
            if fq_seq: # skip first sequence
                fq_seq_dic[fq_id] = fq_seq # add sequence and its id in dictionary
            fq_id = line # set new sequence id
            fq_id = fq_id.replace('@',">") # change id to fasta format
            fq_seq = '' # renew a new sequence
            
        elif count%4 == 2: # find the sequence line
            line = line.upper() # change all nucleotide to uppercase
            fq_seq = line
            
    fq_seq_dic[fq_id] = fq_seq # add last sequence and its id into dictionary
    return(fq_seq_dic)


# function that can combine all the sequence in file into one sequence.

def to_one_seq(seq_dic):
    seq = ''
    for key in seq_dic.keys():
        seq += seq_dic[key] # add all sequence into one sequence
    return(seq)



# Toolbox catalog page --------------------------------------------------------

# the first web page - tollbox catalog page
@app.route('/')
def Toolbox():
    return render_template('Toolbox.html')



# GC-content calculation -----------------------------------------------------


# Function part 


## function that can calculate GC content of one sequence

def GC_content(seq):
    seq = seq.upper().strip() # change all nucleotide into uppercase
    length = len(seq) # find the length of sequence
    count = 0
    for nu in seq: # calculate number of G and C
        if nu == "C":
            count += 1
        if nu == "G":
            count += 1    
    GC_con = count/length*100 # calculate GC-content
    GC_con = '%.2f'%GC_con # set decimal points as two
    GC_content = str(GC_con) # change number to string
    return GC_content


## function that can calculate GC1, GC2 and GC3 content of one sequence

def GC123_content(seq):
    GC1_seq = ''
    GC2_seq = ''
    GC3_seq = ''
    for i in range(1,len(seq)+1):
        if i%3 == 1: # find nucleotide in position 1
            GC1_seq += seq[i-1]
        elif i%3 == 2: # find nucleotide in position 2
            GC2_seq += seq[i-1]
        else : # find nucleotide in position 3
            GC3_seq += seq[i-1]
    GC1_content = GC_content(GC1_seq) # calculate GC-content of these positions separately
    GC2_content = GC_content(GC2_seq)
    GC3_content = GC_content(GC3_seq)
    return [GC1_content, GC2_content, GC3_content]


## function that can combine all GC, GC1, GC2 and GC3 results of one sequence into a file

def GC_one_seq(seq):
    GC_c = GC_content(seq) # calculate GC-content
    GC123_c = GC123_content(seq) # calculate GC1, GC2, ad GC3 content
    content_all = []
    content = []
    content.append('one_seq') # add all GC-content result into one set
    content.append(GC_c)
    content.append(GC123_c[0])
    content.append(GC123_c[1])
    content.append(GC123_c[2])
    content_all.append(content)
    
    # print all results into a output file
    output = open('./static/GC_content.tsv','w')
    print('seqID\tGC(%)\tGC1(%)\tGC2(%)\tGC3(%)',file=output) # print header
    print('one_seq\t{}\t{}\t{}\t{}'.format(GC_c,GC123_c[0],GC123_c[1],GC123_c[2]),file=output)
    output.close()
    
    return content_all


## function that can combine all GC, GC1, GC2 and GC3 results of all sequence in one file

def GC_content_all(seq_dic):
    GC_dic = {}   
    for key in seq_dic.keys():
        seq = seq_dic[key] # extract each sequence
        GC_c = GC_content(seq) # calculate GC content
        GC123_c = GC123_content(seq)  # calculate GC1, GC2, ad GC3 content
        key = key.split()[0][1:] # extract sequence id
        GC_dic[key] = {}
        GC_dic[key]['GC'] = GC_c # extract GC content results of each sequence in a dictionary  
        GC_dic[key]['GC1'] = GC123_c[0]
        GC_dic[key]['GC2'] = GC123_c[1]
        GC_dic[key]['GC3'] = GC123_c[2]
    
    # create a list that contains all GC content results.
    content_all = []    
    for key in GC_dic.keys():
        content = []
        content.append(key)
        content.append(GC_dic[key]['GC'])
        content.append(GC_dic[key]['GC1'])
        content.append(GC_dic[key]['GC2'])
        content.append(GC_dic[key]['GC3'])
        content_all.append(content)
    
    # print all results into a output file
    output = open('./static/GC_content.tsv','a')
    for key in GC_dic.keys():
        print('{}\t{}\t{}\t{}\t{}'.format(key,GC_dic[key]['GC'],GC_dic[key]['GC1'],GC_dic[key]['GC2'],GC_dic[key]['GC3']),file=output)
    output.close()
        
    return content_all


# web part 

## function that can redirect to GC_content calculation web page.
@app.route('/GC-content')
def GC_content_input():
    return render_template('GC_content_input.html')

## function that can redirect to web page that can show results of GC-content calculation. 
@app.route('/GC-content', methods=['POST'])
def show_GC_content() :
    
    # specifiy the HTTP method to use when sending form-data
    if request.method == 'POST': 
        
        # get submit data of user
        seq = request.form.get('seq')
        file = request.files.get('file')   
        header = ['seqID','GC(%)','GC1(%)','GC2(%)','GC3(%)'] # header of output result 
       
        # if user upload a fasta file
        if file and file_fasta(file.filename): # check the file and file format
            if seq: # prevent conflicts if users submit sequence and file at the same time
                return render_template('GC_content_input.html') 
            else:
                seq_dic = sequence_fa(file)
                one_seq = to_one_seq(seq_dic)
                content = GC_one_seq(one_seq)
                content_2 = GC_content_all(seq_dic)
                content.extend(content_2)
                return render_template('GC_content.html', content=content, header=header)
        
    
        if file and file_fastq(file.filename): # check the file and file format
            if seq: # prevent conflicts if users submit sequence and file at the same time
                return render_template('GC_content_input.html') 
            else:            
                seq_dic = sequence_fq(file)
                one_seq = to_one_seq(seq_dic)
                content = GC_one_seq(one_seq)
                content_2 = GC_content_all(seq_dic)
                content.extend(content_2)
                return render_template('GC_content.html',content=content, header=header)
        
        if seq :
            if file: # prevent conflict if users submit sequence and file at the same time
                return render_template('GC_content_input.html') 
            else:            
                 content = GC_one_seq(seq)
                 return render_template('GC_content.html',content=content, header=header)
   
    return render_template('GC_content_input.html') # redirect to result web page


# Length detection-------------------------------------------------------------


# function part

## function that can find the average, maximum and minimum length of sequence.
def length_detection(seq_dic):
    length_dic = {}
    for key in seq_dic.keys():
        seq_id = key.split()[0][1:] # set new sequence id
        length_dic[seq_id] = len(seq_dic[key]) # add id and corresponding length into a dictionary 
    
    length_avg = mean(list(length_dic.values())) # find the average length
    length_max = max(length_dic.values()) # find the maximum length
    length_min = min(length_dic.values()) # find the minimum length
    
    # find ids whose length of sequence is the maximum or minimum 
    seq_length_max = ''
    seq_length_min = ''
    for key in length_dic.keys():
        if length_dic[key] == length_max: # find sequence whose length is the maximum
            if seq_length_max:
                seq_length_max += ','
                seq_length_max += key
            else:
                seq_length_max += key
        if length_dic[key] == length_min: # find sequence whose length is the minimum
            if seq_length_min:
                seq_length_min += ','
                seq_length_min += key
            else:
                seq_length_min += key
    
    # set two decimal points
    length_avg = '%.2f'%length_avg
    length_max = '%.2f'%length_max
    length_min = '%.2f'%length_min
    
    # set the structure of output results
    content = [['average_length',length_avg,'-']]
    content.extend([['maximum_length',length_max,seq_length_max]])
    content.extend([['minimum_length',length_min, seq_length_min]])
    
    # print all results into a output file
    output = open('./static/length_detection.tsv','w')
    print('character\tlength(pb)\tseqID',file=output)
    print('average_length\t{}\t-'.format(length_avg),file=output)
    print('maximum_length\t{}\t{}'.format(length_max,seq_length_max),file=output)
    print('minimum_length\t{}\t{}'.format(length_min, seq_length_min),file=output)
    output.close()
  
    return content


# web part

## redirect to length detection input web page
@app.route('/Length')
def length_input():
    return render_template('Length_detection_input.html')

## redirect to length detection results web page
@app.route('/Length', methods=['POST'])
def length() :
    
    # specifiy the HTTP method to use when sending form-data
    if request.method == 'POST':
        file = request.files.get('file') # get the upload file
        header = ['character','length(pb)','seqID'] # set header of result table
       
        # if user upload a fasta file
        if file and file_fasta(file.filename): # check the file and file format        
            seq_dic = sequence_fa(file) # get the sequence and id dictionary
            content = length_detection(seq_dic) # get the output result
            return render_template('Length_detection.html',content=content, header=header)
        
        # if user upload a fastq file
        if file and file_fastq(file.filename): # check the file and file format                
            seq_dic = sequence_fq(file)
            content = length_detection(seq_dic)
            return render_template('Length_detection.html',content=content, header=header)
   
    return render_template('Length_detection_input.html') 


# Reversed complements --------------------------------------------------------

# Step 1: complements
# step 2: reverse 

# function part

## function that can find the reversed complements sequence of a sequence
def comp_rev(seq):
    seq = seq.upper().strip()
    seq_comp = ''
    
    # find the complements sequence
    for nu in seq:
        if nu == 'A':
            seq_comp += 'T'
        elif nu == 'T':
            seq_comp += 'A'
        elif nu == 'C':
            seq_comp += 'G'
        elif nu == 'G':
            seq_comp += 'C'
        else:
            seq_comp += '?'
    
    # reverse the complements sequence
    seq_rev = seq_comp[::-1]
    
    # print all results into a output file
    output = open('./static/reversed_complements.fasta','w')
    print('>sequence',file=output)
    print(seq_rev,file=output)
    output.close()
    
    return(seq_rev)


## function that can find all reversed complements sequence in a file
def comp_rev_content(seq_dic):
    seq_rev = {}
    for key in seq_dic.keys():
        seq_rev[key] = comp_rev(seq_dic[key])
    
    # set the structure of result table
    content=''
    for key in seq_rev.keys():
        content += '{}\n{}\n'.format(key,seq_rev[key])
    
    # print all results into a output file
    output = open('./static/reversed_complements.fasta','w')
    print(content,file=output)
    output.close()
    
    content = content.replace('\n','<br/>')    
    return content


# web part

## function that can redirect to reversed complements input web page
@app.route('/Rev_comp')
def reversed_complements_input():
    return render_template('Reversed_complements_input.html')

## function that can redirect to reversed complements result web page
@app.route('/Rev_comp', methods=['POST'])
def reversed_complements() :
    
    # specifiy the HTTP method to use when sending form-data
    if request.method == 'POST':
        
        # get the submited sequence or file
        seq = request.form.get('seq')
        file = request.files.get('file')     
       
        # if user upload a fasta file
        if file and file_fasta(file.filename): # check the file and file format           
            if seq: # prevent conflicts if users submit sequence and file at the same time
                return render_template('Reversed_complements_input.html') 
            else:  
                seq_dic = sequence_fa(file) # find the id and sequence dictionary
                content = comp_rev_content(seq_dic) # find the reversed complementary sequences
                return render_template('Reversed_complements.html',content=content)
        
        # if user upload a fastq file
        if file and file_fastq(file.filename): # check the file and file format
            if seq: # prevent conflicts if users submit sequence and file at the same time
                return render_template('Reversed_complements_input.html') 
            else:                 
                seq_dic = sequence_fq(file)
                content = comp_rev_content(seq_dic)
                return render_template('Reversed_complements.html',content=content)
        
        # if user input a pure sequence
        if seq :
            if file: # prevent conflicts if users submit sequence and file at the same time
                return render_template('Reversed_complements_input.html') 
            else: 
                content = comp_rev(seq)
                return render_template('Reversed_complements.html',content=content)
   
    return render_template('Reversed_complements_input.html') 


# Translation -----------------------------------------------------------------


# function part 

## dictionary: codon as key and protein as value.
codons_dic = { 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'TGT': 'C', 'TGC': 'C',
                'GAT': 'D', 'GAC': 'D',
                'GAA': 'E', 'GAG': 'E',
                'TTT': 'F', 'TTC': 'F',
                'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
                'CAT': 'H', 'CAC': 'H',
                'ATT': 'I',  'ATA': 'I', 'ATC': 'I',
                'AAA': 'K', 'AAG': 'K',
                'TTA': 'L', 'TTG': 'L',  'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                'ATG': 'M',
                'AAT': 'N', 'AAC': 'N',
                'CCT': 'P', 'CCC': 'P',  'CCA': 'P',  'CCG': 'P',
                'CAA': 'Q',  'CAG': 'Q',
                'CGT': 'R',  'CGC': 'R', 'CGA': 'R',  'CGG': 'R', 'AGA': 'R',   'AGG': 'R',
                'TCT': 'S',  'TCC': 'S',  'TCA': 'S',  'TCG': 'S', 'AGT': 'S',  'AGC': 'S',
                'ACT': 'T', 'ACC': 'T',  'ACA': 'T',  'ACG': 'T',
                'GTT': 'V',  'GTC': 'V',  'GTA': 'V',  'GTG': 'V',
                'TGG': 'W',
                'TAT': 'Y','TAC': 'Y', 
                'TAA': '*', 'TAG': '*',  'TGA': '*'} # stop codon showed as *


## function that can find translation sequence of a nucleotide sequence
def translation(seq):
    aa_seq = ''
    
    # extract codon which has three nucleotides
    for i in range(0,len(seq)-2,3):
        codon = seq[i:(i+3)]
        if len(codon) < 3: # remove last nucleotide which can not be a codon
            break
        
        if codon in codons_dic.keys():
            aa = codons_dic[codon] # find corresponding protein
            aa_seq += aa
        else:
            aa_seq += '-' # if three nucleotides can not be a codon
    
    # print all results into a output file
    output = open('./static/translation.fasta','w')
    print('>sequence',file=output)
    print(aa_seq,file=output)
    output.close()
    
    return(aa_seq)
    

## function that can find all translation sequence in a file
def translation_content(seq_dic):    
    aa_dic = {}
    for key in seq_dic.keys():
        aa_dic[key] = translation(seq_dic[key])    
    content = ''
    for key in aa_dic.keys():
        content += '{}\n{}\n'.format(key,aa_dic[key])
    
    # print all results into a output file
    output = open('./static/translation.fasta','w')
    print(content,file=output)
    output.close()
    
    content = content.replace('\n','<br/>')    
    return content


# web part

## function that can redirect to translation input web page
@app.route('/Translation')
def translation_input():
    return render_template('Translation_input.html')

## function that can redirect to translation results web page
@app.route('/Translation', methods=['POST'])
def translation_output() :
    
    # specifiy the HTTP method to use when sending form-data
    if request.method == 'POST':
        
        seq = request.form.get('seq')
        file = request.files.get('file')     
       
        # if user upload a fasta file
        if file and file_fasta(file.filename): # check the file and file format           
            if seq: # prevent conflicts if users submit sequence and file at the same time
                return render_template('Translation_input.html') 
            else:
                seq_dic = sequence_fa(file)
                content = translation_content(seq_dic) # get the results
                return render_template('Translation.html',content=content)
        
        # if user upload a fastq file
        if file and file_fastq(file.filename): # check the file and file format     
            if seq: # prevent conflicts if users submit sequence and file at the same time
                return render_template('Translation_input.html') 
            else:
                seq_dic = sequence_fq(file)
                content = translation_content(seq_dic)
                return render_template('Translation.html',content=content)
       
        if seq :
            if file: # prevent conflicts if users submit sequence and file at the same time
                return render_template('Translation_input.html') 
            else:
                content = translation(seq)
                return render_template('Translation.html',content=content)
       
    return render_template('Translation_input.html') 


# Gc content figure -----------------------------------------------------------

# import some function that can help to draw the figure
import matplotlib.pyplot as plt
import matplotlib
import io
import base64
matplotlib.use('Agg')

# function that can find the GC content plot
def gc_plot(one_seq,window,step):
    
    window = int(float(window))
    step = int(float(step))
    GC_figure_dic = {}
    
    for i in range(0,(len(one_seq)-window+1),step):
        sub_seq = one_seq[i:(i+window)] # find the sequence in a window
        min_sub_seq = (window-1)/2 + i # find the midpoint of this window
        GC_sub_seq = float(GC_content(sub_seq))/100 # find the GC-content of this window
        GC_figure_dic[min_sub_seq] = GC_sub_seq # add midpoint and GC-content in a dictionry
    
    # read and write UTF-8 encoded picture bytes in memory
    img = io.BytesIO()
    
    # draw the figure
    x = list(GC_figure_dic.keys()) # midpoint of window as x
    y = list(GC_figure_dic.values()) # GC-content value as y
    plt.cla()
    plt.plot(x,y,'b--') # plot the line
    plt.plot(x,y,'bo-')
    plt.title('GC-content') # add the title of figure
    plt.xlabel('position(bp) on sequence') # add the x lable
    plt.ylabel('GC-content') # add the y lable
    plt.savefig(img) # save the img in memory
    img.seek(0) # back to the first frame of the figure
    
    # read figure from bytes data
    plot_url = base64.b64encode(img.getvalue()).decode() # decode the binary data of the image
    return(plot_url)


# function that can redirect to GC plot input web page
@app.route('/GC-figure')
def gc_figure_input():
    return render_template('GC_figure_input.html')

# function that can redirect to GC plot result page
@app.route('/GC-figure', methods=['POST'])
def gc_figure():
    
     # specifiy the HTTP method to use when sending form-data
     if request.method == 'POST':
        
        seq = request.form.get('seq')
        file = request.files.get('file')   
        window = request.form.get('window')
        step = request.form.get('step')
        
        # if user upload a fasta file
        if file and file_fasta(file.filename) and window and step: # check value of window and step not as blank
            if seq: # prevent conflicts if users submit sequence and file at the same time
                return render_template('GC_figure_input.html') 
            else:
                seq_dic = sequence_fa(file)
                one_seq = to_one_seq(seq_dic) # combine all sequences in file as one sequence
                plot_url = gc_plot(one_seq,window,step) # get the GC plot 
                return render_template('GC_figure.html',plot_url=plot_url)
        
        # if user upload a fastq file
        elif file and file_fastq(file.filename) and window and step:
            if seq: # prevent conflicts if users submit sequence and file at the same time
                return render_template('GC_figure_input.html') 
            else:
                seq_dic = sequence_fq(file)
                one_seq = to_one_seq(seq_dic)
                plot_url = gc_plot(one_seq,window,step)
                return render_template('GC_figure.html',plot_url=plot_url)
        
        # if user submit a pure sequence
        elif seq and window and step:
            if file: # prevent conflicts if users submit sequence and file at the same time
                return render_template('GC_figure_input.html') 
            else:
                plot_url = gc_plot(seq,window,step)
                return render_template('GC_figure.html',plot_url=plot_url)
        else:
            return render_template('GC_figure_input.html') 

     return render_template('GC_figure_input.html') 


# Run the flask ---------------------------------------------------------------

if __name__ == '__main__':
    app.run(debug=True)