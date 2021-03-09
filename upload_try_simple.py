# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 11:49:01 2021

@author: ZYL
"""


from flask import Flask,request,render_template
from numpy import mean

app = Flask(__name__)


# preprocessing functions of file ---------------------------------------------


fasta_format = ['fasta','fna','fa']
fastq_format = ['fastq','fq']

def file_fasta(filename):
    return '.' in filename and filename.rsplit('.',1)[1] in fasta_format
 
def file_fastq(filename):
    return '.' in filename and filename.rsplit('.',1)[1] in fastq_format

def sequence(fasta):
    seq_dic = {}
    seq = ''
    seq_id = ''
    for line in fasta:
        line = line.decode('utf-8')
        line = line.strip()
        if line.startswith(">"): # find id line
            if seq: # skip first sequence
                seq_dic[seq_id] = seq        
            seq_id = line # set new sequence id
            seq = ''
        else:
            line = line.upper()
            seq += line # parse sequence in one line
    seq_dic[seq_id] = seq
    return(seq_dic)

def sequence_fq(fastq):
    fq_seq_dic = {}
    fq_seq = ''
    fq_id = ''
    count = 0
    for line in fastq:
        line = line.decode('utf-8')
        line = line.strip()
        count += 1
        if count%4 == 1:
            
            if fq_seq: # skip first sequence
                fq_seq_dic[fq_id] = fq_seq        
            fq_id = line # set new sequence id
            fq_seq = ''
            
        elif count%4 == 2:
            line = line.upper()
            fq_seq = line
            
    fq_seq_dic[fq_id] = fq_seq
    
    return(fq_seq_dic)


# toolbox catalog page --------------------------------------------------------


@app.route('/')
def Toolbox():
    # return render_template('GC_content_input.html')
    return render_template('Toolbox.html')


# GC-content calculation -----------------------------------------------------


# function part 

def GC_content(seq):
    seq = seq.upper().strip()
    length = len(seq)
    count = 0
    for nu in seq:
        if nu == "C":
            count += 1
        if nu == "G":
            count += 1    
    GC_con = count/length*100
    GC_con = '%.2f'%GC_con
    GC_content = str(GC_con)
    return GC_content

def GC123_content(seq):
    GC1_seq = ''
    GC2_seq = ''
    GC3_seq = ''
    for i in range(1,len(seq)+1):
        if i%3 == 1:
            GC1_seq += seq[i-1]
        elif i%3 == 2:
            GC2_seq += seq[i-1]
        else :
            GC3_seq += seq[i-1]
    GC1_content = GC_content(GC1_seq)
    GC2_content = GC_content(GC2_seq)
    GC3_content = GC_content(GC3_seq)
    return [GC1_content, GC2_content, GC3_content]


def GC_content_all(seq_dic):
    GC_dic = {}   
    for key in seq_dic.keys():
        seq = seq_dic[key]
        GC_c = GC_content(seq)
        GC123_c = GC123_content(seq)
        key = key.split()[0][1:] # set new sequence id
        GC_dic[key] = {}
        GC_dic[key]['GC'] = GC_c # add last sequence information into dictionary
        GC_dic[key]['GC1'] = GC123_c[0]
        GC_dic[key]['GC2'] = GC123_c[1]
        GC_dic[key]['GC3'] = GC123_c[2]
            
        content = 'seqID\tGC(%)\tGC1(%)\tGC2(%)\tGC3(%)'
        for key in GC_dic.keys():
            content += '\n{}\t{}\t{}\t{}\t{}'.format(key,GC_dic[key]['GC'],GC_dic[key]['GC1'],GC_dic[key]['GC2'],GC_dic[key]['GC3'])
            
        content = content.replace('\n','<br/>')    
        return content


# web part 

@app.route('/GC-content')
def GC_content_input():
    return render_template('GC_content_input.html')

@app.route('/GC-content', methods=['POST'])
def show_GC_content() :
    
    if request.method == 'POST':
    
        
        seq = request.form.get('seq')
        file = request.files.get('file')     
       
       
        if file and file_fasta(file.filename):                
            seq_dic = sequence(file)
            content = GC_content_all(seq_dic)
            return render_template('GC_content.html',content=content)
        
        if file and file_fastq(file.filename):                
            seq_dic = sequence_fq(file)
            content = GC_content_all(seq_dic)
            return render_template('GC_content.html',content=content)
        
        if seq :
             GC_c = GC_content(seq)
             GC123_c = GC123_content(seq)
             GC_dic = {}  
             GC_dic['GC'] = GC_c # add last sequence information into dictionary
             GC_dic['GC1'] = GC123_c[0]
             GC_dic['GC2'] = GC123_c[1]
             GC_dic['GC3'] = GC123_c[2]
             content = 'GC(%)\tGC1(%)\tGC2(%)\tGC3(%)\n{}\t{}\t{}\t{}'.format(GC_dic['GC'],GC_dic['GC1'],GC_dic['GC2'],GC_dic['GC3'])
             
             content = content.replace('\n','<br/>')
             return render_template('GC_content.html',content=content)
   
    return render_template('GC_content_input.html') 


# length ----------------------------------------------------------------------


# function part

def length_detection(seq_dic):
    
    length_dic = {}
    
    for key in seq_dic.keys():
        GC_id = key.split()[0][1:] # set new sequence id
        length_dic[GC_id] = len(seq_dic[key]) 
    
    length_avg = mean(list(length_dic.values()))
    length_max = max(length_dic.values())
    length_min = min(length_dic.values())
    
    seq_length_max = ''
    seq_length_min = ''
    for key in length_dic.keys():
        if length_dic[key] == length_max:
            if seq_length_max:
                seq_length_max += ','
                seq_length_max += key
            else:
                seq_length_max += key
        if length_dic[key] == length_min:
            if seq_length_min:
                seq_length_min += ','
                seq_length_min += key
            else:
                seq_length_min += key
    
    content = 'character\tlength(pb)\tseqID\naverage_length\t{}\t-\nmaximum_lengths\t{}\t{}\nminimum_lengths\t{}\t{}'.format('%.2f'%length_avg,'%.2f'%length_max,
                                                                                        seq_length_max,'%.2f'%length_min, seq_length_min)
    content = content.replace('\n','<br/>')
    return content


# web part

@app.route('/Length')
def length_input():
    return render_template('Length_detection_input.html')

@app.route('/Length', methods=['POST'])
def length() :
    
    if request.method == 'POST':
        file = request.files.get('file')     
       
        if file and file_fasta(file.filename):                
            seq_dic = sequence(file)
            content = length_detection(seq_dic)
            return render_template('Length_detection.html',content=content)
        
        if file and file_fastq(file.filename):                
            seq_dic = sequence_fq(file)
            content = length_detection(seq_dic)
            return render_template('Length_detection.html',content=content)
   
    return render_template('Length_detection_input.html') 


# reversed complements --------------------------------------------------------

# Step 1: complements
# step 2: reverse 

# function part

def comp_rev(seq):
    seq = seq.upper().strip()
    seq_comp = ''
    for nu in seq:
        if nu == 'A':
            seq_comp += 'T'
        if nu == 'T':
            seq_comp += 'A'
        if nu == 'C':
            seq_comp += 'G'
        if nu == 'G':
            seq_comp += 'C'
    seq_rev = seq_comp[::-1]
    return(seq_rev)

def comp_rev_content(seq_dic):
    seq_rev = {}
    for key in seq_dic.keys():
        seq_rev[key] = comp_rev(seq_dic[key])
    
    content=''
    for key in seq_rev.keys():
        content += '{}\n{}\n'.format(key,seq_rev[key])
        
    content = content.replace('\n','<br/>')    
    return content


# web part

@app.route('/Rev_comp')
def reversed_complements_input():
    return render_template('Reversed_complements_input.html')

@app.route('/Rev_comp', methods=['POST'])
def reversed_complements() :
    
    if request.method == 'POST':
        
        seq = request.form.get('seq')
        file = request.files.get('file')     
       
        if file and file_fasta(file.filename):                
            seq_dic = sequence(file)
            content = comp_rev_content(seq_dic)
    
            return render_template('Reversed_complements.html',content=content)
        
        if file and file_fastq(file.filename):                
            seq_dic = sequence_fq(file)
            content = comp_rev_content(seq_dic)
            return render_template('Reversed_complements.html',content=content)
   
        if seq :
             content = comp_rev(seq)
             return render_template('Reversed_complements.html',content=content)
   
    return render_template('Reversed_complements_input.html') 


# Translation -----------------------------------------------------------------


# function part 

codons_dic = { 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'TGT': 'C', 'TGC': 'C',
                'GAT': 'D', 'GAC': 'D',
                'GAA': 'E', 'GAG': 'E',
                'TTT': 'F', 'TTC': 'F',
                'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
                'CAT': 'H', 'CAC': 'H',
                'ATT': 'I',  'ATA': 'I',
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
                'TAA': '*', 'TAG': '*',  'TGA': '*'}


def translation(seq):
    aa_seq = ''
    for i in range(0,len(seq)-2,3):
        codon = seq[i:(i+3)]
        if len(codon) < 3:
            break
        
        if codon in codons_dic.keys():
            aa = codons_dic[codon]
            aa_seq += aa
        else:
            aa_seq += '-'
    return(aa_seq)
    

def translation_content(seq_dic):    
    aa_dic = {}
    for key in seq_dic.keys():
        aa_dic[key] = translation(seq_dic[key])    
    content = ''
    for key in aa_dic.keys():
        content += '{}\n{}\n'.format(key,aa_dic[key])
    
    content = content.replace('\n','<br/>')    
    return content


# web part

@app.route('/Translation')
def translation_input():
    return render_template('Translation_input.html')

@app.route('/Translation', methods=['POST'])
def translation_output() :
    
    if request.method == 'POST':
        
        seq = request.form.get('seq')
        file = request.files.get('file')     
       
        if file and file_fasta(file.filename):                
            seq_dic = sequence(file)
            content = translation_content(seq_dic)
            return render_template('Translation.html',content=content)
        
        if file and file_fastq(file.filename):                
            seq_dic = sequence_fq(file)
            content = translation_content(seq_dic)
            return render_template('Translation.html',content=content)
   
        if seq :
             content = translation(seq)
             return render_template('Translation.html',content=content)
   
    return render_template('Translation_input.html') 


# Gc content figure -----------------------------------------------------------

import matplotlib.pyplot as plt
import io
import base64

fasta = open('test.fasta','r')

def one_fa(file):
    
    seq = ''
    
    for line in file:
        line = line.decode('utf-8')
        line = line.strip()
        if line.startswith(">"): # find id line
            continue
        else:
            line = line.upper()
            seq += line # parse sequence in one line
    return(seq)


def one_fq(file):

    fq_seq = ''
    
    count = 0
    for line in file:
        line = line.decode('utf-8')
        line = line.strip()
        count += 1
        if count%4 == 2:
            line = line.upper()
            fq_seq += line
    return(fq_seq)


def gc_plot(one_seq,window,step):
    
    window = int(float(window))
    step = int(float(step))
    GC_figure_dic = {}
    
    for i in range(0,(len(one_seq)-window+1),step):
        
        sub_seq = one_seq[i:(i+window)]
       
        min_sub_seq = (window-1)/2 + i
        GC_sub_seq = float(GC_content(sub_seq))/100
     
        GC_figure_dic[min_sub_seq] = GC_sub_seq
    
    img = io.BytesIO()
        
    x = list(GC_figure_dic.keys())
    y = list(GC_figure_dic.values())
    plt.cla()
    plt.plot(x,y,'b--')
    plt.plot(x,y,'bo-')
    plt.title('GC-content')
    plt.xlabel('position(bp) on sequence')
    plt.ylabel('GC-content')
        
    plt.savefig(img,formate='png')
    img.seek(0)
        
    plot_url = base64.b64encode(img.getvalue()).decode()
    return(plot_url)


@app.route('/GC-figure')
def gc_figure_input():
    return render_template('GC_figure_input.html')

@app.route('/GC-figure', methods=['POST'])
def gc_figure():
    
     if request.method == 'POST':
        
        seq = request.form.get('seq')
        file = request.files.get('file')   
        window = request.form.get('window')
        step = request.form.get('step')
        
        if file and file_fasta(file.filename) and window and step:
            if seq:
                return render_template('GC_figure_input.html') 
            else:
                one_seq = one_fa(file)
                plot_url = gc_plot(one_seq,window,step)
                return render_template('GC_figure.html',plot_url=plot_url)
        
        
        elif file and file_fastq(file.filename) and window and step:
            if seq:
                return render_template('GC_figure_input.html') 
            else:
                one_seq = one_fq(file)
                plot_url = gc_plot(one_seq,window,step)
                return render_template('GC_figure.html',plot_url=plot_url)
        
        
        elif seq and window and step:
            if file:
                return render_template('GC_figure_input.html') 
            else:
                plot_url = gc_plot(seq,window,step)
                return render_template('GC_figure.html',plot_url=plot_url)
        else:
            return render_template('GC_figure_input.html') 

     return render_template('GC_figure_input.html') 







if __name__ == '__main__':
    app.run(debug=True)