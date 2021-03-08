# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 11:49:01 2021

@author: ZYL
"""

import os
from flask import Flask,request,render_template, redirect, url_for, send_from_directory
from werkzeug.utils import secure_filename


app = Flask(__name__)


@app.route('/')
def Toolbox():
    # return render_template('GC_content_input.html')
    return render_template('Toolbox.html')

@app.route('/functionn_list')
def GC_content_input():
    return render_template('GC_content_input.html')

@app.route('/functionn_list', methods=['POST'])
def show_GC_content() :
    
    if request.method == 'POST':
    
        
        seq = request.form.get('seq')
        file = request.files.get('file')
       
        # if file:
        #     filename = secure_filename(file.filename)
        #     file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        #     file_url = url_for('uploaded_file', filename=filename)
       
        if file :                
            seq_dic = sequence(file)
        
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


if __name__ == '__main__':
    app.run(debug=True)