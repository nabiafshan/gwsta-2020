# data from this website
# http://www.rnabinding.com/CPPred/

from Bio import SeqIO
import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import re


# Feature extraction using CPPred paper
def CTD(seq):
	n=len(seq)-1
	n=float(n)
	num_A,num_T,num_G,num_C=0,0,0,0
	AT_trans,AG_trans,AC_trans,TG_trans,TC_trans,GC_trans=0,0,0,0,0,0
	for i in range(len(seq)-1):
		if seq[i]=="A":
			num_A=num_A+1
		if seq[i]=="T":
			num_T=num_T+1
		if seq[i]=="G":
			num_G=num_G+1
		if seq[i]=="C":
			num_C=num_C+1 
		if (seq[i]=="A" and seq[i+1]=="T") or (seq[i]=="T" and seq[i+1]=="A"):
			AT_trans=AT_trans+1
		if (seq[i]=="A" and seq[i+1]=="G") or (seq[i]=="G" and seq[i+1]=="A"):
			AG_trans=AG_trans+1
		if (seq[i]=="A" and seq[i+1]=="C") or (seq[i]=="C" and seq[i+1]=="A"):
			AC_trans=AC_trans+1
		if (seq[i]=="T" and seq[i+1]=="G") or (seq[i]=="G" and seq[i+1]=="T"):
			TG_trans=TG_trans+1
		if (seq[i]=="T" and seq[i+1]=="C") or (seq[i]=="C" and seq[i+1]=="T"):
			TC_trans=TC_trans+1
		if (seq[i]=="G" and seq[i+1]=="C") or (seq[i]=="C" and seq[i+1]=="G"):
			GC_trans=GC_trans+1

	a,t,g,c=0,0,0,0
	A0_dis,A1_dis,A2_dis,A3_dis,A4_dis=0.0,0.0,0.0,0.0,0.0
	T0_dis,T1_dis,T2_dis,T3_dis,T4_dis=0.0,0.0,0.0,0.0,0.0
	G0_dis,G1_dis,G2_dis,G3_dis,G4_dis=0.0,0.0,0.0,0.0,0.0
	C0_dis,C1_dis,C2_dis,C3_dis,C4_dis=0.0,0.0,0.0,0.0,0.0
	for i in range(len(seq)-1):
		if seq[i]=="A":
			a=a+1
			if a == 1:
				A0_dis=((i*1.0)+1)/n
			if a == int(round(num_A/4.0)):
				A1_dis=((i*1.0)+1)/n
			if a == int(round(num_A/2.0)):
				A2_dis=((i*1.0)+1)/n
			if a == int(round((num_A*3/4.0))):
				A3_dis=((i*1.0)+1)/n
			if a == num_A:
				A4_dis=((i*1.0)+1)/n
		if seq[i]=="T":
			t=t+1
			if t == 1:
				T0_dis=((i*1.0)+1)/n
			if t == int(round(num_T/4.0)):
				T1_dis=((i*1.0)+1)/n
			if t == int(round((num_T/2.0))):
				T2_dis=((i*1.0)+1)/n
			if t == int(round((num_T*3/4.0))):
				T3_dis=((i*1.0)+1)/n
			if t == num_T:
				T4_dis=((i*1.0)+1)/n
		if seq[i]=="G":
			g=g+1
			if g == 1:
				G0_dis=((i*1.0)+1)/n
			if g == int(round(num_G/4.0)):
				G1_dis=((i*1.0)+1)/n
			if g == int(round(num_G/2.0)):
				G2_dis=((i*1.0)+1)/n
			if g == int(round(num_G*3/4.0)):
				G3_dis=((i*1.0)+1)/n
			if g == num_G:
				G4_dis=((i*1.0)+1)/n
		if seq[i]=="C":
			c=c+1
			if c == 1:
				C0_dis=((i*1.0)+1)/n
			if c == int(round(num_C/4.0)):
				C1_dis=((i*1.0)+1)/n
			if c == int(round(num_C/2.0)):
				C2_dis=((i*1.0)+1)/n
			if c == int(round(num_C*3/4.0)):
				C3_dis=((i*1.0)+1)/n
			if c == num_C:
				C4_dis=((i*1.0)+1)/n
	return(str(num_A/n),str(num_T/n),str(num_G/n),str(num_C/n),str(AT_trans/(n-1)),str(AG_trans/(n-1)),str(AC_trans/(n-1)),str(TG_trans/(n-1)),str(TC_trans/(n-1)),str(GC_trans/(n-1)),str(A0_dis),str(A1_dis),str(A2_dis),str(A3_dis),str(A4_dis),str(T0_dis),str(T1_dis),str(T2_dis),str(T3_dis),str(T4_dis),str(G0_dis),str(G1_dis),str(G2_dis),str(G3_dis),str(G4_dis),str(C0_dis),str(C1_dis),str(C2_dis),str(C3_dis),str(C4_dis))

"""##ProtParams"""

'''
Extract the most probable ORF in a given sequence 
The most probable ORF is the longest open reading frame found in the sequence
When having same length, the upstream ORF is selected
modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/v1.2.1/
'''

class ExtractORF:
	def __init__(self,seq):
		self.seq = seq
		self.result=(0,0,0,0)
		self.longest = 0
	
	def codons(self,frame):
		start_coord = frame
		while start_coord + 3 <= len(self.seq):
			yield (self.seq[start_coord:start_coord+3],start_coord)
			start_coord += 3
	def longest_orf_in_seq(self,frame_number,start_codon,stop_codon):
		codon_posi = self.codons(frame_number)
		start_codons = start_codon
		stop_codons = stop_codon
		while True:
			try:
				codon,index = next(codon_posi)
			except StopIteration:
				break
			if codon in start_codons and codon not in stop_codons:
				ORF_start = index
				end = False
				while True:
					try:
						codon,index = next(codon_posi)
					except StopIteration:
						end = True
						integrity = -1
					if codon in stop_codons:
						integrity = 1
						end = True
					if end:
						ORF_end = index+3
						ORF_Length=(ORF_end-ORF_start)
						if ORF_Length > self.longest:
							self.longest = ORF_Length
							self.result = (integrity,ORF_start,ORF_end,ORF_Length)
						if ORF_Length == self.longest and ORF_start < self.result[1]:
							self.result = (integrity,ORF_start,ORF_end,ORF_Length)
						break
	def longest_ORF(self,start=['ATG'],stop=['TAA','TAG','TGA']):
		orf_seq = ""
		for frame in range(3):
			self.longest_orf_in_seq(frame,start,stop)
		orf_seq = self.seq[self.result[1]:self.result[2]]
		ORF_integrity = self.result[0]
		ORF_length = self.result[3]
		return ORF_length,ORF_integrity,orf_seq

import sys
import re
from Bio.Seq import Seq
# from ORF import ExtractORF
from Bio.SeqUtils import ProtParam

def mRNA_translate(mRNA):
	return Seq(mRNA).translate()

def protein_param(putative_seqprot):
	return (putative_seqprot.instability_index(),putative_seqprot.isoelectric_point(),putative_seqprot.gravy())

def param(seq):
	strinfoAmbiguous = re.compile("X|B|Z|J|U",re.I)
	ptU = re.compile("U",re.I)
	seqRNA = ptU.sub("T",str(seq).strip())
	seqRNA = seqRNA.upper()
	CDS_size1,CDS_integrity,seqCDS= ExtractORF(seqRNA).longest_ORF(start=['ATG'],stop=['TAA','TAG','TGA'])
	seqprot = mRNA_translate(seqCDS)
	pep_len = len(seqprot.strip("*"))
	newseqprot = strinfoAmbiguous.sub("",str(seqprot))
	protparam_obj = ProtParam.ProteinAnalysis(str(newseqprot.strip("*")))
	if pep_len > 0:
		Instability_index,PI,Gravy = protein_param(protparam_obj)
	else:
		Instability_index = 0.0
		PI=0.0
		Gravy=0.0
	return(Instability_index,PI,Gravy)

"""##Fickett"""
'''the python script is downloaded from https://sourceforge.net/projects/rna-cpat/files/v1.2.2/'''
'''calculate coding potential'''

# Fickett TESTCODE data
# NAR 10(17) 5303-531
position_prob ={
'A':[0.94,0.68,0.84,0.93,0.58,0.68,0.45,0.34,0.20,0.22],
'C':[0.80,0.70,0.70,0.81,0.66,0.48,0.51,0.33,0.30,0.23],
'G':[0.90,0.88,0.74,0.64,0.53,0.48,0.27,0.16,0.08,0.08],
'T':[0.97,0.97,0.91,0.68,0.69,0.44,0.54,0.20,0.09,0.09]
}
position_weight={'A':0.26,'C':0.18,'G':0.31,'T':0.33}
position_para  =[1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,0.0]

content_prob={
'A':[0.28,0.49,0.44,0.55,0.62,0.49,0.67,0.65,0.81,0.21],
'C':[0.82,0.64,0.51,0.64,0.59,0.59,0.43,0.44,0.39,0.31],
'G':[0.40,0.54,0.47,0.64,0.64,0.73,0.41,0.41,0.33,0.29],
'T':[0.28,0.24,0.39,0.40,0.55,0.75,0.56,0.69,0.51,0.58]
}
content_weight={'A':0.11,'C':0.12,'G':0.15,'T':0.14}
content_para  =[0.33,0.31,0.29,0.27,0.25,0.23,0.21,0.17,0]

def look_up_position_prob(value, base):
	'''look up positional probability by base and value'''
	if float(value)<0:
		return None
	for idx,val in enumerate (position_para):
		if (float(value) >= val):
			return float(position_prob[base][idx]) * float(position_weight[base])

def look_up_content_prob(value, base):
	'''look up content probability by base and value'''
	if float(value)<0:
		return None
	for idx,val in enumerate (content_para):
		if (float(value) >= val):
			return float(content_prob[base][idx]) * float(content_weight[base])

def fickett_value(dna):
	'''calculate Fickett value. Input is DNA sequence'''
	if len(dna)<2:
		return 0
	fickett_score=0
	dna=dna.upper()
	total_base = len(dna)
	A_content = float(dna.count('A'))/total_base
	C_content = float(dna.count('C'))/total_base
	G_content = float(dna.count('G'))/total_base
	T_content = float(dna.count('T'))/total_base
	#print "A content\t" + str(A_content)
	#print "C content\t" + str(C_content)
	#print "G content\t" + str(G_content)
	#print "T content\t" + str(T_content)
	
	phase_0 = [dna[i] for i in range(0,len(dna)) if i % 3==0]
	phase_1 = [dna[i] for i in range(0,len(dna)) if i % 3==1]
	phase_2 = [dna[i] for i in range(0,len(dna)) if i % 3==2]
	
	A_position=max(phase_0.count('A'),phase_1.count('A'),phase_2.count('A'))/(min(phase_0.count('A'),phase_1.count('A'),phase_2.count('A')) +1.0)
	C_position=max(phase_0.count('C'),phase_1.count('C'),phase_2.count('C'))/(min(phase_0.count('C'),phase_1.count('C'),phase_2.count('C')) +1.0)
	G_position=max(phase_0.count('G'),phase_1.count('G'),phase_2.count('G'))/(min(phase_0.count('G'),phase_1.count('G'),phase_2.count('G')) +1.0)
	T_position=max(phase_0.count('T'),phase_1.count('T'),phase_2.count('T'))/(min(phase_0.count('T'),phase_1.count('T'),phase_2.count('T')) +1.0)
	#print "A position\t" + str(A_position)
	#print "C position\t" + str(C_position)
	#print "G position\t" + str(G_position)
	#print "T position\t" + str(T_position)

	
	#for i (A_content,C_content,G_content,T_content):
	fickett_score += look_up_content_prob(A_content,'A')
	fickett_score += look_up_content_prob(C_content,'C')
	fickett_score += look_up_content_prob(G_content,'G')
	fickett_score += look_up_content_prob(T_content,'T')
	
	fickett_score += look_up_position_prob(A_position,'A')
	fickett_score += look_up_position_prob(C_position,'C')
	fickett_score += look_up_position_prob(G_position,'G')
	fickett_score += look_up_position_prob(T_position,'T')

	return fickett_score

"""##FrameKmer"""

'''the python script is downloaded from https://sourceforge.net/projects/rna-cpat/files/v1.2.2/'''
'''deal with Kmer. DNA sequence should only A, C, G, T. python2.7 or newer'''

import os,sys
import math

def word_generator(seq,word_size,step_size,frame=0):
	'''generate DNA word from sequence using word_size and step_size. Frame is 0, 1 or2'''
	for i in range(frame,len(seq),step_size):
		word =  seq[i:i+word_size]
		if len(word) == word_size:
			yield word

def kmer_ratio(seq,word_size,step_size,coding,noncoding):
	if len(seq) < word_size:
		return 0
		
	sum_of_log_ratio_0 = 0.0
	sum_of_log_ratio_1 = 0.0
	sum_of_log_ratio_2 = 0.0	
	frame0_count=0.0
	frame1_count=0.0
	frame2_count=0.0
	for k in word_generator(seq=seq, word_size = word_size, step_size=step_size,frame=0):	
		if (not k in coding) or (not k in noncoding):
			continue
		if coding[k]>0 and noncoding[k] >0:
			sum_of_log_ratio_0  +=  math.log( coding[k] / noncoding[k])
		elif coding[k]>0 and noncoding[k] == 0:
			sum_of_log_ratio_0 += 1
		elif coding[k] == 0 and noncoding[k] == 0:
			continue
		elif coding[k] == 0 and noncoding[k] >0 :
			sum_of_log_ratio_0 -= 1
		else:
			continue
		frame0_count += 1
	'''	
	for k in word_generator(seq=seq, word_size = word_size, step_size=step_size,frame=1):
		if (not coding.has_key(k)) or (not noncoding.has_key(k)):
			continue
		if coding[k]>0 and noncoding[k] >0:
			sum_of_log_ratio_1  +=  math.log( coding[k] / noncoding[k])
		elif coding[k]>0 and noncoding[k] == 0:
			sum_of_log_ratio_1 += 1
		elif coding[k] == 0 and noncoding[k] == 0:
			continue
		elif coding[k] == 0 and noncoding[k] >0 :
			sum_of_log_ratio_1 -= 1
		else:
			continue
		frame1_count += 1
	
	for k in word_generator(seq=seq, word_size = word_size, step_size=step_size,frame=2):
		if (not coding.has_key(k)) or (not noncoding.has_key(k)):
			continue
		if coding[k]>0 and noncoding[k] >0:
			sum_of_log_ratio_2  +=  math.log( coding[k] / noncoding[k])
		elif coding[k]>0 and noncoding[k] == 0:
			sum_of_log_ratio_2 += 1
		elif coding[k] == 0 and noncoding[k] == 0:
			continue
		elif coding[k] == 0 and noncoding[k] >0 :
			sum_of_log_ratio_2 -= 1
		else:
			continue
		frame2_count += 1
	return max(sum_of_log_ratio_0/frame0_count, sum_of_log_ratio_1/frame1_count,sum_of_log_ratio_2/frame2_count)	
	'''
	try:
		return sum_of_log_ratio_0/frame0_count
	except:
		return -1

"""##ORF Length"""

import string
import sys
# from string import maketrans
# import ORF

def extract_feature_from_seq(seq,stt,stp):
	'''extract features of sequence from fasta entry'''
	
	stt_coden = stt.strip().split(',')
	stp_coden = stp.strip().split(',')
	transtab = str.maketrans("ACGTNX","TGCANX")
	mRNA_seq = seq.upper()
	mRNA_size = len(seq)
	tmp = ExtractORF(mRNA_seq)
	(CDS_size1, CDS_integrity, CDS_seq1) = tmp.longest_ORF(start=stt_coden, stop=stp_coden)
	return (mRNA_size, CDS_size1,CDS_integrity)

start_codons = 'ATG'
stop_codons = 'TAG,TAA,TGA'
Coverage = 0

def len_cov(seq):
	(mRNA_size, CDS_size,CDS_integrity) = extract_feature_from_seq(seq = seq, stt = start_codons,stp = stop_codons)
	mRNA_len = mRNA_size
	CDS_len = CDS_size
	Coverage = float(CDS_len)/mRNA_len
	Integrity = CDS_integrity
	return(CDS_len,Coverage,Integrity)

"""##Hexamer"""

def coding_nocoding_potential(input_file):
	coding={}
	noncoding={}
	for line in open(input_file).readlines():
		fields = line.split()
		if fields[0] == 'hexamer':continue
		coding[fields[0]] = float(fields[1])
		noncoding[fields[0]] =  float(fields[2])
	return coding,noncoding

"""
## Run CPPred functions
"""

def output_feature(seq_file,hex_file,species):
	tmp = open('test.f_svm','w')
	feature = open('test.feature','w')
	out_label = 1
	coding,noncoding = coding_nocoding_potential(hex_file)
	if species == "Human":
		feature.write("\t".join(map(str,["#ID","ORF-integrity","ORF-coverage","Instability","T2","C0","PI","ORF-length","AC","T0","G0","C2","A4","G2","TG","A0","TC","G1","C3","T3","A1","GC","T1","G4","C1","G3","A3","Gravy","Hexamer","C4","AG","Fickett","A2","T4","C","G","A","T"]))+"\n")
	if species == "Integrated":
		feature.write("\t".join(map(str,["#ID","ORF-coverage","ORF-integrity","GC","Instability","ORF-length","T0","Fickett","G2","C3","PI","A3","C1","G3","Hexamer","TG","G1","TC","A0","A1","AC","C2","G0","T4","C0","A4","G","A2","T","T3","G4","C4","Grary","T2","AG","AT","T1","A","C"]))+"\n")
	for seq in Seq.parse(seq_file,'fasta'):
		seqid = seq.id
		A,T,G,C,AT,AG,AC,TG,TC,GC,A0,A1,A2,A3,A4,T0,T1,T2,T3,T4,G0,G1,G2,G3,G4,C0,C1,C2,C3,C4 = CTD.CTD(seq.seq)
		insta_fe,PI_fe,gra_fe = PP.param(seq.seq)
		fickett_fe = fickett.fickett_value(seq.seq)
		hexamer = FrameKmer.kmer_ratio(seq.seq,6,3,coding,noncoding)
		Len,Cov,inte_fe = len.len_cov(seq.seq)
		if species == "Human":
			tem = [inte_fe,Cov,insta_fe,T2,C0,PI_fe,Len,AC,T0,G0,C2,A4,G2,TG,A0,TC,G1,C3,T3,A1,GC,T1,G4,C1,G3,A3,gra_fe,hexamer,C4,AG,fickett_fe,A2,T4,C,G,A,T]
			feature.write("\t".join(map(str,[seqid,inte_fe,Cov,insta_fe,T2,C0,PI_fe,Len,AC,T0,G0,C2,A4,G2,TG,A0,TC,G1,C3,T3,A1,GC,T1,G4,C1,G3,A3,gra_fe,hexamer,C4,AG,fickett_fe,A2,T4,C,G,A,T]))+"\n")
		if species == "Integrated":
			tem = [Cov,inte_fe,GC,insta_fe,Len,T0,fickett_fe,G2,C3,PI_fe,A3,C1,G3,hexamer,TG,G1,TC,A0,A1,AC,C2,G0,T4,C0,A4,G,A2,T,T3,G4,C4,gra_fe,T2,AG,AT,T1,A,C]
			feature.write("\t".join(map(str,[seqid,Cov,inte_fe,GC,insta_fe,Len,T0,fickett_fe,G2,C3,PI_fe,A3,C1,G3,hexamer,TG,G1,TC,A0,A1,AC,C2,G0,T4,C0,A4,G,A2,T,T3,G4,C4,gra_fe,T2,AG,AT,T1,A,C]))+"\n")
		print >> tmp, out_label,
		for label,item in enumerate(tem):
			print >> tmp, str(label+1)+':'+str(item),
		print >> tmp
	tmp.close()

seq_file = 'Human.coding_RNA_training.fa'
hex_file = 'Human_Hexamer.tsv'
coding,noncoding = coding_nocoding_potential(hex_file)
i = 0
cDict = {}
parser = SeqIO.parse(open(seq_file), 'fasta')
for seq in parser:
  seqid = seq.id
  seq = seq.seq  
  hexamer = kmer_ratio(seq,6,3,coding,noncoding)
  # ------------------------------------------------------------------------------
  A,T,G,C,AT,AG,AC,TG,TC,GC,A0,A1,A2,A3,A4,T0,T1,T2,T3,T4,G0,G1,G2,G3,G4,C0,C1,C2,C3,C4 = CTD(seq)
  insta_fe,PI_fe,gra_fe = param(seq)
  fickett_fe = fickett_value(seq)
  Len,Cov,inte_fe = len_cov(seq)
  tem = [Cov,inte_fe,GC,insta_fe,Len,T0,fickett_fe,G2,C3,
         PI_fe,A3,C1,G3,hexamer,TG,G1,TC,A0,A1,AC,C2,G0,T4,
         C0,A4,G,A2,T,T3,G4,C4,gra_fe,T2,AG,AT,T1,A,C]
  tem = [float(i) for i in tem]
  cDict[seqid] = tem

"""---
# Read Pickles
"""

# with open('cDict.pickle', 'wb') as handle:
#     pickle.dump(cDict, handle)

# with open('ncDict.pickle', 'wb') as handle:
#     pickle.dump(ncDict, handle)

with open('cDict.pickle', 'rb') as handle:
    cDict = pickle.load(handle)

with open('ncDict.pickle', 'rb') as handle:
    ncDict = pickle.load(handle)

# Number of coding and non-coding samples
print(len(cDict), len(ncDict))

cdf = pd.DataFrame.from_dict(cDict, orient='index')
cdf.columns = ["#ID","ORF-integrity","ORF-coverage","Instability","T2","C0","PI","ORF-length","AC","T0","G0","C2","A4","G2","TG","A0","TC","G1","C3","T3","A1","GC","T1","G4","C1","G3","A3","Gravy","Hexamer","C4","AG","Fickett","A2","T4","C","G","A","T"]
cdf['label'] = 1
ncdf = pd.DataFrame.from_dict(ncDict, orient='index')
ncdf.columns = ["#ID","ORF-integrity","ORF-coverage","Instability","T2","C0","PI","ORF-length","AC","T0","G0","C2","A4","G2","TG","A0","TC","G1","C3","T3","A1","GC","T1","G4","C1","G3","A3","Gravy","Hexamer","C4","AG","Fickett","A2","T4","C","G","A","T"]
ncdf['label'] = 0

mydf = pd.concat([cdf, ncdf])
# print(cdf.shape, ncdf.shape, mydf.shape)

"""#SVM"""

from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC

Y = mydf['label']
X = mydf.drop(['label'], axis=1)
Xtrain, Xtest, Ytrain, Ytest = train_test_split(X, Y, test_size=0.2)

clf = SVC(gamma='scale')
clf.fit(Xtrain, Ytrain) 

Ypred_train = clf.predict(Xtrain)
print(confusion_matrix(Ytrain, Ypred_train))
print('')

print('Test Classification report:')
print(classification_report(Ytrain, Ypred_train))
print('')

print('-------------------------------------------------------')

Ypred = clf.predict(Xtest)
print('Test Confusion Matrix: ')
print(confusion_matrix(Ytest, Ypred))
print('')

print('Test Classification report:')
print(classification_report(Ytest, Ypred))

clf = SVC()
clf.fit(Xtrain, Ytrain) 
Ypred = clf.predict(Xtest)
print('Confusion Matrix: ')
print(confusion_matrix(Ytest, Ypred))
print('')

print('classification report:')
print(classification_report(Ytest, Ypred))

"""#Confusion Matrix Plot"""

# https://scikit-learn.org/stable/auto_examples/model_selection/plot_confusion_matrix.html#sphx-glr-auto-examples-model-selection-plot-confusion-matrix-py
def plot_confusion_matrix(y_true, y_pred, classes,
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    classes = classes[unique_labels(y_true, y_pred)]
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    return ax

# Plot non-normalized confusion matrix
plot_confusion_matrix(y_test, y_pred, classes=class_names,
                      title='Confusion matrix, without normalization')

# Plot normalized confusion matrix
plot_confusion_matrix(y_test, y_pred, classes=class_names, normalize=True,
                      title='Normalized confusion matrix')

plt.show()

# ROC Curve 
# https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html#sphx-glr-auto-examples-model-selection-plot-roc-py
  
# ROC Curve With Crossvalidation
# https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html#sphx-glr-auto-examples-model-selection-plot-roc-crossval-py

"""#Grid Search"""

# Commented out IPython magic to ensure Python compatibility.
# Set the parameters by cross-validation
tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-2, 1e-3, 1e-4, 1e-5],
                     'C': [0.001, 0.10, 0.1, 10, 25, 50, 100, 1000]},
                    {'kernel': ['sigmoid'], 'gamma': [1e-2, 1e-3, 1e-4, 1e-5],
                     'C': [0.001, 0.10, 0.1, 10, 25, 50, 100, 1000]},
                    {'kernel': ['linear'], 'C': [0.001, 0.10, 0.1, 10, 25, 50, 100, 1000]}
                   ]

scores = ['f1']

print("# Tuning hyper-parameters for %s" % scores)
print()

clf = GridSearchCV(SVC(C=1), tuned_parameters, cv=5,
                   scoring='f1', return_train_score=True)
clf.fit(Xtrain, Ytrain)

print("Best parameters set found on development set:")
print()
print(clf.best_params_)
print()
print("Grid scores on development set:")
print()
means = clf.cv_results_['mean_test_score']
stds = clf.cv_results_['std_test_score']
for mean, std, params in zip(means, stds, clf.cv_results_['params']):
    print("%0.3f (+/-%0.03f) for %r"
#           % (mean, std * 2, params))
