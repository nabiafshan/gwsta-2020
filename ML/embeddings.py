import pickle
import random
import numpy as np
import textwrap
from dna2vec.multi_k_model import MultiKModel

# with open('modif_cdict.pickle', 'wb') as handle:
#     pickle.dump(modif_cdict, handle)

# with open('modif_ncdict.pickle', 'wb') as handle:
#     pickle.dump(modif_ncdict, handle)

with open('modif_cdict.pickle', 'rb') as handle:
    modif_cdict = pickle.load(handle)

with open('modif_ncdict.pickle', 'rb') as handle:
    modif_ncdict = pickle.load(handle)

# https://stackoverflow.com/a/10125602/11444192

def pickSampFromDict(aDict, size = 1000):
  '''
  Takes as input a dict aDict
  Returns a sample(dict) of given size from values in aDict
  '''
  sampDict = {}
  keys = list(aDict.keys())
  samp_key = random.sample(keys, size)
  for key in samp_key:
    sampDict[key] = aDict[key]
  return sampDict

cSamp = pickSampFromDict(modif_cdict, 3000)
print(len(cSamp))
ncSamp = pickSampFromDict(modif_ncdict, 3000)
print(len(ncSamp))


with open('nt_samp_cdict.pickle', 'wb') as handle:
    pickle.dump(cSamp, handle)

with open('nt_samp_ncdict.pickle', 'wb') as handle:
    pickle.dump(ncSamp, handle)

# with open('nt_samp_cdict.pickle', 'rb') as handle:
#     cSamp = pickle.load(handle)

# with open('nt_samp_ncdict.pickle', 'rb') as handle:
#     ncSamp = pickle.load(handle)

"""

https://blog.keras.io/using-pre-trained-word-embeddings-in-a-keras-model.html

##Preparing text data
"""

text = []
for key in list(cSamp.keys()):
  text.append(str(cSamp[key]))
for key in list(ncSamp.keys()):
  text.append(str(ncSamp[key]))
print(len(text))

with open('nt_6000_c_nc_list.pickle', 'wb') as handle:
    pickle.dump(text, handle)
# with open('nt_6000_c_nc_list.pickle', 'rb') as handle:
#     text = pickle.load(handle)





filepath = 'pretrained/dna2vec-20161219-0153-k3to8-100d-10c-29320Mbp-sliding-Xat.w2v'
mk_model = MultiKModel(filepath)



"""#Get kmers"""

def ntSeqToKmer(aSeq):
  '''
  Given: a nt (ATCG) sequence, 
    get disjoint 8-mers
    get (100, ) embeddings for each 8-mer
  Return: a (100, len(seq)//8) array 
    Each column of array is embedding for a kmer.
  '''
  mer8 = textwrap.wrap(aSeq, 8) #returns list of 8-mer strings
  for i, mer in enumerate(mer8):
    if i == 0:
      a = mk_model.vector(mer).reshape((1,100))
    elif i == len(mer8)-1 and len(mer8[i]) < 3:
      pass
    elif i > 0:
      b = mk_model.vector(mer).reshape((1, 100))
      a = np.concatenate((a,b), axis = 0)
  return a

def ntSeqToKmerDict(aDict):
  '''
  Given: A dict of sequences
  Return: A new dict containing disjoint kmer(8-mer) 
    embeddings for each seq
  '''
  bDict = {}
  for key in list(aDict.keys()):
    seq = str(aDict[key])
    bDict[key] = ntSeqToKmer(seq)
  return bDict

# 6mins 34secs
samp_embed_cdict = ntSeqToKmerDict(cSamp)
samp_embed_ncdict = ntSeqToKmerDict(ncSamp)

with open('samp_embed_cdict.pickle', 'wb') as handle:
    pickle.dump(samp_embed_cdict, handle)

with open('samp_embed_ncdict.pickle', 'wb') as handle:
    pickle.dump(samp_embed_ncdict, handle)

# with open('embed_cdict.pickle', 'rb') as handle:
#     embed_cdict = pickle.load(handle)

# with open('embed_ncdict.pickle', 'rb') as handle:
#     embed_ncdict = pickle.load(handle)
