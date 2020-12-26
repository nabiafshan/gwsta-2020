
import pickle
import numpy as np
from keras.preprocessing.text import Tokenizer
from keras.preprocessing.sequence import pad_sequences
from keras.utils import to_categorical
import itertools

# https://github.com/hzy95/EPIVAN/blob/master/sequence_processing.py
import itertools
import numpy as np
from keras.preprocessing.text import Tokenizer
from keras.preprocessing.sequence import pad_sequences


def sentence2word(str_set):
    word_seq=[]
    for sr in str_set:
        tmp=[]
        for i in range(len(sr)-5):
            if('N' in sr[i:i+6]):
                tmp.append('null')
            else:
                tmp.append(sr[i:i+6])
        word_seq.append(' '.join(tmp))
    return word_seq

def word2num(wordseq,tokenizer,MAX_LEN):
    sequences = tokenizer.texts_to_sequences(wordseq)
    numseq = pad_sequences(sequences, maxlen=MAX_LEN)
    return numseq

def sentence2num(str_set,tokenizer,MAX_LEN):
    wordseq=sentence2word(str_set)
    numseq=word2num(wordseq,tokenizer,MAX_LEN)
    return numseq

def get_tokenizer():
    f= ['a','c','g','t']
    c = itertools.product(f,f,f,f,f,f)
    res=[]
    for i in c:
        temp=i[0]+i[1]+i[2]+i[3]+i[4]+i[5]
        res.append(temp)
    res=np.array(res)
    NB_WORDS = 4097
    tokenizer = Tokenizer(num_words=NB_WORDS)
    tokenizer.fit_on_texts(res)
    acgt_index = tokenizer.word_index
    acgt_index['null']=0
    return tokenizer

def get_data(cRNA, ncRNA):
    tokenizer=get_tokenizer()
    MAX_LEN=2000
    X_c=sentence2num(cRNA,tokenizer,MAX_LEN)
    MAX_LEN=2000
    X_nc=sentence2num(ncRNA,tokenizer,MAX_LEN)
    return X_c, X_nc


# with open('X_c.pickle', 'wb') as handle:
#     pickle.dump(X_c, handle)
# with open('X_nc.pickle', 'wb') as handle:
#     pickle.dump(X_nc, handle)

with open('X_c.pickle', 'rb') as handle:
    X_c = pickle.load(handle)
with open('X_nc.pickle', 'rb') as handle:
    X_nc = pickle.load(handle)

len_cRNA  = 33360
len_ncRNA = 24163
Y = np.concatenate([ np.ones((len_cRNA), dtype=int), np.zeros((len_ncRNA), dtype=int) ])

VALIDATION_SPLIT = 0.2
data = np.vstack((X_c, X_nc))

labels = Y
labels = to_categorical(np.asarray(labels))
print('Shape of data tensor:', data.shape)
print('Shape of label tensor:', labels.shape)

# split the data into a training set and a validation set
indices = np.arange(data.shape[0])
np.random.shuffle(indices)
data = data[indices]
labels = labels[indices]
nb_validation_samples = int(VALIDATION_SPLIT * data.shape[0])

x_train = data[:-nb_validation_samples]
y_train = labels[:-nb_validation_samples]
x_val = data[-nb_validation_samples:]
y_val = labels[-nb_validation_samples:]

from keras.layers import *
from keras.models import *
import numpy as np
from keras.engine.input_layer import Input
from keras.layers import Conv1D, MaxPooling1D, Flatten, Dense
from keras.models import Model
from keras import optimizers
rmsprop = optimizers.RMSprop(lr=0.0001)



MAX_SEQUENCE_LENGTH = 2000
NB_WORDS = 65
EMBEDDING_DIM = 100
embedding_matrix = np.load('embedding_matrix.npy')



# Tokenize-----------------------------
f= ['a','c','g','t']
c = itertools.product(f,f,f,f,f,f)
res=[]
for i in c:
    temp=i[0]+i[1]+i[2]+i[3]+i[4]+i[5]
    res.append(temp)
res=np.array(res)
NB_WORDS = 4097
tokenizer = Tokenizer(num_words=NB_WORDS)
tokenizer.fit_on_texts(res)
word_index = tokenizer.word_index
word_index['null']=0
# ------------------------------------
 

from keras.layers import Embedding

embedding_layer = Embedding(len(word_index),
                            EMBEDDING_DIM,
                            weights=[embedding_matrix],
                            input_length=MAX_SEQUENCE_LENGTH,
                            trainable=True)

sequence_input = Input(shape=(MAX_SEQUENCE_LENGTH, ), dtype='int32')
embedded_sequences = embedding_layer(sequence_input)
x = Conv1D(128, 5, activation='relu')(embedded_sequences)
x = MaxPooling1D(5)(x)
x = Conv1D(128, 5, activation='relu')(x)
x = MaxPooling1D(5)(x)
x = Conv1D(128, 5, activation='relu')(x)
x = MaxPooling1D(35)(x)  # global max pooling
x = Flatten()(x)
x = Dense(128, activation='relu')(x)
preds = Dense(2, activation='softmax')(x)

model = Model(sequence_input, preds)
model.compile(loss='binary_crossentropy',
              optimizer=rmsprop,
              metrics=['acc'])

print(model.summary())

# happy learning!
history = model.fit(x_train, y_train, validation_data=(x_val, y_val), epochs=50, batch_size=512)
