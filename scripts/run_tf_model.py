import argparse 
import os 
import pathlib 
import pickle 
import numpy as np 
import pandas as pd 
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
plt.ioff() 
from tensorflow.distribute import MirroredStrategy
from tensorflow import keras
from data_objs import TfKmerObj, TfSeqvector
from tensorflow.keras.layers import Embedding, LSTM, Dense, Dropout, Activation, Bidirectional, Input, Masking
from tensorflow.keras.mixed_precision import experimental as mixed_precision
policy = mixed_precision.Policy('mixed_float16')
mixed_precision.set_policy(policy)


parser=argparse.ArgumentParser()
parser.add_argument('--workingDir', action = 'store', dest = 'working_dir')
parser.add_argument('--lsfFile', action = 'store', dest = 'lsf_file')
parser.add_argument('--lsfIDfile', action = 'store', dest = 'lsf_id_file')
parser.add_argument('--labFile', action = 'store', dest='lab_file')
parser.add_argument('--batchSize', action = 'store', type = int, dest='batch_size', default = 100)
parser.add_argument('--nepochs', action = 'store', type = int, default = 1)
parser.add_argument('--nproc', action = 'store', type=int, default=2)
parser.add_argument('--dataType', action = 'store',choices=['kmer', 'seqvector'], dest = 'data_type')
parser.add_argument('--model', action = 'store')
parser.add_argument('--outputDir', action = 'store',  dest='output_dir')
args=parser.parse_args()
os.chdir(args.working_dir)

outdir = args.output_dir
if outdir[-1] is not '/':
    outdir+='/'
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

#no mirror 226ms/step
def base_LSTM(obj):
    mirrored_strategy = MirroredStrategy()
    with mirrored_strategy.scope():
        METRICS = [ 
        keras.metrics.CategoricalAccuracy(name='accuracy'),
        keras.metrics.Precision(name='precision'),
        keras.metrics.Recall(name='recall'),
        keras.metrics.AUC(name='auc')]
        model= keras.models.Sequential(name = 'base_LSTM')
        model.add(Embedding(input_dim = obj.num_words, output_dim = 128,input_length = obj.max_seq_size, mask_zero=True ))
        model.add(Bidirectional(LSTM(units = 128,dropout = .2)))
        model.add(Dense(units = 128, activation = 'relu'))
        model.add(Dropout(.4))
        model.add(Dense(units = 64, activation = 'relu'))
        model.add(Dropout(.2))
        model.add(Dense(units = 21, activation = 'softmax',  dtype='float32'))
        model.compile(optimizer = keras.optimizers.Adam(), loss = 'categorical_crossentropy', metrics = METRICS)
        model.summary()
    return model


def dynamic_LSTM(obj):
    mirrored_strategy = MirroredStrategy()
    with mirrored_strategy.scope():
        METRICS = [ 
        keras.metrics.CategoricalAccuracy(name='accuracy'),
        keras.metrics.Precision(name='precision'),
        keras.metrics.Recall(name='recall'),
        keras.metrics.AUC(name='auc')]
        model= keras.models.Sequential(name = 'dynamic_LSTM')
        model.add(Embedding(input_dim = obj.num_words, output_dim = 128,input_length = None, mask_zero=True ))
        model.add(LSTM(units = 128,dropout = .2))
        model.add(Dense(units = 128, activation = 'relu'))
        model.add(Dropout(.4))
        model.add(Dense(units = 64, activation = 'relu'))
        model.add(Dropout(.2))
        model.add(Dense(units = 21, activation = 'softmax',  dtype='float32'))
        model.compile(optimizer = keras.optimizers.Adam(), loss = 'categorical_crossentropy', metrics = METRICS)
        model.summary()
    return model

def seqvector_LSTM(obj):
    mirrored_strategy = MirroredStrategy()
    with mirrored_strategy.scope():
        METRICS = [ 
        keras.metrics.CategoricalAccuracy(name='accuracy'),
        keras.metrics.Precision(name='precision'),
        keras.metrics.Recall(name='recall'),
        keras.metrics.AUC(name='auc')]
        model = keras.models.Sequential(name = 'seqvector_LSTM')
        model.add(Masking(mask_value = 0., input_shape = obj.train_set[0].shape[1:]))
        model.add(Bidirectional(LSTM(units = 128,dropout = .2)) )
        model.add(Dense(units = 128, activation = 'relu'))
        model.add(Dropout(.4))
        model.add(Dense(units = 64, activation = 'relu'))
        model.add(Dropout(.2))
        model.add(Dense(units = 21, activation = 'softmax',  dtype='float32'))
        model.compile(optimizer = keras.optimizers.Adam(), loss = 'categorical_crossentropy', metrics = METRICS)
        model.summary()
        return model

        

model_dict = {'base_LSTM' : base_LSTM, 'dynamic_LSTM' : dynamic_LSTM, 'LSTM_seqvector': LSTM_seqvector}




with open(outdir+'model_results.csv', 'w+') as model_res_file:
    if args.data_type == 'kmer':
        data_obj = TfKmerObj(args.lsf_file,args.lsf_id_file, args.lab_file, args.batch_size,args.nepochs, outdir)
    elif args.data_type == 'seqvector':
        data_obj = TfSeqvectorObj(args.lsf_file, args.lab_file, args.batch_size,args.nepochs, outdir)
        
    model = model_dict[args.model]
    data_obj.model = model(data_obj)
    data_obj.fit_model()
    
        









