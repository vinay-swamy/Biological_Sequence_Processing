#%%
import numpy as np 
import pandas as pd 
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
plt.ioff() 
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, auc
from sklearn.model_selection import train_test_split
from tensorflow import keras
import pathlib
import pickle 
from Bio import SeqIO

    
class TfDataObj:
    ''' container for data formatted the way tensorflow likes'''
    def __init__(self,tokenized_data, lsf_ids,  label_file, batch_size, nepochs, outdir):
        
        self.max_seq_size = max([len(i) for i in tokenized_data])
        lab_df = pd.read_csv(label_file).rename({'family_id' : 'target_label'},axis =1)
        lab_df = lab_df[~lab_df.duplicated('embl_id')].set_index('embl_id')
        family_ids = [int(lab_df.loc[i]['target_label']) for i in lsf_ids]
        targ_labs = keras.utils.to_categorical(family_ids)
        tokenized_data=np.asarray(tokenized_data)
        X_train, X_val, Y_train, Y_val = train_test_split(tokenized_data, targ_labs, test_size = .2, stratify = family_ids, random_state=32)
        X_val, X_test, Y_val, Y_test = train_test_split(X_val, Y_val, test_size = .5, stratify = [np.argmax(i) for i in Y_val ], random_state = 42)
        X_val = np.array([np.array(i) for i in X_val])
        X_test = np.array([np.array(i) for i in X_test])
        # even_batch_size = int(batch_size * int(len(X_test)/ batch_size))# even out batch_sizes
        # X_test=X_test[:even_batch_size]
        # Y_test=Y_test[:even_batch_size]
        
        self.train_set = (X_train, Y_train)
        self.val_set = (X_val, Y_val)
        self.test_set = (X_test,Y_test)
        self.batch_size = batch_size
        self.nepochs=nepochs
        self.trained_model=None
        self.outdir=outdir
        self.model = None
        self.history=None
    
    def plot_metrics(self):
        history = self.history
        outdir = self.outdir
        metrics =  ['loss', 'auc', 'accuracy']
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for n, metric in enumerate(metrics):
            name = metric.replace("_"," ").capitalize()
            plt.subplot(2,2,n+1)
            plt.plot(history.epoch,  history.history[metric], color=colors[0], label='Train')
            plt.plot(history.epoch, history.history['val_'+metric],
                     color=colors[0], linestyle="--", label='Val')
            plt.xlabel('Epoch')
            plt.ylabel(name)
            if metric == 'loss':
                plt.ylim([0, plt.ylim()[1]])
            elif metric == 'auc':
                plt.ylim([0.7,1])
            else:
                plt.ylim([.5,1.1])

            plt.legend()
        plt.suptitle('training metrics')
        plt.savefig(outdir + f'training_metrics.png')
        plt.close()

    def fit_model(self):
        '''pad all tokens upto the max length in the data set'''
        history = self.model.fit(x = self.train_set[0], y = self.train_set[1] ,
                            epochs=self.nepochs, batch_size = self.batch_size, 
                            validation_data=self.val_set,
                            verbose = 1)
        self.history = history
        self.model.save(self.outdir + 'trained_model.h5')
        self.plot_metrics()


    def model_results(self, Y_true, Y_pred_class, Y_pred_prob):
        data_name  = f'{self.name}_{self.model.name}'
        Y_true=np.asarray(Y_true)
        Y_pred_class = np.asarray(Y_pred_class)
        print(Y_pred_prob.shape)
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        mpl.rcParams['figure.figsize'] = (12, 10)


        print(classification_report(y_pred=Y_pred_class, y_true=Y_true))
        cr_dict = classification_report(
            y_pred=Y_pred_class, y_true=Y_true, output_dict=True)
    def predict_model(self):
        pred_probs=self.model.predict(self.test_set[0])
        pred_class = np.asarray([np.argmax(p) for p in pred_probs])
        true_class = np.asarray([np.argmax(p) for p in self.test_set[1]])
        self.model_results( true_class, pred_class, pred_probs)
    
                
##### class for dynamic RNN, unfinished. 
# class TfGeneratorObj()            
#     def fit_model(self):
#         '''make variable length batches for trainging using a generator '''
        
        
#         batch_size = self.batch_size
#         nepochs = self.nepochs
#         training_generator = make_padded_batches(self.train_set[0], self.train_set[1],batch_size, nepochs)
#         validation_generator = make_padded_validation_batch(self.val_set[0], self.val_set[1], batch_size)
#         testing_generator = make_padded_batches(self.test_set[0], self.test_set[1],3, 1)
#         trainsteps = int(len(self.train_set[0]) / batch_size)
#         valsteps = int(len(self.val_set[0]) / batch_size)
#         history = self.model.fit(x = training_generator, 
#                             epochs=nepochs, steps_per_epoch = trainsteps,  
#                             validation_data=validation_generator, validation_steps=valsteps,  
#                             verbose = 1)
#         self.model.save(self.outdir + 'trained_model.h5')
#         plot_metrics(history, self.outdir)
#         pred_probs= self.model.predict(testing_generator)
#         pred_class = np.asarray([np.argmax(p) for p in pred_probs])
        
#         true_class = np.asarray([np.argmax(p) for p in self.test_set[1]])
#         model_results( true_class, pred_class, pred_probs, f'{self.name}_{self.model.name}', self.outdir)
    
#     def load_model(self, path_to_model):
#         self.model = keras.models.load_model(path_to_model)
#         return model_res_line
#     def zip_sort_unzip(X, Y):
#         XY = list(zip(X, Y))
#         XY.sort(key=lambda x: len(x[0]))
#         X = [i[0] for i in XY]
#         Y = np.array([i[1] for i in XY])
#         return X, Y



#     def make_padded_batches(X, Y, batch_size, number_of_epochs):
#         #keep batch sizes to an even size, because batch_sizes<=3 can some times cause code to fail
#         X,Y = zip_sort_unzip(X, Y)
#         even_batch_size = int(batch_size * int(len(X)/ batch_size))# even out batch_sizes
#         X=X[:even_batch_size]
#         Y=Y[:even_batch_size]
#         for _ in range(number_of_epochs):
#             for i in range(0, len(X), batch_size):
#                 c_X = X[i:i+batch_size]
#                 c_Y = Y[i:i+batch_size]
#                 max_len=max([len(i) for i in c_X]) 
#                 c_X_msk = keras.preprocessing.sequence.pad_sequences(c_X, max_len, padding = 'post')
#                 #dont run batches less than batch size
#                 yield c_X_msk, c_Y
#
#     def make_padded_validation_batch(X, Y, batch_size):
#         '''
#         I cant exactly figure out why but the model consumes more batches than what its supposed to, so just going to keep giving it the data set until it finishes 
#         '''
#         X,Y = zip_sort_unzip(X,Y)
#         even_batch_size = int(batch_size * int(len(X)/ batch_size))# even out batch_sizes
#         X=X[:even_batch_size]
#         Y=Y[:even_batch_size]
#         while True:# infinitelt produces batches to 
#             for i in range(0, len(X), batch_size):
#                 c_X = X[i:i+batch_size]
#                 c_Y = Y[i:i+batch_size]
#                 max_len=max([len(i) for i in c_X]) 
#                 c_X_msk = keras.preprocessing.sequence.pad_sequences(c_X, max_len, padding = 'post')
#                 yield c_X_msk, c_Y


class TfKmerObj(TfDataObj):
    def __init__(self,lsf_file, lsf_id_file,lable_file, batch_size, nepochs, outdir, token_dir = 'data/tokenized_data/' ):
        fn_pfx= lsf_file.split('/')[-1].split('.')[0]
        tokenized_lsf_file = token_dir + fn_pfx + '_data.pickle'
        tokenizer_file = token_dir + fn_pfx +'_tokenizer.pickle'
        pickled_lsf_ids = token_dir + fn_pfx + '_txids.pickle'
        if pathlib.Path(tokenized_lsf_file).is_file():# only tokenize data once 
            print('loading tokenized data\n') 
            with open(tokenized_lsf_file, 'rb') as tlf, open(tokenizer_file, 'rb') as tzr, open(pickled_lsf_ids, 'rb') as ptx:
                tokenized_data = pickle.load(tlf)
                tokenizer = pickle.load(tzr)
                lsf_ids = pickle.load(ptx)
            
        else:
            print('tokenizing and saving data')
            tokenized_data, tokenizer = self.tokenize_keras(lsf_file)
            ### need to remove any transcripts that are completely out of vocab https://github.com/tensorflow/tensorflow/issues/33148

            tokenized_data = np.array([np.array(i) for i in tokenized_data])
            tsums = np.array([sum(i) for i in tokenized_data])
            tokenized_data=tokenized_data[tsums!=0]    
            
            with open(lsf_id_file) as f :
                lsf_ids= np.asarray([line.strip('\n') for line in f ])
            lsf_ids = lsf_ids[tsums!=0]
            ####
        
            with open(tokenized_lsf_file, 'wb+') as tlf, open(tokenizer_file, 'wb+') as tzr, open(pickled_lsf_ids, 'wb+') as ptx:
                pickle.dump(tokenized_data, tlf)
                pickle.dump(tokenizer, tzr)
                pickle.dump(lsf_ids, ptx)
                print('tokenization Finished\n')
        ##finished loading
        self.name = f'{fn_pfx}_bs-{batch_size}_ne-{nepochs}'
        self.num_words = len(tokenizer.word_index) +1 #total number of words  + the oov token 
        print(self.num_words)
        TfDataObj.__init__(self,tokenized_data, lsf_ids, lable_file,  batch_size, nepochs, outdir)
        # pad sequences 
        self.train_set[0] = keras.preprocessing.sequence.pad_sequences(self.train_set[0], self.max_seq_size, padding = 'post')
        self.val_set[0] = keras.preprocessing.sequence.pad_sequences(self.val_set[0], self.max_seq_size, padding = 'post')
        self.test_set[0] = keras.preprocessing.sequence.pad_sequences(self.test_set[0], self.max_seq_size, padding = 'post')
    def tokenize_keras(self, file):
        ''' 
        convert lsf based kmers into tokenized kmers 
        remove rare kmers to improve generalizablity
        '''
        with open(file) as  infile:
            sentences = [line.strip() for line in infile ]
        tokenizer=keras.preprocessing.text.Tokenizer(oov_token='MSG')
        tokenizer.fit_on_texts(sentences)
        count_thres = 30
        low_count_words = [w for w,c in tokenizer.word_counts.items() if c < count_thres]
        for w in low_count_words:
            del tokenizer.word_index[w]
            del tokenizer.word_docs[w]
            del tokenizer.word_counts[w]
        res=tokenizer.texts_to_sequences(sentences)
        return(res, tokenizer)

class TfSeqvectorObj(TfDataObj):
    def __init__(self,fasta_file, label_file, batch_size, nepochs, outdir, token_dir ='data/vector_seq_data/'  ):
        fn_pfx= fasta_file.split('/')[-1].split('.')[0]
        encoded_seq_file = token_dir + fn_pfx + '_data.pickle'
        id_file = token_dir + fn_pfx + '_txids.pickle'
        if pathlib.Path(encoded_seq_file).is_file():
            print('loading Data')
            with open(encoded_seq_file, 'rb') as seqfile, open(id_file, 'rb') as idfile:
                padded_encoded_seqs = pickle.load(seqfile)
                ids = pickle.load(idfile)
        else:
            print('encoding Data')
            fasta = [record for record in SeqIO.parse(fasta_file, format = 'fasta')]
            max_seq_size = max([ len(seq) for seq in  fasta])
            ids = [str(record.id) for record in fasta]
            padded_encoded_seqs = np.asarray([ self.encode_and_pad_seq(str(record.seq), max_seq_size) for record in fasta])
            with open(encoded_seq_file, 'wb+') as seqfile, open(id_file, 'wb+') as idfile:
                pickle.dump(padded_encoded_seqs, seqfile, protocol = pickle.HIGHEST_PROTOCOL)
                pickle.dump(ids, idfile,  protocol = pickle.HIGHEST_PROTOCOL)
                

        TfDataObj.__init__(self, padded_encoded_seqs, ids, label_file, batch_size, nepochs, outdir)
        print('Data Loaded')
    def encode_and_pad_seq(self, seq, max_seq_size):
        res = [None]*len(seq)
        '''one hot encod nucleotide array. Any nonstandard nucleotide will be set to n
            https://www.bioinformatics.org/sms/iupac.html
        '''
        for i, nt in enumerate(seq.lower()):
            if nt == 'a':
                res[i]=np.array([1,0,0,0,0])
            elif nt == 'c':
                res[i] = np.array([0,1,0,0,0])
            elif nt == 'g':
                res[i] = np.array([0,0,1,0,0])
            elif nt == 'n':
                res[i] = np.array([0,0,0,1,0])
            elif nt == 't':
                res[i] = np.array([0,0,0,0,1])
            else:
                res[i] = np.array([0,0,0,1,0])
        res = np.array(res)
        pad_size = max_seq_size- len(res) 
        pad = np.zeros((pad_size, res.shape[1]))
        return np.concatenate([res, pad], axis = 0)
