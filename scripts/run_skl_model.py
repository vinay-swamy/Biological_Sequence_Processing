import argparse 
import os 
import pathlib 
import numpy as np 
import pandas as pd 
from sklearn.ensemble import RandomForestClassifier 
from sklearn.model_selection import train_test_split
from data_objs import model_results




parser=argparse.ArgumentParser()
parser.add_argument('--workingDir', action = 'store', dest = 'working_dir')
parser.add_argument('--inputFile', action = 'store', dest = 'input_file')
parser.add_argument('--labFile', action = 'store', dest='lab_file')
parser.add_argument('--nproc', action = 'store', type=int, default=2)
parser.add_argument('--outputDir', action = 'store',  dest='output_dir')
args=parser.parse_args()
os.chdir(args.working_dir)

   
class SklDataObj:
    def __init__(self,x_file_name, lab_file):
        fn=x_file_name.split('/')[-1]
        dm=int(fn.split('_')[1].split('-')[1])
        wc=int(fn.split('_')[2].split('-')[1])
        kmer_size=int(fn.split('_')[3].split('-')[1])
        edim=int(fn.split('_')[4].split('-')[1].split('.')[0] )
        data_name=f'dm-{dm}_wc-{wc}_kmer-{str(kmer_size)}_edim-{str(edim)}'
        X_df=pd.read_csv(x_file_name,names=['embl_id']+ list(range(edim)))
        labs=pd.read_csv(lab_file ).rename({'family_id' : 'target_label'}, axis =1)
        X_df_labeled=pd.merge(left=labs, how='inner', right=X_df, left_on='embl_id', right_on='embl_id')
        X_data=np.asarray(X_df_labeled.iloc[:,3:])#drop the first 3 columns
        labs = X_df_labeled.iloc[:,:3]
        Y_vec=np.asarray(X_df_labeled['target_label'])
        X_train, X_test, train_labs, test_labs =train_test_split(X_data,labs,test_size=.2, 
                                                            random_state=42,stratify=Y_vec)
        self.X_train=X_train
        self.Y_train=np.asarray(train_labs.iloc[:,2])
        self.train_labs = train_labs
        self.X_test=X_test
        self.Y_test=np.asarray(test_labs.iloc[:,2])
        self.test_labs = test_labs
        self.name=data_name
        self.model=None
    def summary(self):
        tr_len=len(self.X_train)
        ts_len=len(self.X_test)
        print(f'Training size: {tr_len}\nClass Counts:')
        print(self.train_labs.protein_families.value_counts())
        print(f'Test size: {ts_len}\nClass Counts:')
        print(self.test_labs.protein_families.value_counts())
    def run_model(self, model, model_name, outdir):
        model.fit(self.X_train, self.Y_train)
        self.model=model 
        Y_pred_class = model.predict(self.X_test)
        Y_pred_prob = model.predict_proba(self.X_test)
        Y_true = self.Y_test
        model_results( Y_true, Y_pred_class, Y_pred_prob, f'{self.name}_{model_name}', outdir)
        #return model_res_line

outdir = args.output_dir
if outdir[-1] is not '/':
    outdir+='/'
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

with open(outdir+'model_results.csv', 'w+') as model_res_file:
    data_obj = SklDataObj(args.input_file, args.lab_file)
    data_obj.summary()
    rf_model = RandomForestClassifier(n_estimators=100, random_state=32, n_jobs=args.nproc)
    data_obj.run_model(rf_model, 'random_forest', outdir)
    












