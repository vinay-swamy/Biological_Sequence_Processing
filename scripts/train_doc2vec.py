#!/usr/bin/env python
# coding: utf-8

# In[1]:
from gensim.models.doc2vec import Doc2Vec, TaggedDocument
import sys
import pickle
import logging
corpus_file=sys.argv[1]
outfile=sys.argv[2]


# In[ ]:
logging.basicConfig(format='%(levelname)s : %(message)s', level=logging.INFO)
logging.root.level = logging.INFO



#ref_data_file='../run_k-6_l-2/ref_transcript_kmers.pydata'
#dummy_data_file='../run_k-6_l-2/dummy_transcript_kmers.pydata'




# In[ ]:


#PV-DBOW == skipgram == `dm=0`
print('training')
model = Doc2Vec(corpus_file=corpus_file, dm=0, vector_size=300, min_count=3, epochs=15, seed=42, workers=64)


# In[ ]:


model.save(outfile)

