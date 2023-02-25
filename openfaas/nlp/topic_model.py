"""
Created on Thu Nov 18 15:06:52 2021

@author: Danielle Lambion
@author: Bob Schmitz
"""

import os
from minio import Minio
import pandas as pd
import pickle
import gensim
from gensim import models
#import nltk
from nltk.stem import WordNetLemmatizer, SnowballStemmer
#import io
# this will be included in the docker image
#nltk.download('wordnet', download_dir='/tmp')
#nltk.download('omw-1.4', download_dir='/tmp')

# Model files
model_files=['/tmp/lda.model',
             '/tmp/lda.model.expElogbeta.npy',
             '/tmp/lda.model.id2word',
             '/tmp/lda.model.state']

class topic_model:
    def __init__(self, mc):
        self.__mc = mc


    def preprocess(self,
                   training_data='/tmp/news_train.csv',
                   bucket_in='topic-modeling',
                   bucket_out='topic-modeling'):
        # =============================================================================
        #     LOAD news_train.csv FROM S3 BUCKET
        #     We will use the last 80% of the dataset for model training
        # =============================================================================
        if not os.path.exists(training_data):
            self.__mc.fget_object(bucket_in, os.path.basename(training_data), training_data)
        df = pd.read_csv(training_data, on_bad_lines='skip',
                         usecols=['publish_date', 'headline_text'])
        df['processed_text'] = df['headline_text'].apply(lambda x: self.process_data(x))
        dictionary = self.create_dict(df['processed_text'])
        corpus_tfidf = self.create_tfidf_model(df['processed_text'], dictionary)
        # =============================================================================
        #     SAVE corpus_tfidf AND dictionary TO S3 BUCKET
        # =============================================================================
        with open('/tmp/dictionary.p', 'wb') as f:
            pickle.dump(dictionary, f)
        # Create an out bucket if it does not exist
        if not self.__mc.bucket_exists(bucket_out):
            self.__mc.make_bucket(bucket_out)
        self.__mc.fput_object(bucket_out, 'dictionary.p', '/tmp/dictionary.p')
        with open('/tmp/corpus_tfidf.p', 'wb') as f:
            pickle.dump(corpus_tfidf, f)
        self.__mc.fput_object(bucket_out, 'corpus_tfidf.p', '/tmp/corpus_tfidf.p')


    def train(self,
              corpus_tfidf='/tmp/corpus_tfidf.p',
              dictionary='/tmp/dictionary.p',
              bucket_out='topic-modeling'):
        # =============================================================================
        #     LOAD corpus_tfidf AND dictionary FROM S3 BUCKET
        # =============================================================================
        if not os.path.exists(corpus_tfidf):
            self.__mc.fget_object(bucket_out, os.path.basename(corpus_tfidf), corpus_tfidf)
        if not os.path.exists(dictionary):
            self.__mc.fget_object(bucket_out, os.path.basename(dictionary), dictionary)
        corpus_tfidf = pickle.load(open(corpus_tfidf, 'rb'))
        dictionary = pickle.load(open(dictionary, 'rb'))
        # DOESN'T WORK IN LAMBDA
        lda_model = models.LdaMulticore(corpus_tfidf, num_topics=5,
                                        id2word=dictionary, passes=2,
                                        workers=8)
        #lda_model = models.LdaModel(corpus_tfidf, num_topics=5, id2word=dictionary)

        # =============================================================================
        #     SAVE lda_model TO S3 BUCKET
        # =============================================================================
        lda_model.save(model_files[0])
        for mfile in model_files:
            self.__mc.fput_object(bucket_out, os.path.basename(mfile), mfile)


    def query(self,
              test_data='/tmp/news_test_smaller.csv',
              dictionary='/tmp/dictionary.p',
              bucket_in='topic-modeling',
              bucket_out='topic-modeling'):
        # =============================================================================
        #     LOAD lda_model AND dictionary AND news_test.csv FROM S3 BUCKET
        #     We will use the last 20% of the dataset to query the model
        # =============================================================================
        if not os.path.exists(test_data):
            self.__mc.fget_object(bucket_in, os.path.basename(test_data), test_data)
        if not os.path.exists(dictionary):
            self.__mc.fget_object(bucket_out, os.path.basename(dictionary), dictionary)
        for mfile in model_files:
            if not os.path.exists(mfile):
                self.__mc.fget_object(bucket_out, os.path.basename(mfile), mfile)
        with open(dictionary, 'rb') as f:
            dictionary = pickle.load(f)
        lda_model = models.LdaModel.load(model_files[0])
        df_query = pd.read_csv(test_data, on_bad_lines='skip',
                               usecols=['publish_date', 'headline_text'])
        df_query['processed_text'] = df_query['headline_text'].apply(lambda x: self.process_data(x))
        query_tfidf = self.create_tfidf_model(df_query['processed_text'], dictionary)
        df_query = self.get_topic(df_query, lda_model, query_tfidf)
        # =============================================================================
        #     SAVE df_query AS A CSV TO S3 BUCKET
        #    (or return it to wherever user might want it)
        # =============================================================================
        results_file = '/tmp/results.csv'
        df_query.to_csv(results_file)
        self.__mc.fput_object(bucket_out, os.path.basename(results_file), results_file)


    # =============================================================================
    # Create a token word dictionary. Tokens that appear in less than 15 headlines
    # are removed. Tokens appearing in more than 50% of the corpus are removed.
    # =============================================================================
    def create_dict(self, docs):
        dictionary = gensim.corpora.Dictionary(docs)
        dictionary.filter_extremes(no_below=15, no_above=0.5)
        return dictionary


    # =============================================================================
    # Create a TFIDF model from a bag-of-words generated by the corpus dictionary.
    # =============================================================================
    def create_tfidf_model(self, docs, dictionary):
        bow_corpus = [dictionary.doc2bow(doc) for doc in docs]
        tfidf = models.TfidfModel(bow_corpus)
        corpus_tfidf = tfidf[bow_corpus]
        return corpus_tfidf


    # =============================================================================
    # Tokenize the String text. Stopwords and words less than 3 characters are
    # removed. Words are stemmed and lemmatized and tokens are returned in
    # their root form.
    # =============================================================================
    def process_data(self, text):
        processed_text = []
        for token in gensim.utils.simple_preprocess(text):
            if token not in gensim.parsing.preprocessing.STOPWORDS and len(token) > 2:
                # lemmatizing verbs
                lemtext = WordNetLemmatizer().lemmatize(token, pos='v')
                # reduce to root form
                stemttext = SnowballStemmer("english").stem(lemtext)
                processed_text.append(stemttext)
        return processed_text


    # =============================================================================
    # Queries the model for the topic number, match score, and topic and appends
    # this information onto the query dataframe.
    # =============================================================================
    def get_topic(self, df, model, tfidf):
        topics_df = pd.DataFrame()
        for tfidf_val in tfidf:
            for index, score in sorted(model[tfidf_val], key=lambda tup: -1*tup[1]):
                topics_df = topics_df.append(pd.Series([index,score,model.print_topic(index, 10)]),
                                             ignore_index=True)
        topics_df.columns = ['topic_number', 'score', 'topic']
        return df.join(topics_df)
