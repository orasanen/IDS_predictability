# These are the scripts used to train LSTM predictive model for F0 contours in Rasanen, Kakouros & Soderstrom (submitted):
# "Is infant-directed speech interesting because it is surprising? –
# Linking properties of IDS to statistical learning and attention at the prosodic level"

import scipy,scipy.io
from keras.models import Sequential, Model
from keras.layers import Dense, Dropout, Activation, Embedding
from keras.layers import LSTM
import keras
from keras.callbacks import Callback
import numpy as np
import os
import sys

# Get path to the training data
datadir = sys.argv[1]

# Load data
loaddata = scipy.io.loadmat('%sLSTM_traindata.mat' % datadir)
train_in = loaddata['train_in']
train_out = loaddata['train_out']
test_in = loaddata['test_in']
test_out = loaddata['test_out']

# early stopping criterion
earlyStopping=keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0.0001, patience=0, verbose=0, mode='auto')

# create LSTM model
model = Sequential()
model.add(Embedding(np.max(train_in)+1, 30, input_length=train_in.shape[1]))
model.add(LSTM(30,activation='sigmoid'))
model.add(Dense(train_out.shape[1], activation='softmax'))

model.compile(loss='categorical_crossentropy', optimizer='rmsprop', metrics=['mean_squared_error'])
print(model.summary())
# Fit model to data
model.fit(train_in, train_out, validation_data=(train_in, train_out), shuffle=True, epochs=150,batch_size=10,callbacks=[earlyStopping])
# get likelihoods for test data
preds = model.predict(test_in)

# save to .mat format for MATLAB reading
scipy.io.savemat('%spreds_out.mat' % datadir,{'preds' : preds})
