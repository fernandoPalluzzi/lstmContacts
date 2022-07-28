#! /usr/bin/python3
# -*- coding: utf-8 -*-

import re
import pandas
from random import randint, uniform
from numpy import median
from numpy import argmax
from numpy import array
from numpy import array_equal
from keras.layers import Input
from keras.layers import LSTM
from keras.layers import Dense
from keras.models import Model
from keras.utils import to_categorical



#filename = "/home/fernando/contactsCore/contactLibrary_t5.txt"
#variant = "beta"
#antibody = "7kmg"
#residue = "h.R50"
#sep = "\t"
#lineterm = "\n"
#group = 1

#clib.loc[(clib['antibody'] == "7kmg") & (clib['variant'] == "wt")]
#clib[["a", "y"]][(clib['antibody'] == "7kmg") & (clib['variant'] == "wt")]
#clib[["a", "y"]][(clib['antibody'].isin(["7kmg"])) & (clib['variant'].isin(["wt"]))]


#filename = "~/contactsCore/contactLibrary_t5.txt"
#variant = "beta, omicron, delta, wt"
#antibody = "7kmg, 7cm4, 7l7d"
#residue = None
#group = None
#sep = "\t"
#lineterm = "\n"


def contactLibrary(clib = None, filename = None, variant = None,
                   antibody = None, residue = None, group = None,
                   sep = "\t", lineterm = "\n", keepAttrs = True):
	if variant != None:
		if ", " in variant:
			variant = re.sub(", ", "\", \"", variant)
		variant = "[\"" + variant + "\"]"
	if antibody != None:
		if ", " in antibody:
			antibody = re.sub(", ", "\", \"", antibody)
		antibody = "[\"" + antibody + "\"]"
	if residue != None:
		if ", " in residue:
			residue = re.sub(", ", "\", \"", residue)
		residue = "[\"" + residue + "\"]"
	if clib == None:
		clib = pandas.read_csv(filename, sep = sep, lineterminator = lineterm)
	c = []
	if variant == None and antibody == None and residue == None and group == None:
		L = clib[['a', 'y']]
		V = None
	else:
		if variant != None:
			c += ["(clib['variant'].isin(" + variant + "))"]
		if antibody != None:
			c += ["(clib['antibody'].isin(" + antibody + "))"]
		if residue != None:
			c += ["(clib['res'].isin(" + residue + "))"]
		if group != None:
			c += ["(clib['group'].isin(" + str(group) + "))"]
			g = "[(clib.iloc[w]['group'].isin(['" + str(group) + "']))]"
		else:
			g = ""
		if keepAttrs:
			V = eval("clib[" + ' & '.join(c) + "]")
			keep = "clib.iloc[w]"
		else:
			V = eval("clib[['a', 'y']][" + ' & '.join(c) + "]")
			keep = "clib.iloc[w][['a', 'y']]"
		u = [i for i in clib.index]
		v = [i for i in V.index]
		w = [i for i in u if i not in v]
		L = eval(keep + g)		
	return L, V


#t = 50
#n0 = 6
#n1 = 3


def randomSequence(length, u):
	return [randint(1, u - 1) for _ in range(length)]


def getRandomset(n0, t, n, n1 = 0, w = 0):
	x0, x1, y = list(), list(), list()
	if n1 == 0:
		n1 = n0
	
	for _ in range(n):
		
		# Generating random data
		source = randomSequence(n0, t)
		if w > 0:
			target = [a + uniform(-w, w) for a in source]
		else:
			target = source[:n1]
			target.reverse()
		lagged = [0] + target[:-1]
		
		# Encoding
		encoded0 = to_categorical([source], num_classes = t)[0]
		encoded1 = to_categorical([lagged], num_classes = t)[0]
		encoded2 = to_categorical([target], num_classes = t)[0]
		
		# Storing
		x0.append(encoded0)
		x1.append(encoded1)
		y.append(encoded2)
		
	return array(x0), array(x1), array(y)


def getDataset(x, b = 100, t = 100, h = 0, rounded = False):
	x0, x1, y = list(), list(), list()
	
	for i in x.index:
		
		# Reading data
		if rounded:
			source = [round(b*float(a) + h) for a in x['a'][i].split(',')]
			target = [round(b*float(y) + h) for y in x['y'][i].split(',')]
		else:
			source = [b*float(a) + h for a in x['a'][i].split(',')]
			target = [b*float(y) + h for y in x['y'][i].split(',')]
		lagged = [0] + target[:-1]
		
		# Encoding
		encoded0 = to_categorical([source], num_classes = t + h + 1)[0]
		encoded1 = to_categorical([lagged], num_classes = t + h + 1)[0]
		encoded2 = to_categorical([target], num_classes = t + h + 1)[0]
		
		# Storing
		x0.append(encoded0)
		x1.append(encoded1)
		y.append(encoded2)
	
	return array(x0), array(x1), array(y)


def onehotEncode(seq, u):
	encoded = list()
	for x in seq:
		v = [0 for _ in range(u)]
		v[x] = 1
		encoded.append(v)
	return array(encoded)


def onehotDecode(encoded):
	return [argmax(x) for x in encoded]


def define_models(n0, n1, u):
	
	# Training encoder
	encoder_inputs = Input(shape = (None, n0))
	encoder = LSTM(u, return_state = True)
	encoder_outputs, state_h, state_c = encoder(encoder_inputs)
	encoder_states = [state_h, state_c]
	
	# Training decoder
	decoder_inputs = Input(shape = (None, n1))
	decoder_lstm = LSTM(u, return_sequences = True, return_state = True)
	decoder_outputs, _, _ = decoder_lstm(decoder_inputs, initial_state = encoder_states)
	decoder_dense = Dense(n1, activation = 'softmax')
	decoder_outputs = decoder_dense(decoder_outputs)
	model = Model([encoder_inputs, decoder_inputs], decoder_outputs)
	
	# Inference
	encoder_model = Model(encoder_inputs, encoder_states)
	decoder_state_input_h = Input(shape = (u,))
	decoder_state_input_c = Input(shape = (u,))
	decoder_states_inputs = [decoder_state_input_h, decoder_state_input_c]
	decoder_outputs, state_h, state_c = decoder_lstm(decoder_inputs, initial_state = decoder_states_inputs)
	decoder_states = [state_h, state_c]
	decoder_outputs = decoder_dense(decoder_outputs)
	decoder_model = Model([decoder_inputs] + decoder_states_inputs, [decoder_outputs] + decoder_states)
	
	return model, encoder_model, decoder_model


def predict_sequence(encoder, infdec, source, t, n):
	
	# Encoding
	state = encoder.predict(source)
	target = array([0.0 for _ in range(n)]).reshape(1, 1, n)
	
	# Predictions
	output = list()
	for i in range(t):
		yhat, h, c = decoder.predict([target] + state)
		output.append(yhat[0,0,:])
		
		# Update
		state = [h, c]
		target = yhat
	
	return array(output)





# -- RUN ------------------------------------------------------------ #

n_features = 101
t0 = 5
t1 = 5

# Model
model, encoder, decoder = define_models(n_features, n_features, 128)
model.compile(optimizer = 'adam', loss = 'categorical_crossentropy',
              metrics = ['accuracy'])


# Training set

"""
X1, X2, y = getRandomset(n0 = t0, n1 = t1, t = n_features, n = 26100)
#X1, X2, y = getRandomset(n0 = t0, n1 = t1, t = n_features, n = 26100, w = 1)
print(X1.shape, X2.shape, y.shape)
"""

#clib = contactLibrary(filename = "~/contactsCore/contactLibrary_t5.txt")
#clib = contactLibrary(filename = "~/contactsCore/contactLibrary_t5.txt", variant = "beta", antibody = "7kmg")

#X1, X2, y = getDataset(clib[0])
#print(X1.shape, X2.shape, y.shape)

#V1, V2, Vy = getDataset(clib[1])
#print(V1.shape, V2.shape, Vy.shape)

#tset = contactLibrary(filename = "~/contactsCore/contactLibrary_t5.txt", group = 1)
#tset = contactLibrary(filename = "~/contactsCore/contactLibrary_t5.txt", variant = "beta")
#tset = contactLibrary(filename = "~/contactsCore/contactLibrary_t5.txt", variant = "beta", antibody = "7kmg")
tset = contactLibrary(filename = "~/contactsCore/contactLibrary_t5.txt", variant = "beta, omicron, delta, wt", antibody = "7kmg, 7cm4, 7l7d")
#tset = contactLibrary(filename = "~/contactsCore/contactLibrary_t5.txt", variant = "wt, delta, omicron", antibody = "7l7d")

X1, X2, y = getDataset(tset[0])
#X1, X2, y = getDataset(tset[1])
print(X1.shape, X2.shape, y.shape)


# Model training

#model.fit([X1, X2], y, epochs = 5)
#model.fit([X1, X2], y, epochs = 10)
model.fit([X1, X2], y, epochs = 100)
#model.fit([V1, V2], Vy, epochs = 80)
#model.fit([V1, V2], Vy, epochs = 100)

# Performances

"""
total, correct = 100, 0

for _ in range(total):
	X1, X2, y = getRandomset(n0 = t0, n1 = t1, t = n_features, n = 1)
	target = predict_sequence(encoder, decoder, X1, t1, n_features)
	if array_equal(onehotDecode(y[0]), onehotDecode(target)):
		correct += 1

print('Accuracy: %.2f%%' % (float(correct)/float(total)*100.0))
"""


vset = contactLibrary(filename = "~/contactsCore/contactLibrary_t5.txt", variant = "beta", antibody = "7kmg")
#vset = contactLibrary(filename = "~/contactsCore/contactLibrary_t5.txt", variant = "wt", antibody = "7l7d")
#vset = contactLibrary(filename = "~/contactsCore/contactLibrary_t5.txt", variant = "delta", antibody = "7l7d")

V1, V2, Vy = getDataset(vset[1])
print(V1.shape, V2.shape, Vy.shape)


N, T = len(V1), 0
maxdist = 5
maxsum = 15
profile = []

for i in range(N):
	target = predict_sequence(encoder, decoder, V1[[i]], t1, n_features)
	y_hat = onehotDecode(target)
	profile += y_hat
	y_encoded = onehotDecode(Vy[[i]][0])
	d = [abs(y_hat[j] - y_encoded[j]) for j in range(len(y_hat))]
	if max(d) <= maxdist and sum(d) <= maxsum:
		T += 1

print('Accuracy:', str(round(100*float(T)/float(N), 3)) + '%')

chunk = 100
profile = [profile[i:i+chunk] for i in range(0, len(profile), chunk)]
profile = array([array(x) for x in profile])
#profile = mean(profile, axis = 0)
profile = median(profile, axis = 0)
















total, correct = len(X1[:100]), 0

for i in range(total):
	target = predict_sequence(infenc, infdec, V1[[i]], n_steps_out, n_features)
	if array_equal(onehotDecode(y[[i]][0]), onehotDecode(target)):
		correct += 1

print('Accuracy: %.2f%%' % (float(correct)/float(total)*100.0))


# Prediction

#for _ in range(10):
#	X1, X2, y = get_dataset(n_steps_in, n_steps_out, n_features, 1)
#	target = predict_sequence(infenc, infdec, X1, n_steps_out, n_features)
#	print('X=%s y=%s, yhat=%s' % (one_hot_decode(X1[0]), one_hot_decode(y[0]), one_hot_decode(target)))



