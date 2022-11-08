#! /usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
from lstmContacts import *

parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("contactlib",
                    help = "Input contact library file path.",
                    type = str)

parser.add_argument("pset",
                    help = "Prediction set file path.",
                    type = str)

parser.add_argument("-n", "--nfeatures",
                    help = "Number of features to be modeled by the by LSTM.",
                    type = int,
                    default = 101)

parser.add_argument("-u", "--units",
                    help = "Number of units of the LSTM.",
                    type = int,
                    default = 128)

parser.add_argument("-e", "--epochs",
                    help = "Number of epochs of the LSTM.",
                    type = int,
                    default = 100)

parser.add_argument("-t0", "--tin",
                    help = "Length of the input time step.",
                    type = int,
                    default = 5)

parser.add_argument("-t1", "--tout",
                    help = "Length of the predicted time step.",
                    type = int,
                    default = 5)

parser.add_argument("-m", "--method",
                    help = "Method used to merge single contacts time series.",
                    type = str,
                    default = "median")

parser.add_argument("-c", "--cutoff",
                    help = "Median affinity score cutoff.",
                    type = int,
                    default = 88)

parser.add_argument("-p", "--plotname",
                    help = "Time series plot destination.",
                    type = str,
                    default = 'affinity_plot.pdf')

args = parser.parse_args()

tset = contactLibrary(filename = contactlib)
pset = contactLibrary(filename = pset)

model, encoder, decoder = lstmTraining(tset, n = args.nfeatures, units = args.units, epochs = args.epochs)

profile = lstmProfile(pset, encoder, decoder, t0 = args.tin, t1 = args.tout, n = args.nfeatures, method = args.method)

plotProfile(profile, where = args.plotname, threshold = args.cutoff)
