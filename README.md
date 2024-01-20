# EEG during Mental Arithmetic
This repository contains the project of **analysis on EEG during mental arithemetic** produced for the Biomedical Signal Processing module of the course BIOMEDICAL SIGNAL PROCESSING AND MEDICAL IMAGES taught by prof. Signorini Maria Gabriella and PHD Steyde Giulio at "Politecnico di Milano" A.Y. 2023/2024.

## Database
The dataset EEGMAT[1] published on PhysioNet [2] is
chosen as a mental arithmetic task. EEGMAT has collected 35
subjects that are during, after three minutes of rest, and one
minute of serial subtraction; this work is focused on the first six
subjects. The obtained signal is composed of 19 channels
sampled with a frequency of 500 Hz

## Goal
Our group amimed to discern the two state of the subject: the resting state and a task(mental arithemetic perdormed) state. The goal must be reached only trough signal elaboration techniques.

## Methods
The analysis was peformed chooseing as feature the area of the Potenrial Spectrum Density(PSD) curve in the theta1, theta2, beta1, and beta2 bands for the CZ, F7, FP2, O2, and P3 EEG channels.
All the pre-processing, the anlysis, and the plotting was carried out through Matlab and for the EEG topoplot was used EEGLab.
For more details and results see [here](docs\abstract.md)

## How to use this application
1. Add your subjects files to [Data](Data\) as `SubjectXX_1.mat` and `SubjectXX_2.mat` where `XX` rapresent the number of the subject, `_1` the rest file, and `_2` the task file.
2. If needed use [filter_creation.m](filter_creation.m) to remake automatically the filters for the different subbands. 
3. Run the [main script](main.m) to performe the analysis and show results.

## References
1. Zyma I, Tukaev S, Seleznov I, Kiyono K, Popov A, Chernykh M,
Shpenkov O. Electroencephalograms during Mental Arithmetic Task
Performance. Data. 2019; 4(1):14. https://doi.org/10.3390/data4010014
2. Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark
RG, Mietus JE, Moody GB, Peng CK, Stanley HE. PhysioBank,
PhysioToolkit, and PhysioNet: Components of a New Research Resource
for Complex Physiologic Signals. Circulation 101(23):e215-e220 [Circulation Electronic Pages;
http://circ.ahajournals.org/content/101/23/e215.full]; 2000 (June 13).
PMID: 10851218; doi: 10.1161/01.CIR.101.23.e215
