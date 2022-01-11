# Protein Liquid-liquid Phase Separation Predictor Website (JFK laboratory, SickKids Hospital, Toronto)

Website URL: http://abragam.med.utoronto.ca/~JFKlab/Software/psp.htm

This program predicts the likelihood of intrinsically disordered protein regions (IDRs) to phase separate based on propensity for long-range planar pi-pi contacts. The phase separation predictor returns a single score per sequence supplied in fasta sequence format, ignoring sequences that are shorter 140 residues (shorter than any sequence in the training/test set) or that have ambigious residues (outside the scope of the training process). This predictor was designed primarily with the goal of testing the relationship between planar pi-pi contacts and IDR phase separation, but validation of the predictor demonstrates that it discriminates known phase-separating IDR-containing proteins from other protein sequences.

For use in predicting sequences note that the contact predictions and validation tests both derive primarily from sequences found in nature, and the supported use of this program is in highlighting full sequences as found in proteomic datasets. Scores for artificial sequences and for arbitrary segments of proteins represent untested extrapolations.

The main algorithm and software is developed and written by Dr. Robert M. Vernon. The implementation and deployment of algorighm, website and server is conducted by Dr. Hao Cai.

Reference: Pi-Pi contacts are an overlooked protein feature relevant to phase separation R. Vernon et al.
