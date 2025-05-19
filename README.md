# drosophila-nsi
This repo contains a model of a Drosophila sensillum containing two insect olfactory neurons intended to investigate non-synaptic 'ephaptic' interactions in the insect olfactory system. The results of using this model are in preparation for publication.  The model was implemented in Brian 2 version 2.6.0 and model equations were integrated using the "milstein" method with $\Delta t = 0.05$ ms timesteps.  

To use this software you will need Brian 2 installed https://brian2.readthedocs.io/en/2.8.0/introduction/install.html

This model simulates two Hodgkin-Huxley-type model neurons (Traub and Miles, 1991).  In order to simulate the electrical separation between the ORN soma and dendrites, each neuron consists of two coupled compartments.  The dendritic compartments are further influenced by an additional transmembrane potential, simulating the arrangement of two ORNs co-housed in a single sensillum.  Spike frequency adaptation is modelled with an M-type adaptation current.

### droso_AB3.py
--
This script is the model with parameters set to emmulate the ORNs in the ab3 sensillum.

### dose_response.py
--
In this computational experiment, designed to assess the dose response of the model, both ORNs are stimulation with 3 second stimulations of increasing intensity.

### NSI_B_to_A_pulses
--
In this computational experiment, designed to examine B-to-A lateral inhibtion, 500 millisecond B pulse stimulations are delivered, offset from the start of 5 second A stimulations, emmulating natural overallping odour filaments.

### droso_AB3_fit.py
--
This script is a version of the model for fitting a number of parameters.

### model_fit_simple.py
--
This script can be used to run further model fitting.

### ab3A_conc_sdf.csv
--
This file contains the data used to fit the model to real odour responses.