import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import pygsti;reload(pygsti)
print 'everything reloaded'
import pickle
import re
#import qutip
reload(mbi)


def seq_append(split):
    """
    Help function to make single fiducial / germ sequences
    """
    string = ''

    for i in range(len(split)):
        if split[i] == "Gi":
            string +='e'
        elif split[i] == "Gx":
            string +='x'
        elif split[i] == "Gy":
            string +='y'   
        elif split[i] == "Gu":
            string +='u' 
        elif split[i] == "Gv":
            string +='v' 
            
    if string=='':
        string +='e'
        
    return [string]
        

def create_experiment_list_pyGSTi(filename):
    """
    Extracting list of experiments from .txt file
    Parameters:
    filename: string
        Name of the .txt file. File must be formatted in the way as done by
        pyGSTi.
        One gatesequence per line, formatted as e.g.:Gx(Gy)^2Gx. This is usually created
        with the simulateGST ipython notebook, and saved as MyDataTemplate.txt
    Returns:
    Nested list of single experiment in a format the code below can understand it 
    """

    try:
        experiments = open(filename)
        sequences = experiments.read().split("\n")
        experiments.close()   
    except:
        print "Could not open the filename you provided."
        return
    
    experimentlist = []

    #A parameter to label our experimental runs 
    a = 1

    for i in range(len(sequences)):
        clean_seq = sequences[i].strip()
        gateseq = []
        fiducial = []
        germs = []
        measfiducial = []
        power = 0
        
        # Special case of no fiducials and no germs
        if "{}" in clean_seq:
            fiducial.extend('e')
            measfiducial.extend('e')
            germs.extend('e')
            power = 1
        
        # Case of a germ sequence with at least one fiducial around
        elif "(" in clean_seq:
            #Case of repeated germs
            if "^" in clean_seq:
                power = int(re.findall("\d+", clean_seq)[0])
                result = re.split("[(]|\)\^\d", clean_seq)
            
            #Only one germ
            else:
                power = 1
                result = re.split("[()]", clean_seq)
                
            fiducial.extend(seq_append(re.findall("G[xyuvi]", result[0])))
            germs.extend(seq_append(re.findall("G[xyuvi]", result[1])))
            measfiducial.extend(seq_append(re.findall("G[xyuvi]", result[2])))

        #Otherwise it is a single germ sequence without fiducials
        elif ("Gi" in clean_seq) or ("Gx" in clean_seq) or ("Gy" in clean_seq) or ("Gu" in clean_seq) or ("Gv" in clean_seq) and ("(" not in clean_seq) :
            power = 1
            
            fiducial.extend('e')
            germs.extend(seq_append(re.findall("G[xyuvi]", clean_seq)))
            measfiducial.extend('e')

        #Make sure we only extend the experimentlist in case we did not look at a comment line or empty line
        if power !=0:
            gateseq.extend(fiducial)           
            gateseq.extend(germs)  
            gateseq.extend(measfiducial)
            gateseq.append(power)
            gateseq.append(a)
            
            a+=1

            experimentlist.append(gateseq)
       
    #Make sure everything worked. -2 for sequences as first line is comments and last line is empty
    if len(experimentlist) < (len(sequences)-2):
        print("Lenght list of experiments too short, something went wrong. Check your mess.")
        return
    
    return experimentlist



def analysis_gateset(contains = '',timestamp=None,ssro_folder = None,
            ssro_tstamp ='',base_folder = None, correct_ionization = True):
    ''' 
    Function to analyze gate set tomography data, acquired with the individual_awg-write function
    '''
    if timestamp != None:
        folder = toolbox.data_from_time(timestamp,folder = base_folder)
    else:
        folder = toolbox.latest_data(contains,folder = base_folder)
    print "analyis", folder
    if ssro_folder == None:
        if ssro_tstamp == '':
            ssro_folder = ssro_tstamp

        else:
            ssro_folder = toolbox.data_from_time(ssro_tstamp)
    
    a = mbi.MBIAnalysis(folder,CR_after_check = False)

    measurement_name = a.g.attrs['all_run_names']   


    sweep_pts = []
    ssro_results_no_correction = []
    ssro_results_w_correction = []
    repetitions = []

    for k in range(0,len(measurement_name)):
        sweep_pts.append(a.get_sweep_pts())
        a.get_readout_results(name=str(measurement_name[k]),CR_after_check = correct_ionization)
        a.get_electron_ROC(ssro_folder)

        result_no_RO_correction = a.ssro_results.flatten()
        result_w_RO_correction  = a.p0.reshape(-1)[:]
        ssro_results_no_correction.extend(result_no_RO_correction)
        ssro_results_w_correction.extend(result_w_RO_correction)

        for i in range(len(result_no_RO_correction)): 
            repetitions.append(a.reps)

    neg_no_correction = np.array(repetitions)-np.array(ssro_results_no_correction)
    neg_w_correction  = (1-np.array(ssro_results_w_correction))*np.array(repetitions)
    ssro_results_w_correction *=np.array(repetitions)
    neg_w_correction = neg_w_correction.astype(int)
    ssro_results_w_correction = ssro_results_w_correction.astype(int)

    insert_counts_into_dataset(str(folder)+"/MyDataTemplate.txt", str(folder)+"/Dataset_no_RO_correction.txt", ssro_results_no_correction, neg_no_correction)
    insert_counts_into_dataset(str(folder)+"/MyDataTemplate.txt", str(folder)+"/Dataset_w_RO_correction.txt", ssro_results_w_correction, neg_w_correction)

    return

def insert_counts_into_dataset(filename_gateseq, filename_dataoutput, pos_counts, neg_counts):
    """
    Function for inserting counts retreived from experiment
    into pyGSTi. Counts are added to .text file which was used
    to create the gate sequences for the experiment.
    Formatted as pyGSTi .txt style
    Parameters:
    ----------------------------------------------------
    filename_gatseq: string referring to a .txt file
    must be a .txt file
    filename_dataoutput: string 
    This is the name of the output file
    positive counts: list of the positive counts
    negative coutns: list of the negative counts

    """

    #Open the experiment list and split it in lines. This thus also includes
    # the first line(title), and last whiteline
    #After that close the file again
    experiment_text_file_gateseq = open(filename_gateseq)
    sequences = experiment_text_file_gateseq.read().split("\n") 
    experiment_text_file_gateseq.close()

    f = open(filename_dataoutput, "w")
    print >> f, "## Columns = plus count, minus count"
    for i in range(len(sequences)-2):
            if i == 0:
                print >> f, "{}"+"  %s" % pos_counts[0]+"  %s" % neg_counts[0]
                # gatesquence + 2 spaces + # count_plus+2 spaces+count_min
            else:
                #The clean sequence removes all the zero entries from the dataset, this is 6 characters long by default
                cleanseq = sequences[i+1][:-6]
                print >> f, cleanseq+'  %s' % pos_counts[i]+'  %s' % neg_counts[i]
                #gatesquence + 2 spaces + count_plus+2 spaces+count_min
    f.close()


def create_pygsti_report(file_target_gateset, file_prep_fiducials, file_meas_fiducials, file_germs, file_max_lengths, file_data, contains = '', timestamp = None, base_folder = None,
                         gauge_ratio = 1.0e-3, gate_ratio = 1, constrainToTP=True, full_pdf = False, easy_pdf = False):
    """
    This function is used to create a gateset tomography report in pdf format using the pygsti package.
    Inputs are:
        -   file_target_gateset text file specifying the target gate sets that we wanted to perform
        -   file_prep_fiducials text file specifying the used preparation fiducials
        -   file_meas_fiducials text file specifying the used measurement fiducials
        -   file_germs          text file specifying the used germs
        -   file_data           text file that cotains the data we want to analyze, created with the analysis_gateset function

    Some notes about the results function:

    Typically we want to constrain the resulting gates to be trace-preserving, so we leave constrainToTP set to True (the default).
    gaugeOptRatio specifies the ratio of the state preparation and measurement (SPAM) weight to the gate weight when performing a gauge
    optimization. When this is set to 0.001, as below, the gate parameters are weighted 1000 times more relative to the SPAM parameters. 
    Typically it is good to weight the gates parameters more heavily since GST amplifies gate parameter errors via long gate sequences
    but cannot amplify SPAM parameter errors. If unsure, 0.001 is a good value to start with.

    By setting full_pdf and/or easy_pdf to true, according pdf results are created by pyGSTi that visualize the results


    TO DO: Currently supports readout corrected data as input, but can not analyze it as bound to be either zero or 1 for readout
            results by pygsti analysis
    """

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp, folder = base_folder)
    else:
        folder = toolbox.latest_data(contains, folder = base_folder)
    print "create folder", folder

    #Load in the parameters how we took the data
    gs_target       =  pygsti.io.load_gateset(str(folder)+str(file_target_gateset))
    fiducials_prep  =  pygsti.io.load_gatestring_list(str(folder)+str(file_prep_fiducials))
    fiducials_meas  =  pygsti.io.load_gatestring_list(str(folder)+str(file_meas_fiducials))
    germs           =  pygsti.io.load_gatestring_list(str(folder)+str(file_germs))
    maxLengths      =  pickle.load(open(str(folder)+str(file_max_lengths)))

    #Calculate the results. 
    if constrainToTP:
        gs_target.set_all_parameterizations("TP")
    results = pygsti.do_long_sequence_gst(str(folder)+str(file_data), gs_target, 
                                        fiducials_prep, fiducials_meas, germs, maxLengths,
                                        gaugeOptParams={'itemWeights': {'spam': gauge_ratio, 'gates': gate_ratio}} )
   
    #Make sure we save the report pdf that we create in the following section with a name that specifies if we used the readout 
    #corrected results or not.
    save_name = "_with_RO_correction"

    if 'no_RO_correction' in file_data:
        save_name = "_no_RO_correction"

    pickle.dump( results, open(str(folder)+"/"+save_name+str(".pkl"), "wb"))
    s = pickle.dumps(results)
    r2 = pickle.loads(s)

    final_estimate =  r2.gatesets['final estimate']
    print "The final estimates are", final_estimate, final_estimate['rho0'], type(final_estimate['rho0']), final_estimate['rho0'][0]



    #Create pdf reports that visualize the results
    if easy_pdf:
        #create a brief GST report (just highlights of full report but fast to generate; best for folks familiar with GST)
        results.create_brief_report_pdf(confidenceLevel=95, filename=str(folder)+str("/brief_report")+save_name+str(".pdf"), verbosity=2)

    if full_pdf:
        #create a full GST report (most detailed and pedagogical; best for those getting familiar with GST)
        results.create_full_report_pdf(confidenceLevel=95, filename=str(folder)+str("/full_report")+save_name+str(".pdf"), verbosity=2)

    return

def simulator(rho, eo, x, y, u, v, experiment):
    exp_i = []
    #Stack together the preparation fiducial, germs and measurement fiducials
    for k in range(3):
        exp_i.extend(experiment[k])

    for l in range(len(exp_i)):
        if l == 0:
            #This will be the case where we need to use rho*first one. Might need subfunction to execute some checking, or can we build up and use a dictionary?
            print l
        else:
            #This will be the case where we multiply with another one
            print l

    return #this will be where we return the final probabilities


def simulate_gateset_results(file_target_gateset, file_prep_fiducials, file_meas_fiducials, file_germs, file_max_lengths, file_data, constrainToTP=True, contains = '', timestamp = None, base_folder = None,
                         gauge_ratio = 1.0e-3, gate_ratio = 1):

    if timestamp != None:
        folder = toolbox.data_from_time(timestamp, folder = base_folder)
    else:
        folder = toolbox.latest_data(contains, folder = base_folder)
    print "simulate folder", folder

    #Make sure we save the report pdf that we create in the following section with a name that specifies if we used the readout 
    #corrected results or not.
    save_name = "_with_RO_correction"

    if 'no_RO_correction' in file_data:
        save_name = "_no_RO_correction"

    try:
        results = pickle.load(open(str(folder)+"/"+save_name+str(".pkl")))

    except: 
        create_pygsti_report(file_target_gateset, file_prep_fiducials, file_meas_fiducials, file_germs, file_max_lengths, file_data, contains = contains, timestamp = None, base_folder = base_folder,
                         gauge_ratio = 1.0e-3, gate_ratio = 1, constrainToTP=True, full_pdf = False, easy_pdf = False)

        results = pickle.load(open(str(folder)+"/"+save_name+str(".pkl")))


    experiments = create_experiment_list_pyGSTi(str(folder)+str(file_data))


    #This is where we will use qutip and the GST results to define the operators


    # rho = qutip.Qobj(qt.Qobj([[1.3, 0.008],
    #                             [0.008, -0.3]]))

    # eo = qutip.Qobj(qt.Qobj([[0.7, 0.008],
    #                             [0.008, 0.18]]))

    # x = qutip.Qobj(qt.Qobj([[1, 0, 0, 0],
    #                         [0, 1, 0, 0],
    #                         [0, 0, 0, -1],
    #                         [0, 0, 1, 0]], type='super'))

    simulated_data = []

    for i in range(len(experiments)):
       
        #Get the experiments, arranged that [0] fiducial prep, [1] germ [2 meas fiducial] 
        simulated_data.append(simulator(rho=0, eo=0, x=0, y=0, u=0, v=0, experiment=0))

    simulated_data = np.array(simulated_data)

    #This will be where we will implement a checking algorithms that checks for outliers between the simulated results and the acquired dataset. Can we also get this x2 values in here. that 
    #would be sweet to facilitate comparison. 



a = np.array([0.7, 0, 0, 0.7])

a = np.transpose(a)

b = np.array([[1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0,-1, 0],
    [0, 0, 0, 1]])

print np.dot(b, a)
