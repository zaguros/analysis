import cPickle as pickle
import time #Used for timing purposes

filenames =  ['simulation_data_weekend/NV_C13_0.011_Raw_Data_20131205_1740',
     'simulation_data_weekend/NV_C13_0.003_Raw_Data_20131206_0032',
     'simulation_data_weekend/NV_C13_0.0011_Raw_Data_20131207_0258']
t_start = time.time()

for idfn, filename in enumerate(filenames) :
    print 'loading data from file: ' +str(filename)
    file = open(filename,'rU')
    B_Fields = pickle.load(file)
    NV_List = pickle.load(file)
    Raw_Data = pickle.load(file)
    filenameConv = filename+ '_bin'
    file= open( filenameConv, 'wb') #Shelving would be cleaner but so far this works fine
    pickle.dump(B_Fields,file)
    pickle.dump(NV_List,file)
    pickle.dump(Raw_Data,file)

    print 'Converting data took' +str(time.time()-t_start)
