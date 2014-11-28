

QEC_dict = {}


for state in ['Z']:
    logic_state = state
    for RO in range(7):
	    for k in range(4):
		    for error_sign in [1,-1]:
		    	for direction in ['positive','negative']:
	    			timestamp, folder = toolbox.latest_data(contains = direction+'_RO'+str(RO)+'_k'+str(k)+'_sign'+ str(error_sign)+'_'+logic_state, older_than = older_than,return_timestamp = True)
	    			SSRO_timestamp, SSRO_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = timestamp,return_timestamp = True)

	    			QEC_dict[state]['TOMO_'+str(RO)][str(error_sign)]['k_'+str(k)]['direction'] = Plot_QEC(timestamp = timestamp, measurement_name = measurement_name, folder_name = folder,
			            ssro_calib_timestamp = SSRO_timestamp) 