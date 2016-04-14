
from numpy import *
import pylab as plt
from matplotlib import rc
#import plots


E6IIno2={
'SIL45':{'sat_power_lt1':241,'sat_cnts_MM_lt1':450},
'SILx1':{'sat_power_lt1':278,'sat_cnts_MM_lt1':387},
'SILx2':{'sat_power_lt1':281,'sat_cnts_MM_lt1':553},
'SILx3':{'sat_power_lt1':405,'sat_cnts_MM_lt1':431},
'SIL35':{'sat_power_lt1':100,'sat_cnts_MM_lt1':789},
'SIL34':{'sat_power_lt1':419,'sat_cnts_MM_lt1':440},
'SIL59':{'sat_power_lt1':165,'sat_cnts_MM_lt1':554},
'SIL58':{'sat_power_lt1':83,'sat_cnts_MM_lt1':596},
'SIL57':{'sat_power_lt1':276,'sat_cnts_MM_lt1':574},
'SIL56':{'sat_power_lt1':121,'sat_cnts_MM_lt1':590},
'SIL38':{'sat_power_lt1':223,'sat_cnts_MM_lt1':618},
'SIL36':{'sat_power_lt1':158,'sat_cnts_MM_lt1':732},
'SIL54':{'sat_power_lt1':319,'sat_cnts_MM_lt1':533},
'SILx4':{'sat_power_lt1':305,'sat_cnts_MM_lt1':411},
}


the111no1_before_dict={
'SIL10':{'sat_power_lt2':141,'sat_cnts_MM_lt2':523,'sat_cnts_SM_lt2':150},
'SIL9':{'sat_power_lt2':90,'sat_cnts_MM_lt2':420,'sat_cnts_SM_lt2':116},
'SIL15':{'sat_power_lt2':135,'sat_cnts_MM_lt2':300},
'SIL18':{'sat_power_lt2':47,'sat_cnts_MM_lt2':396},
'SIL7':{'sat_power_lt2':120,'sat_cnts_MM_lt2':303},
'SIL8':{'sat_power_lt2':92,'sat_cnts_MM_lt2':356},
'SIL12':{'sat_cnts_MM_lt2':350},
'SIL1':{'sat_cnts_MM_lt2':420},
}
the111no1_after_dict={
'SIL1':{'sat_power_lt3':70,'sat_cnts_MM_lt3':690,'sat_cnts_SM_lt3':160},
'SIL7':{'sat_power_lt3':160,'sat_cnts_MM_lt3':636,'sat_cnts_SM_lt3':130},
'SIL8':{'sat_power_lt3':30,'sat_cnts_MM_lt3':730,'sat_cnts_SM_lt3':280},
'SIL9':{'sat_power_lt3':55,'sat_cnts_MM_lt3':676,'sat_cnts_SM_lt3':200},
'SIL12':{'sat_power_lt3':25,'sat_cnts_MM_lt3':770,'sat_cnts_SM_lt3':260},
'SIL18':{'sat_power_lt3':50,'sat_cnts_MM_lt3':771,'sat_cnts_SM_lt3':280},
}

the111no2_before_dict={
'SIL1':{'sat_power_lt1':475,'sat_cnts_MM_lt1':278},
'SIL2':{'sat_power_lt1':200,'sat_cnts_MM_lt1':530,'sat_cnts_MM_lt3':1150,'sat_cnts_SM_lt3':350},
'SIL3':{'sat_power_lt1':90,'sat_cnts_MM_lt1':643},
'SIL4':{'sat_power_lt1':0,'sat_cnts_MM_lt1':0},
'SIL5':{'sat_power_lt1':231,'sat_cnts_MM_lt1':470},
'SIL6':{'sat_power_lt1':1521,'sat_cnts_MM_lt1':140},
'SIL7':{'sat_power_lt1':187,'sat_cnts_MM_lt1':532},
'SIL8':{'sat_power_lt1':120,'sat_cnts_MM_lt1':571},
'SIL9':{'sat_power_lt1':730,'sat_cnts_MM_lt1':306},
'SIL10':{'sat_power_lt1':1233,'sat_cnts_MM_lt1':53},
'SIL11':{'sat_power_lt1':0,'sat_cnts_MM_lt1':0},
'SIL12':{'sat_power_lt1':79,'sat_cnts_MM_lt1':592},
'SIL13':{'sat_power_lt1':232,'sat_cnts_MM_lt1':312},
'SIL14':{'sat_power_lt1':91,'sat_cnts_MM_lt1':590},
'SIL15':{'sat_power_lt1':1163,'sat_cnts_MM_lt1':104},
'SIL16':{'sat_power_lt1':266,'sat_cnts_MM_lt1':405},
'SIL17':{'sat_power_lt1':806,'sat_cnts_MM_lt1':211},
'SIL18':{'sat_power_lt1':685,'sat_cnts_MM_lt1':241},
}

the111no2_after_dict={
'SIL1':{'sat_power_rt2':40,'sat_cnts_SM_rt2':280,'sat_cnts_MM_rt2':990,'sat_cnts_MM_lt3':1200},
'SIL2':{'sat_power_rt2':30,'sat_cnts_SM_rt2':280,'sat_cnts_MM_rt2':970,'sat_cnts_MM_lt3':1300,'sat_cnts_SM_lt3':650},
'SIL5':{'sat_power_rt2':35,'sat_cnts_SM_rt2':0,'sat_cnts_MM_rt2':790,'sat_cnts_MM_lt3':1200},
'SIL14':{'sat_power_rt2':0,'sat_cnts_SM_rt2':0,'sat_cnts_MM_rt2':0},
'SIL16':{'sat_power_rt2':120,'sat_cnts_SM_rt2':170,'sat_cnts_MM_rt2':920,'sat_cnts_MM_lt3':80},
}


Pippin_noAR={
'SIL1':{'sat_power_rt2':130,'sat_cnts_SM_rt2':310,'sat_cnts_MM_rt2':1170},
'SIL2':{'sat_power_rt2':60,'sat_cnts_SM_rt2':360,'sat_cnts_MM_rt2':1020},
'SIL3':{'sat_power_rt2':65,'sat_cnts_SM_rt2':290,'sat_cnts_MM_rt2':910},
'SIL4':{'sat_power_rt2':170,'sat_cnts_SM_rt2':150,'sat_cnts_MM_rt2':1120},
'SIL5':{'sat_power_rt2':245,'sat_cnts_SM_rt2':190,'sat_cnts_MM_rt2':1070},
'SIL6':{'sat_power_rt2':210,'sat_cnts_SM_rt2':250,'sat_cnts_MM_rt2':1130},
'SIL7':{'sat_power_rt2':440,'sat_cnts_SM_rt2':160,'sat_cnts_MM_rt2':1010},
'SIL8':{'sat_power_rt2':40,'sat_cnts_SM_rt2':300,'sat_cnts_MM_rt2':800},
'SIL9':{'sat_power_rt2':110,'sat_cnts_SM_rt2':210,'sat_cnts_MM_rt2':880},
'SIL10':{'sat_power_rt2':315,'sat_cnts_SM_rt2':190,'sat_cnts_MM_rt2':1100},
'SIL11':{'sat_power_rt2':240,'sat_cnts_SM_rt2':70,'sat_cnts_MM_rt2':440},
'SIL12':{'sat_power_rt2':120,'sat_cnts_SM_rt2':230,'sat_cnts_MM_rt2':1030},
'SIL13':{'sat_power_rt2':0,'sat_cnts_SM_rt2':0,'sat_cnts_MM_rt2':0},
'SIL14':{'sat_power_rt2':800,'sat_cnts_SM_rt2':50,'sat_cnts_MM_rt2':790},
'SIL15':{'sat_power_rt2':125,'sat_cnts_SM_rt2':130,'sat_cnts_MM_rt2':790},
'SIL16':{'sat_power_rt2':210,'sat_cnts_SM_rt2':210,'sat_cnts_MM_rt2':1050},
}



Gretel_dict={
'SIL7':{'z':0.65	,'z_sil':0.85	,'sat_cnts_MM_rt2':0	,'sat_power_rt2':0,'sat_power_lt1':524,'sat_cnts_MM_lt1':0,'sat_cnts_SM_lt1':171},
'SIL9':{'z':0.75	,'z_sil':0.85	,'sat_cnts_MM_rt2':0	,'sat_power_rt2':0,'sat_power_lt1':0,'sat_cnts_MM_lt1':0,'sat_cnts_SM_lt1':0},
'SIL24':{'z':	0.85	,'z_sil':0.85	,'sat_cnts_MM_rt2':0	,'sat_power_rt2':0,'sat_power_lt1':0,'sat_cnts_MM_lt1':0,'sat_cnts_SM_lt1':0},
'SIL4':{'z':0.85	,'z_sil':0.85	,'sat_cnts_MM_rt2':0	,'sat_power_rt2':0,'sat_power_lt1':109,'sat_cnts_MM_lt1':1205,'sat_cnts_SM_lt1':261},
'SIL20':{'z':	0.95	,'z_sil':0.85	,'sat_cnts_MM_rt2':655	,'sat_power_rt2':127,'sat_power_lt1':290,'sat_cnts_MM_lt1':1501,'sat_cnts_SM_lt1':372},
'SIL1':{'z':1.05	,'z_sil':1.1	,'sat_cnts_MM_rt2':684	,'sat_power_rt2':34,'sat_power_lt1':49,'sat_cnts_MM_lt1':1093,'sat_cnts_SM_lt1':441	},
'SIL12':{'z':	1.05	,'z_sil':1.1	,'sat_cnts_MM_rt2':830	,'sat_power_rt2':16,'sat_power_lt1':34,'sat_cnts_MM_lt1':1018,'sat_cnts_SM_lt1':429},
'SIL6':{'z':1.15	,'z_sil':1.1	,'sat_cnts_MM_rt2':891	,'sat_power_rt2':171,'sat_power_lt1':198,'sat_cnts_MM_lt1':1184,'sat_cnts_SM_lt1':234},	
'SIL8':{'z':1.15	,'z_sil':1.1	,'sat_cnts_MM_rt2':0	,'sat_power_rt2':0,'sat_power_lt1':0,'sat_cnts_MM_lt1':0,'sat_cnts_SM_lt1':0	},
'SIL11':{'z':	1.35	,'z_sil':1.45	,'sat_cnts_MM_rt2':742	,'sat_power_rt2':41,'sat_power_lt1':58,'sat_cnts_MM_lt1':1327,'sat_cnts_SM_lt1':507	},
'SIL23':{'z':	1.45	,'z_sil':1.45	,'sat_cnts_MM_rt2':1070,'sat_power_rt2':	49,'sat_power_lt1':48,'sat_cnts_MM_lt1':1304,'sat_cnts_SM_lt1':539},
'SIL2':{'z':1.55	,'z_sil':1.45	,'sat_cnts_MM_rt2':893	,'sat_power_rt2':69,'sat_power_lt1':74,'sat_cnts_MM_lt1':1330,'sat_cnts_SM_lt1':442	},
'SIL15':{'z':	1.85	,'z_sil':2.05	,'sat_cnts_MM_rt2':1019,'sat_power_rt2':	147,'sat_power_lt1':0,'sat_cnts_MM_lt1':0,'sat_cnts_SM_lt1':0},
'SIL14':{'z':	1.95	,'z_sil':2.05	,'sat_cnts_MM_rt2':0	,'sat_power_rt2':0,'sat_power_lt1':0,'sat_cnts_MM_lt1':0,'sat_cnts_SM_lt1':0	},
'SIL3':{'z':2.15	,'z_sil':2.05	,'sat_cnts_MM_rt2':905	,'sat_power_rt2':100,'sat_power_lt1':183,'sat_cnts_MM_lt1':1473,'sat_cnts_SM_lt1':243},
'SIL10':{'z':	2.35	,'z_sil':2.5	,'sat_cnts_MM_rt2':840	,'sat_power_rt2':112,'sat_power_lt1':195,'sat_cnts_MM_lt1':1060,'sat_cnts_SM_lt1':278},
'SIL18':{'z':	2.45	,'z_sil':2.5	,'sat_cnts_MM_rt2':809	,'sat_power_rt2':180,'sat_power_lt1':0,'sat_cnts_MM_lt1':0,'sat_cnts_SM_lt1':0},
'SIL19':{'z':	2.55	,'z_sil':2.5	,'sat_cnts_MM_rt2':862	,'sat_power_rt2':103,'sat_power_lt1':0,'sat_cnts_MM_lt1':0,'sat_cnts_SM_lt1':0},
'SIL16':{'z':	2.65	,'z_sil':2.5	,'sat_cnts_MM_rt2':725,'sat_power_rt2':	65,'sat_power_lt1':107,'sat_cnts_MM_lt1':953,'sat_cnts_SM_lt1':223},
'SIL21':{'z':	2.75	,'z_sil':2.9	,'sat_cnts_MM_rt2':562	,'sat_power_rt2':309,'sat_power_lt1':228,'sat_cnts_MM_lt1':0,'sat_cnts_SM_lt1':234},
'SIL17':{'z':	2.85	,'z_sil':2.9	,'sat_cnts_MM_rt2':723	,'sat_power_rt2':174,'sat_power_lt1':156,'sat_cnts_MM_lt1':793,'sat_cnts_SM_lt1':331},
'SIL5':{'z':3.15	,'z_sil':2.9	,'sat_cnts_MM_rt2':0	,'sat_power_rt2':0,'sat_power_lt1':0,'sat_cnts_MM_lt1':0,'sat_cnts_SM_lt1':0,'sat_power_lt1':0,'sat_cnts_MM_lt1':0,'sat_cnts_SM_lt1':0}}

Fritz={
'SIL5':{'sat_power_lt3':28,'sat_cnts_SM_lt3':630,'sat_cnts_MM_lt3':1200},
'SIL1':{'sat_power_lt3':60,'sat_cnts_SM_lt3':450,'sat_cnts_MM_lt3':560},
'SIL2':{'sat_power_lt3':80,'sat_cnts_MM_lt3':615},
}

Fritz_noAR={
'SIL1':{'sat_power_rt2':130,'sat_cnts_SM_rt2':230,'sat_cnts_MM_rt2':820},
'SIL2':{'sat_power_rt2':130,'sat_cnts_SM_rt2':165,'sat_cnts_MM_rt2':719},
'SIL3':{'sat_power_rt2':600,'sat_cnts_SM_rt2':50},
'SIL5':{'sat_power_rt2':35,'sat_cnts_SM_rt2':358,'sat_cnts_MM_rt2':948},
'SIL6':{'sat_power_rt2':167,'sat_cnts_SM_rt2':167},
'SIL7':{'sat_power_rt2':0,'sat_cnts_SM_rt2':0,'sat_cnts_MM_rt2':0},
'SIL8':{'sat_power_rt2':0,'sat_cnts_SM_rt2':0,'sat_cnts_MM_rt2':0},
'SIL4':{'sat_power_rt2':537,'sat_cnts_SM_rt2':89}
}

Sam_noAR={
'SIL1':{'sat_power_rt2':130,'sat_cnts_SM_rt2':230,'sat_cnts_MM_rt2':820},
'SIL2':{'sat_power_rt2':0,'sat_cnts_SM_rt2':0,'sat_cnts_MM_rt2':0},
'SIL3':{'sat_power_rt2':0,'sat_cnts_SM_rt2':0,'sat_cnts_MM_rt2':0},
'SIL4':{'sat_power_rt2':35,'sat_cnts_SM_rt2':358,'sat_cnts_MM_rt2':948},
'SIL5':{'sat_power_rt2':35,'sat_cnts_SM_rt2':440,'sat_cnts_MM_rt2':940},
'SIL6':{'sat_power_rt2':445,'sat_cnts_SM_rt2':95,'sat_cnts_MM_rt2':1000}
}

#Gretel_bad=['SIL7','SIL14','SIL16','SIL17','SIL9','SIL24','SIL4','SIL8','SIL5','SIL18']
Gretel_bad=[]
def get_SIL_data(dict,bad_sil):
	z_SIL=[]
	sat_cnts=[]
	sat_power=[]
	z_rel=[]
	sat_cnts_MM_lt1=[]
	sat_cnts_SM_lt1=[]
	sat_ratio=[]
	sat_power_lt1=[]
	nr_of_zero=0
	for i in dict:
		if i not in bad_sil:
			z_SIL.append(7.5*dict[i]['z_sil']/2.9)
			z_rel.append(dict[i]['z']/dict[i]['z_sil'])
			sat_cnts.append(dict[i]['sat_cnts'])
			sat_cnts_MM_lt1.append(dict[i]['sat_cnts_MM_lt1'])
			sat_cnts_SM_lt1.append(dict[i]['sat_cnts_SM_lt1'])
			if (dict[i]['sat_cnts_SM_lt1'] !=0):
				sat_ratio.append(float(dict[i]['sat_cnts_MM_lt1'])/float(dict[i]['sat_cnts_SM_lt1']))
			else: 
				sat_ratio.append(0)	
				nr_of_zero+=1
			sat_power.append(dict[i]['sat_power'])
			sat_power_lt1.append(dict[i]['sat_power_lt1'])
	
	return z_SIL,z_rel,sat_cnts,sat_power,sat_power_lt1,sat_cnts_MM_lt1,sat_cnts_SM_lt1,sat_ratio,sat_power_lt1, nr_of_zero

def get_SIL_prop(dict,property):
	r=[]
	for i in dict:
		if property in dict[i]:
			r.append(dict[i][property])
	return r			

def get_setup_prop(dict_list,property):
	Prop=[]
	
	for i in dict_list:
    		p = get_SIL_prop(i,property)
    		Prop=Prop+p
	return Prop

before_dicts=[E6IIno2,the111no1_before_dict,the111no2_before_dict]
MM_rt2_before=get_setup_prop(before_dicts,'sat_cnts_MM_rt2')
MM_lt1_before=get_setup_prop(before_dicts,'sat_cnts_MM_lt1')
MM_lt2_before=get_setup_prop(before_dicts,'sat_cnts_MM_lt2')
MM_lt3_before=get_setup_prop(before_dicts,'sat_cnts_MM_lt3')
MM_all_before=MM_lt1_before+MM_lt2_before+MM_lt3_before+MM_rt2_before

after_dicts=[Fritz,Gretel_dict,Pippin_noAR,the111no1_after_dict,the111no2_after_dict] #Sam_noAR,Fritz_noAR,
MM_rt2_after=get_setup_prop(after_dicts,'sat_cnts_MM_rt2')
MM_lt1_after=get_setup_prop(after_dicts,'sat_cnts_MM_lt1')
MM_lt2_after=get_setup_prop(after_dicts,'sat_cnts_MM_lt2')
MM_lt3_after=get_setup_prop(after_dicts,'sat_cnts_MM_lt3')
MM_all_after=MM_lt1_after+MM_lt2_after+MM_lt3_after+MM_rt2_after

plt.figure(1)
plt.hist(MM_all_before,facecolor='Crimson',edgecolor='Crimson',alpha=0.75,label='Normal Trench')
#plt.hist(MM_all_after,facecolor='RoyalBlue',edgecolor='RoyalBlue',alpha=0.75,label='Enlarged Trench and AR')
#plt.legend(prop={'size':12})
#plt.hist(MM_all_before,histtype='step',edgecolor='Crimson',alpha=0.75)
#plt.hist(MM_all_after,histtype='step',edgecolor='RoyalBlue',alpha=0.75)
plt.xlabel('MM countrate [Kcnts/s]',fontsize=14)
plt.ylabel('Occurence',fontsize=14)
plt.show()

'''
z_SIL,z_rel,sat_cnts,sat_power,sat_power_lt1,sat_cnts_MM_lt1,sat_cnts_SM_lt1,sat_ratio,sat_power_lt1,nr_of_zero = get_SIL_data(Gretel_dict,Gretel_bad)
plt.figure(1)
#sat_cnts=[950,0,900,160,720,140,700,0,600,750,963,250,0,750,908,1000]
#z_SIL=arange(len(sat_cnts))
plt.plot(z_SIL,sat_cnts,'bo')
plt.xlabel ('SIL Radius [um]', fontsize = 16)
plt.ylabel ('SILnr [arb. unit]', fontsize = 16)   
plt.ylim([0,1200])
#plt.ylim ([0, 1])
plt.show()

plt.figure(2)
plt.plot(z_rel,sat_cnts,'bo')
plt.xlabel ('Conversion Factor [a.u.]', fontsize = 16)
plt.ylabel ('Saturation counts [Kcnts]', fontsize = 16)   
#plt.ylim ([0, 1])
plt.show()

plt.figure(3)
plt.plot(z_SIL,sat_power,'bo')
plt.xlabel ('SIL Radius [um]', fontsize = 16)
plt.ylabel ('Saturation power [uW]', fontsize = 16)   
#plt.ylim ([0, 1])
plt.show()

plt.figure(3)
plt.plot(z_rel,sat_power,'bo')
plt.xlabel ('Conversion Factor[a.u.]', fontsize = 16)
plt.ylabel ('Saturation power [uW]', fontsize = 16)   
#plt.ylim ([0, 1])
plt.show()

plt.figure(4)
plt.plot(sat_power,sat_ratio,'ko')
plt.xlabel ('Saturation power rt2[uW]', fontsize = 16)
plt.ylabel ('MM/SM ratio lt1', fontsize = 16)  
plt.xlim([0,300])
plt.show()

plt.figure(5)
plt.plot(sat_power_lt1,sat_cnts_SM_lt1,'bo',label='SM')
plt.plot(sat_power_lt1,sat_cnts_MM_lt1,'ro',label='MM')
#plt.plot(sat_power_lt1,sat_ratio,'ko')
plt.xlabel ('Saturation power [uW]', fontsize = 16)
plt.ylabel ('Saturation counts lt1 [Kcnts]', fontsize = 16)
plt.legend()
#plt.ylabel ('MM/SM ratio lt1', fontsize = 16)  
plt.xlim([0,300])
#plt.ylim ([0, 1])
plt.show()

plt.figure(5)
plt.plot(z_SIL,sat_ratio,'ro')
plt.xlabel ('SIL radius', fontsize = 16)
plt.ylabel ('Sat power lt1 [uW]', fontsize = 16)  
#plt.xlim([0,500])
#plt.ylim ([0, 1])
plt.show()

plt.figure(6)
plt.plot(sat_power,sat_power_lt1,'ro')
x=arange(0,500)
y=x
plt.plot(y,x,'b-')
plt.xlabel ('Sat power rt2 [uW]', fontsize = 16)
plt.ylabel ('Sat power lt1 [uW]', fontsize = 16)  
plt.xlim([10,500])
plt.ylim ([10, 500])
plt.show()

print nr_of_zero
'''

