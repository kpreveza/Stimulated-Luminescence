		#########################################
		#                                       #
		#   Deconvolution of LM-OSL curves      #
		#                                       #
		#   Author: Prevezanou Konstantina      #
		#                                       #
		#   E-mail: kpreveza@physics.auth.gr    #
		#                                       #
		#########################################

#LIBRARIES

'''We import the libraries that we are going to use'''
import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.special import lambertw
from lmfit import Parameters, report_fit, Parameter, Model
import pybroom as br
import pandas as pd
from easygui import *
from tkinter import *
from tkinter import filedialog,ttk
import os
from os.path import exists


#FUNCTIONS

def lm(x,Im,tm,R):
	c=(1-R)/R
	z=1/c-np.log(c)+((x**2)/(tm**2*(1-R)*(1+0.534156*R**0.7917)))
	z_m=1/c-np.log(c)+(1/((1-R)*(1+0.534156*R**0.7917)))
	#w=np.real(lambertw(np.exp(z)))
	w_m=np.real(lambertw(np.exp(z_m)))

	for i in range(len(z)):
		if z[i]<710:
			w=np.real(lambertw(np.exp(z)))
		else:
			w=z-np.log(z)
		return x*(Im/tm)*((w_m+w_m**2)/(w+w**2))


def make_model(num):
	pref = "lm{}_".format(num)
	peak = list_peaks[num]
	model = Model(lm, prefix = pref)
	model.set_param_hint(pref+'Im', value = a[0][num], min = a[1][num], max = a[2][num])
	model.set_param_hint(pref+'tm', value = b[0][num], min = b[1][num], max = b[2][num])
	model.set_param_hint(pref+'R', value = c[0][num], min = c[1][num], max = c[2][num])
	return model


def loop_init(name1,name2,name3,num):
	'''Function to make the list for each parameter including the initial parameters'''
	name1=[]
	name2=[]
	name3=[]
	for i in range(num,len(datap),4):
		name1.append(datap[i][0])
		name2.append(datap[i][1])
		name3.append(datap[i][2])
	name1=np.array(name1,dtype=float)
	name2=np.array(name2,dtype=float)
	name3=np.array(name3,dtype=float)
	return name1,name2,name3


def file_exist(f):
	'''Function to check if a file exists and if it does it adds a number at the end of the name'''
	fnew = f
	root, ext = os.path.splitext(f)
	i = 0
	while os.path.exists(fnew):
		i += 1
		fnew = '%s_%i%s' % (root, i, ext)
	return fnew


#INPUT DATA

'''Info and input from user about experimental data, we use the tkinter method to make a pop-up box'''
root=Tk()
root.title('Data for the deconvolution')
root.geometry('700x370')

global b_sequence, files, fileswithoutpath
b_sequence=[]
files=[]
fileswithoutpath=[]


def addfile():
	'''Function to go through computer files'''
	filename=filedialog.askopenfilename(initialdir=os.getcwd(),title="Select File",filetypes=(('TXT files','.txt'),('All files','.')))
	files.append(filename)
	filenamewithout=os.path.split(filename)[1]
	fileswithoutpath.append(filenamewithout)

def which_button(t):
	'''Function to save the order of button push'''
	b_sequence.append(t)

def clickonce(k):
	'''Function to allow the user to press the button only once'''
	k.configure(state=DISABLED)

notice=Label(root,text="What you should know to run the script properly:",font=('Arial',10, 'underline'))
notice.pack()

info=Label(root,text="\n1. You need one file with experimental data (first column index in ascending order, second column x-axis data,\nthird column y-axis data) in ‘txt (tab delimited)’ formating.\n\n2. You need a second file with the initial guess of the parameters in ‘txt (tab delimited)’ formating.\nThe format should be: three columns - first contains the values, second the minimum and third the maximum for each value.\nAlso, for each peak the parameters (Im, t, R) will be given and between the different peaks should be a blank line.\n\n3. Read carefully the instructions of each pop-up box to correctly complete the gaps.\n\n4. At the end, fitted parameters will be stored in a new directory, ‘Fitted parameters’,\nand theoretical data in new directory, ‘Theoretical data’.\n\n5. To save the graph just click on ‘Save’ (the button is on the down left corner of the graph).\n\n")
info.pack()

out=Button(root,text='Continue',padx=5,pady=4,fg='white',bg='black',command=root.destroy ,font=('Arial', 11))
out.pack(side=BOTTOM)

frame1=Frame(root)
frame1.place(relwidth=0.41,relheight=0.1,relx=0.3,rely=0.73)

expdata=Button(frame1,text='Experimental data',padx=6,pady=5,fg='white',bg='black',command=lambda:[addfile(),which_button('exp'),clickonce(expdata)],font=('Arial', 11))
expdata.pack(side=LEFT)

expparams=Button(frame1,text='Initial parameters',padx=6,pady=5,fg='white',bg='black',command=lambda:[addfile(),which_button('para'),clickonce(expparams)],font=('Arial', 11))
expparams.pack(side=RIGHT)

root.mainloop()

'''We import the experimental data'''
if b_sequence[0]=='exp':
	#Import experimental data
	with open (files[0], 'r') as file:
		data=list(csv.reader(file,delimiter='\t'))

	usedata=np.array(data[1:], dtype=float)
	x=usedata[:,1]
	y=usedata[:,2]


	#Import parameters
	with open (files[1], 'r') as filep:
		datap=list(csv.reader(filep,delimiter='\t'))
		
	a = loop_init('Iinit','minIm','maxIm',1)
	b = loop_init('tinit','mintm','maxtm',2)
	c = loop_init('Rinit','minR','maxR',3)

	'''We create the names of the output files'''
	file1=fileswithoutpath[0].replace('.txt','')

	newfile1=file1+'_out.txt'
	newfile2=file1+'_param.txt'
else:
	#Import experimental data
	with open (files[1], 'r') as file:
		data=list(csv.reader(file,delimiter='\t'))

	usedata=np.array(data[1:], dtype=float)
	x=usedata[:,1]
	y=usedata[:,2]


	#Import parameters
	with open (files[0], 'r') as filep:
		datap=list(csv.reader(filep,delimiter='\t'))
	
	a = loop_init('Iinit','minIm','maxIm',1)
	b = loop_init('tinit','mintm','maxtm',2)
	c = loop_init('Rinit','minR','maxR',3)

	'''We create the names of the output files'''
	file1=fileswithoutpath[1].replace('.txt','')

	newfile1=file1+'_out.txt'
	newfile2=file1+'_param.txt'


'''User tells the number of peaks, the start+end of the deconvolution and the heating rate for the calculation of the frequency factor'''
msg = "Write down how many peaks you want to use, the start and end number."
title = "Data for deconvolution"
fieldNames = ['Number of peaks',"Start number",'End number']
fieldValues = []
fieldValues = multenterbox(msg,title,fieldNames)

peaks=int(fieldValues[0])
list_peaks=range(peaks)
startfit = int(fieldValues[1])
endfit = int(fieldValues[2])

'''We define the start and end of the experimental data that we are going to use'''
x=np.array(x[startfit:endfit])
y=np.array(y[startfit:endfit])


#DECONVOLUTION

'''Run the model to make the total peak'''
mod = None
for i in range(len(list_peaks)):
	this_mod = make_model(i)
	if mod is None:
		mod = this_mod
	else:
		mod = mod + this_mod


'''Choose method for deconvolution'''
msg = 'Choose one method to perform the deconvolution.\n\n1. Levenberg Marquardt - leastsq\n2. Nelder Mead - nelder\n3. Powell - powell'
title = "Deconvolution method"
choices = ['leastsq','nelder','powell']
choice = choicebox(msg, title, choices)

howtofit=choice

'''Deconvolution'''
out=mod.fit(y, x=x, method=howtofit,nan_policy='omit')


#NEW DIRECTORIES

'''Creating new directories to save separately the output files'''
if not os.path.exists('Fitted parameters'):
	os.mkdir('Fitted parameters') 
else:    
	pass

if not os.path.exists('Theoretical data'):
	os.mkdir('Theoretical data')
else:    
	pass


#FILES

'''File with theoretical data for each and total peak'''
da = br.augment(out)
da=da.round(3)
da.to_csv(file_exist("Theoretical data/{}".format(newfile1)),sep='\t')


'''Calculation of total FOM'''
residual = da['residual'].to_numpy()
yfitdata= da['best_fit'].to_numpy()
fom=round(100*(np.sum(abs(residual))/(np.sum(yfitdata))),3)
Fom=['' for i in range(peaks)]
Fom[0]=fom


'''Parameters file creation'''
dt=br.tidy(out)
values=dt['value'].to_numpy()

If=[]
tf=[]
Rf=[]
name_peaks=[]

for i in range(1,peaks+1):
	name_peaks.append('Peak_{}'.format(i))
for i in range(0,len(dt),3):
	If.append(values[i])
for i in range(1,len(dt),3):
	Rf.append(values[i])
for i in range(2,len(dt),3):
	tf.append(values[i])


df = pd.DataFrame(list(zip(name_peaks,If,tf,Rf,Fom)),columns =['Peaks', 'Imax','tmax','R','Total Fom'])
df=df.round(3)
df.to_csv(file_exist("Fitted parameters/{}".format(newfile2)),sep='\t')

'''New dataframe with fitted parameters thw will be used for the separate peak files'''
dk=pd.DataFrame(list(zip(name_peaks,If,tf,Rf)),columns =['Peaks', 'Imax','tmax','R'])
dk=dk.round(3)


'''Peaks to separate files'''
msg = 'Do you want to save each peak in a seperate file?'
title = 'Seperate files'
choices = ['Yes','No']
choice = choicebox(msg, title, choices)

seper_file=choice

if seper_file=='No':
	pass
else:
	'''Extra condition to import'''
	msg = "What is the experimental condition (ex. Dose, Temperature) ?"
	title = "Condition"
	fieldNames = ['Value']
	fieldValues = []
	fieldValues = multenterbox(msg,title,fieldNames)

	column_vl=fieldValues[0]

	dk['Condition'] = column_vl

	'''Separate files creation'''
	for i in range(peaks):
		if exists('Fitted parameters/Peak_{}_{}.txt'.format(i+1,file1))==True:
			with open ('Fitted parameters/Peak_{}_{}.txt'.format(i+1,file1), 'a') as filep:
				filep.write('\n')
				filep.write('{}'.format(i)+'  '+str(dk.loc[i,'Peaks'])+'  '+str(dk.loc[i,'Imax'])+'  '+str(dk.loc[i,'tmax'])+'  '+str(dk.loc[i,'R'])+'  '+str(dk.loc[i,'Condition']))
		else:
			with open ('Fitted parameters/Peak_{}_{}.txt'.format(i+1,file1), 'w') as filep:
				filep.write(str(dk.loc[[i]]))


#OUTPUT

'''Print fitted parameters and FOM'''
print('\nTotal FOM is: ', fom, ' %.\n')
out.params.pretty_print(columns=['value', 'min', 'max','init_value'])


'''Plot each peak'''
i=0
for i in range(peaks):
	ytl=da["Model(lm, prefix='lm{}_')".format(i)].to_numpy()
	plt.plot(x,ytl,label='Peak {}'.format(i+1))
	

'''Plot theoretical and experimental total peak'''
plt.plot(x, out.best_fit,label='Theoretical Data, \nFOM = {}'.format(fom),color='g')
plt.scatter(x,y,label='Experimental Data',color='r',marker='.')
plt.title('LM-OSL Data')
plt.xlabel(data[0][1])
plt.ylabel(data[0][2])
plt.xscale('log')
plt.grid(True)
plt.tight_layout()
#plt.legend()
plt.show()

#plt.savefig('G:/My Drive/Paper LED 2021/LM OSL/lmosl_data.eps', format='eps')