##### Author : Liao Herui ######
##### E-mail : liaoherui@mail.dlut.edu.cn #######
import os
import re
import getopt
import sys
#############  User Info #######################################
################################################################
## How to use this pipeline	            		      ##		
##					    		      ##	 
## python 0.overall.py -l [contig.list] -s [sample_name list] ##
#############	Get  Option#####################################
opts,args=getopt.getopt(sys.argv[1:],"l:s:")
fa_list=''
now=os.getcwd()
output=now+'/overall'  #[default]
sl=''

for opt,arg in opts:
	if opt=="-l":
		fa_list=arg
	elif opt=='-s':
		sl=arg

#### sample name list ####
sn=open(sl,'r')
snl=[]
while True:
	line=sn.readline().strip()
	if not line:break
	snl.append(line)


#############	Make  the output dir and the shell dir#############
if not (os.path.exists(output)):
	os.makedirs(output,0755)


#om=open(output+'/'+'quast.sh')

############	Evaluation	##############
f=open(fa_list,'r')
om=open(output+'/quast.sh','w+')
name=''
while True:
	line=f.readline().strip()
	if not line:break
	for s in snl:
		if re.search(s,line):name=s
	postfix=re.split('/',line)[-1]
	postfix=re.split('\.',postfix)[0]
	if not os.path.exists(output+'/quast/'+name+'/'+postfix):
		os.makedirs(output+'/quast/'+name+'/'+postfix,0755)
	om.write('quast -o '+output+'/quast/'+name+'/'+postfix+' '+line+'\n')

oq=open(output+'/qsub.sh','w')
oq.write('qsub -cwd  -V -q bigmem.q -l h_vmem=15G -pe smp 4 '+output+'/quast.sh')

#os.chdir(output)
#os.system('sh qsub.sh')






	

