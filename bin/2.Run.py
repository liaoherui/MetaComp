import os
import re
import getopt
import sys
import time

opts,args=getopt.getopt(sys.argv[1:],"o:a:")
pwd=os.getcwd()
output=pwd+'/Result'  #[default]
run='Y'  #[default]
gl='N'   #[default]
all_big='N'  #[default]
for opt,arg in opts:
        if opt=="-o":
                output=arg
	elif opt=='-a':
		all_big=arg
if not re.search('/',output):
	output=pwd+'/'+output
	
if not os.path.exists('Submit'):
        os.makedirs('Submit',0755)
oz=open('Submit/qsub.sh','w+')
if run=='Y':
	for filename in os.listdir('Merge_shell'):
		os.chdir(pwd+'/Merge_shell')
		if re.search('checkm',filename):
			name=re.sub('.sh','',filename)
			name='MAEP_'+name
			oz.write('qsub -cwd  -V -q bigmem.q -l h_vmem=100G -pe smp 10 -N '+name+' '+pwd+'/Merge_shell/'+filename+'\n')
			#os.system('qsub -cwd  -V -q bigmem.q -l h_vmem=100G -pe smp 10 -N '+name+' '+pwd+'/Merge_shell/'+filename)
		else:
			name=re.sub('.sh','',filename)
			name='MAEP_'+name
			if all_big=='N':
				#os.system('qsub -cwd  -V -q beta.q -l h_vmem=80G -pe smp 10 -N '+name+' '+pwd+'/Merge_shell/'+filename)
				oz.write('qsub -cwd  -V -q beta.q -l h_vmem=80G -pe smp 10 -N '+name+' '+pwd+'/Merge_shell/'+filename+'\n')
			else:
				oz.write('qsub -cwd  -V -q bigmem.q -l h_vmem=80G -pe smp 10 -N '+name+' '+pwd+'/Merge_shell/'+filename+'\n')
				#os.system('qsub -cwd  -V -q bigmem.q -l h_vmem=80G -pe smp 10 -N '+name+' '+pwd+'/Merge_shell/'+filename)

