### Author : Liao Herui ###
### Version : MAEP 3.0 ###
import re
import os
import getopt
import sys
###### Get opt ########
opts,args=getopt.getopt(sys.argv[1:],"l:o:")
contig_raw_dir=''
#sample_dir=''
output_dir=''

for opt,arg in opts:
	if opt=='-l':
		contig_raw_dir=arg
	elif opt=='-o':
		output_dir=arg

if output_dir=='' or contig_raw_dir=='' :
	print 'Error:\n\tYou require to give the complete parameters!\n\tPlease check:\n\t\t -l [contig and raw reads list]\n\t\t -s [sample name list]\n\t\t -o [output dir] '
	exit()
### Current working dir ###
pwd=os.getcwd()

if not re.search('/',output_dir):
	#pwd=os.getcwd()
	output_dir=pwd+'/'+output_dir
main_bash_dir=re.split('/',output_dir)[-1]
#### initialize data structure ####
#1.sample name array
sample_arr=[]

#2.contig and reads dict
main={}  #this is the main dict including pre / contig /raw reads information
fl=open(contig_raw_dir,'r')
while True:
	line=fl.readline().strip()
	if not line:break
	ele=line.split('\t')
	if ele[0] not in sample_arr:
		sample_arr.append(ele[0])
	if ele[0] not in main:
		main[ele[0]]={}
		main[ele[0]][ele[1]]={}
		main[ele[0]][ele[1]]['contig']=ele[2]
		if len(ele)==4:
			main[ele[0]][ele[1]]['type']='SE'
			main[ele[0]][ele[1]]['raw_reads']=ele[3]
		elif len(ele)==5:
			main[ele[0]][ele[1]]['type']='PE'
			main[ele[0]][ele[1]]['raw_reads1']=ele[3]
			main[ele[0]][ele[1]]['raw_reads2']=ele[4]
		else:
			print 'Error:\n\tYour contig and raw reads list is wrong!\n\tEach line should include 4(SE) or 5(PE) elements.However,your list file doesn\'t meet the conditions.\n\tPlease check,correct and run again!\n\tThe error happens in the below line:\n\t\t'+line+'\n'
			exit()
	else:
		if ele[1] in main[ele[0]]:
			print 'Error:\n\tYou have used the same strategy with one sample!\n\tPlease check!\n\tThe error happens in the below line:\n\t\t'+line+'\n'
			exit()
		main[ele[0]][ele[1]]={}
                main[ele[0]][ele[1]]['contig']=ele[2]
                if len(ele)==4:
                        main[ele[0]][ele[1]]['type']='SE'
                        main[ele[0]][ele[1]]['raw_reads']=ele[3]
                elif len(ele)==5:
                        main[ele[0]][ele[1]]['type']='PE'
                        main[ele[0]][ele[1]]['raw_reads1']=ele[3]
			main[ele[0]][ele[1]]['raw_reads2']=ele[4]
                else:
                        print 'Error:\n\tYour contig and raw reads list is wrong!\n\tEach line should include 4(SE) or 5(PE) elements.However,your list file doesn\'t meet the conditions.\n\tPlease check,correct and run again!\n\tThe error happens in the below line:\n\t\t'+line+'\n'
                        exit()
#test or see main structure
#print main
#exit()

#### ** Pipeline: Go ! ** ####

#### ** Build the merge bash script ** ####
lazy_o=open('Auto_Lazy_bash.sh','w+')
#### ** Make  dir and bin/coverage/map list generate and Main bash Script generate ** #####
step_arr=['1.Contig_Quast','2.Binning_Maxbin','3.Map_BWA','4.Evaluate_MAEP','5.Table','6.Plot','7.Report','8.Log_Tem_File','9.Shell']
for arr in sample_arr:
	if arr not in main:
		print 'Error:\n\tSample : '+arr+' is not in your contig and raw reads list.\n\tPlease check!'
	for srr in step_arr:
		if not os.path.exists(output_dir+'/'+arr+'/'+srr):
			os.makedirs(output_dir+'/'+arr+'/'+srr,0755)
	#Bin parameter generate and Shell of 1,2,3 generte
	lazy_o.write('cd '+output_dir+'/'+arr+'/9.Shell \n')
	lazy_o.write('qsub -cwd  -V -q beta.q -l h_vmem=15G -pe smp 2 shell_quast1.sh \n')
	#obt=open(output_dir+'/'+arr+'/8.Log_Tem_File/bin.list','w+')
	# ** List (Coverage/Map_rate) Generate ** #
	cov=open(output_dir+'/'+arr+'/8.Log_Tem_File/coverage.list','w+')
	map_rate=open(output_dir+'/'+arr+'/8.Log_Tem_File/map_rate.list','w+')
	oqt=open(output_dir+'/'+arr+'/9.Shell/shell_quast1.sh','w+')
        oqt.write('quast -o '+output_dir+'/'+arr+'/1.Contig_Quast')
	for key in sorted(main[arr].keys()):
		if not os.path.exists(output_dir+'/'+arr+'/2.Binning_Maxbin/'+arr+'_'+key):
			os.makedirs(output_dir+'/'+arr+'/2.Binning_Maxbin/'+arr+'_'+key,0755)
		#obt.write(arr+'_'+key+'\t'+output_dir+'/'+arr+'/2.Binning_Maxbin/'+arr+'_'+key+'\n')
		#map_rate info and coverage list #
		if main[arr][key]['type']=='PE':
			cov.write(output_dir+'/'+arr+'/2.Binning_Maxbin/'+arr+'_'+key+'/'+arr+'_'+key+'.abundance\n')
		else:
			cov.write(output_dir+'/'+arr+'/2.Binning_Maxbin/'+arr+'_'+key+'/'+arr+'_'+key+'.abund1\n')
		map_rate.write(arr+'_'+key+'\t'+output_dir+'/'+arr+'/3.Map_BWA/'+arr+'_'+key+'_map_rate.info\t150\n')
		## Main bash script generate ##
		if re.match('\d',key):
			name=arr+'_'+key
			ost=open(output_dir+'/'+arr+'/9.Shell/'+name+'_shell.sh','w+')
			lazy_o.write('qsub -cwd  -V -q bigmem.q -l h_vmem=100G -pe smp 32 '+name+'_shell.sh\n')
		else:
			ost=open(output_dir+'/'+arr+'/9.Shell/'+key+'_shell.sh','w+')
			lazy_o.write('qsub -cwd  -V -q bigmem.q -l h_vmem=100G -pe smp 32 '+key+'_shell.sh\n')
		ost.write('#### ** Step 2: Binning_Maxbin ** ####\n')
		ost.write('echo ==========['+arr+'_'+key+': Binning_Maxbin ...] Start at : `date` ==========\n')
		if main[arr][key]['type']=='PE':
			ost.write('perl /mnt/osf2/user/liaoherui/Maxbin/MaxBin-2.2.5/run_MaxBin.pl -contig '+main[arr][key]['contig']+' -reads '+main[arr][key]['raw_reads1']+' -reads2 '+main[arr][key]['raw_reads2']+' -out '+output_dir+'/'+arr+'/2.Binning_Maxbin/'+arr+'_'+key+'/'+arr+'_'+key+' -thread 32 &&\\\n')
		else:
			ost.write('perl /mnt/osf2/user/liaoherui/Maxbin/MaxBin-2.2.5/run_MaxBin.pl -contig '+main[arr][key]['contig']+' -reads '+main[arr][key]['raw_reads']+' -out '+output_dir+'/'+arr+'/2.Binning_Maxbin/'+arr+'_'+key+'/'+arr+'_'+key+' -thread 32 &&\\\n')
		ost.write('echo ==========['+arr+'_'+key+':  Binning_Maxbin ...] End at : `date` ==========\n\n')
		ost.write('#### ** Step 3: Map_Bwa ** ####\n')
		ost.write('echo ==========['+arr+'_'+key+': Map_Bwa ...] Start at : `date` ==========\n')
		ost.write('cp '+main[arr][key]['contig']+' '+output_dir+'/'+arr+'/3.Map_BWA  &&\\\n')	
		fa_name=re.split('/',main[arr][key]['contig'])[-1]
		ost.write('bwa index '+output_dir+'/'+arr+'/3.Map_BWA/'+fa_name+' &&\\\n')
		if main[arr][key]['type']=='PE':
			ost.write('bwa mem -t 32 -m 100G '+output_dir+'/'+arr+'/3.Map_BWA/'+fa_name+' '+main[arr][key]['raw_reads1']+' '+main[arr][key]['raw_reads2']+' | samtools sort -m 50G -o '+output_dir+'/'+arr+'/3.Map_BWA/'+arr+'_'+key+'.bam - &&\\\nsamtools flagstat '+output_dir+'/'+arr+'/3.Map_BWA/'+arr+'_'+key+'.bam > '+output_dir+'/'+arr+'/3.Map_BWA/'+arr+'_'+key+'_map_rate.info  &&\\\n')
		else:
			ost.write('bwa mem -t 32 -m 100G '+output_dir+'/'+arr+'/3.Map_BWA/'+fa_name+' '+main[arr][key]['raw_reads']+' | samtools sort -m 50G -o '+output_dir+'/'+arr+'/3.Map_BWA/'+arr+'_'+key+'.bam - &&\\\nsamtools flagstat '+output_dir+'/'+arr+'/3.Map_BWA/'+arr+'_'+key+'.bam > '+output_dir+'/'+arr+'/3.Map_BWA/'+arr+'_'+key+'_map_rate.info  &&\\\n')
		ost.write('echo ==========['+arr+'_'+key+': Map_Bwa ...] End at : `date` ==========\n\n')
		ost.write('#### ** Step 4: Evaluate_MAEP ** ####\n')
		ost.write('echo ==========['+arr+'_'+key+': Evaluate_MAEP ...] Start at : `date` ==========\n')
		ost.write('cd '+output_dir+'/'+arr+' \n')
		ost.write('python '+pwd+'/bin/1.AR.py -l \''+arr+'_'+key+'\t'+output_dir+'/'+arr+'/2.Binning_Maxbin/'+arr+'_'+key+'\' -o '+output_dir+'/'+arr+'/4.Evaluate_MAEP &&\\\n')
		ost.write('sh '+output_dir+'/'+arr+'/MAEP_Shell/all_checkm_'+arr+'_'+key+'.sh &&\\\n')
		ost.write('sh '+output_dir+'/'+arr+'/MAEP_Shell/'+arr+'_'+key+'.sh  &&\\\n')
		#judge if the last step is all over
		ost.write('#### ** Judge if the last step is over ** ####\n')
		file_num=len(main[arr].keys())
                ost.write('s=`ls -l '+output_dir+'/'+arr+'/4.Evaluate_MAEP |grep "^d"|wc -l'+'`\n')
		ost.write('if [ \"$s\" -ne '+str(file_num)+' ] ; then\n\techo \"'+arr+'\t'+key+':$s'+'\"\n\texit 1\nfi\n\n')
		ost.write('python '+pwd+'/bin/2.5.generate_list.py -o '+output_dir+'/'+arr+'/4.Evaluate_MAEP'+' -g Y &&\\\n')
		ost.write('echo ==========['+arr+'_'+key+': Evaluate_MAEP ...] End at : `date` ==========\n\n')
		ost.write('#### ** Step 5: Tabel_Plot ** ####\n')
		ost.write('echo ==========['+arr+'_'+key+': Tabel_Plot ...] Start at : `date` ==========\n')
		if main[arr][key]['type']=='PE':
			ost.write('python '+pwd+'/bin/3.ER.py -o 4.Evaluate_MAEP '+' -g Y '+' -c '+output_dir+'/'+arr+'/8.Log_Tem_File/coverage.list'+' -m '+output_dir+'/'+arr+'/8.Log_Tem_File/map_rate.list -n Y -t '+main[arr][key]['type']+' &&\\\n')
		else:
			ost.write('python '+pwd+'/bin/3.ER.py -o 4.Evaluate_MAEP '+' -g Y '+' -c '+output_dir+'/'+arr+'/8.Log_Tem_File/coverage.list'+' -m '+output_dir+'/'+arr+'/8.Log_Tem_File/map_rate.list -n Y -t '+main[arr][key]['type']+' &&\\\n')
		ost.write('python '+pwd+'/bin/4.Plot.py -v N '+' -s '+arr+' -o '+arr+' &&\\\n')
		ost.write('mv Tabel_4.Evaluate_MAEP/* 5.Table \n')		
		ost.write('#mv Plot/'+arr+'/* 6.Plot \n')
		ost.write('#mv Report_'+arr+' 7.Report\n')
		ost.write('echo ==========['+arr+'_'+key+': Tabel_Plot ...] End at : `date` ==========\n')
		oqt.write(' '+main[arr][key]['contig'])
	
	
#Module test	
#exit()







	

	
	
	




