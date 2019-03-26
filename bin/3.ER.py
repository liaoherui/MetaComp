#### Author : Liao Herui ######
#### E-mail : liaoherui@mail.dlut.edu.cn ########
import os
import re
import gzip
import random
import getopt
import time
import sys
#######sub function###############
def motify(line,which):
	level3=''
        ele=re.split('/',line)
	if which=='checkm':
		name=ele[-5]
	        return name
	else:
		name=ele[-4]
		return name

def suiji(s,d):
        strs=[]
        strs.append(s)
        strs.append(d)
        n=random.randint(0,len(strs)-1)
        return strs[n]
##################################
##### Get Option #################
opts,args=getopt.getopt(sys.argv[1:],"o:g:c:m:n:t:")
pwd=os.getcwd()
output=pwd+'/Result'  #[default]
cl=''   #[default with coverage]
gl='Y'   #[defalut]
map_file=''	#map_rate_file list
lazy='N'
fq_type=''
for opt,arg in opts:
        if opt=="-o":
                output=arg
	elif opt=='-g':
		gl=arg
	elif opt=='-c':
		cl=arg
	elif opt=='-m':
		map_file=arg
	elif opt=='-n':
		lazy=arg
	elif opt=='-t':
if not re.search('/',output):
	tname=output
else:
	tname=re.split('/',output)[-1]
##### Main ######################
result={}
class1=[]
pwd=os.getcwd()
for filename1 in os.listdir(output):
	if re.search('_',filename1):
		class1.append(filename1)
		result[filename1]={}
for c in class1:
	for filename2 in os.listdir(output+'/'+c):
		result[c][filename2]={}
	
		
###	checkm	########
f1=open('list/checkm.list','r')
while True:
	line=f1.readline()
	if not line:break
	line=line.strip()
	c=motify(line,'checkm')
	f1_s=open(line,'r')
	while True:
		line=f1_s.readline()
		if not line:break
		line=line.strip()	
		if re.search('-',line) or re.search('Marker',line):
			#print line
			continue
		#if re.search(post,line.lower()) or re.search(post,line):
		else:
			ele=line.split()
			result[c][ele[0]]['completeness']=ele[12]	
			result[c][ele[0]]['contamination']=ele[13]
#print result
###	aragorn	#########
f2=open('list/aragorn.list','r')
while True:
	line=f2.readline()
	if not line:break
	line=line.strip()
	ele=line.split('/')
	b=ele[-3]
	c=motify(line,'aragorn')	
	f2s=open(line+'/result.txt','r')
	while True:
		line2=f2s.readline()
		if not line2:break
		if re.search('Total',line2) and re.search('tRNA',line2):
			nums=re.split('=',line2)
			num=nums[1].strip()
			result[c][b]['tRNA_Num']=num
		elif re.search('Number of tRNA genes',line2):
			nums=re.split('=',line2)
                        num=nums[1].strip()
                        result[c][b]['tRNA_Num']=num
		else:
			continue
	if 'tRNA_Num' not in result[c][b]:
		print c+' '+b+' miss its tRNA'
		result[c][b]['tRNA_Num']='NA'
###	barrnap	#########
f3=open('list/barrnap.list','r')
while True:
	line=f3.readline()
	if not line:break
	line=line.strip()
	ele=line.split('/')
	b=ele[-3]
	c=motify(line,'barrnap')
	f3_s=open(line+'/result.txt','r')
	s5=[]
	s23=[]
	s16=[]
	while True:
		line=f3_s.readline()
		if not line:break
		line=line.strip()
		if re.search('23S_rRNA',line):
			s23.append(line)
		elif re.search('16S_rRNA',line):
			s16.append(line)
		elif re.search('5S_rRNA',line):
			s5.append(line)
	if len(s5)>1 or len(s5)==1:
		result[c][b]['5S_RNA']='yes'
	else:
		result[c][b]['5S_RNA']='no'
	if len(s23)>1 or len(s23)==1:
                result[c][b]['23S_RNA']='yes'
        else:
                result[c][b]['23S_RNA']='no'
	if len(s16)>1 or len(s16)==1:
                result[c][b]['16S_RNA']='yes'
        else:
                result[c][b]['16S_RNA']='no'
####	kraken	############################		
f4=open('list/kraken.list','r')
while True:
	line=f4.readline()
	line=line.strip()
	test=line
	line=re.sub('\n','',line)
	if not line:break
	ele=line.split('/')
	b=ele[-3]
	c=motify(line,'kraken')
	f4_s=open(line+'/result.report','r')
	#species variable
	percent=[]
	ps={}
	species=[]
	cv={}
	th={}
	#genus variable
	percent2=[]
	ps2={}
	species2=[]
        cv2={}
        th2={}
	#phylum variable
	percent3=[]
        ps3={}
        species3=[]
        cv3={}
        th3={}
	while True:
		line=f4_s.readline()
		if not line:break
		ele=line.split('\t')	
		if ele[3]=='S':
			cv[float(ele[0])]=ele[0]
			if ele[5].strip() not in th:
                                th[ele[5].strip()]=int(ele[2])
                        else:
                                print ele[5]+' appears more than one time!!!'
			if ele[0] not in percent:
				percent.append(float(ele[0]))
			if ele[5] not in species:
				species.append(ele[5].strip())
			if ele[0] not in ps:
				ps[ele[0]]=ele[5].strip()
			else:
				if int(ele[2])>th[ps[ele[0]]]:
					ps[ele[0]]=ele[5].strip()
					continue
				elif int(ele[2])<th[ps[ele[0]]]:
					continue
				else:
					ps[ele[0]]=suiji(ele[5].strip(),ps[ele[0]])
		elif ele[3]=='G':
			cv2[float(ele[0])]=ele[0]
                        if ele[5].strip() not in th2:
                                th2[ele[5].strip()]=int(ele[2])
                        else:
                                print ele[5]+' appears more than one time!!!'
                        if ele[0] not in percent2:
                                percent2.append(float(ele[0]))
                        if ele[5] not in species2:
                                species2.append(ele[5].strip())
                        if ele[0] not in ps2:
                                ps2[ele[0]]=ele[5].strip()
                        else:
                                if int(ele[2])>th2[ps2[ele[0]]]:
                                        ps2[ele[0]]=ele[5].strip()
                                        continue
                                elif int(ele[2])<th2[ps2[ele[0]]]:
                                        continue
                                else:
                                        ps2[ele[0]]=suiji(ele[5].strip(),ps2[ele[0]])
		elif ele[3]=='P':
			cv3[float(ele[0])]=ele[0]
                        if ele[5].strip() not in th3:
                                th3[ele[5].strip()]=int(ele[2])
                        else:
                                print ele[5]+' appears more than one time!!!'
                        if ele[0] not in percent3:
                                percent3.append(float(ele[0]))
                        if ele[5] not in species3:
                                species3.append(ele[5].strip())
                        if ele[0] not in ps3:
                                ps3[ele[0]]=ele[5].strip()
                        else:
                                if int(ele[2])>th3[ps3[ele[0]]]:
                                        ps3[ele[0]]=ele[5].strip()
                                        continue
                                elif int(ele[2])<th3[ps3[ele[0]]]:
                                        continue
                                else:
                                        ps3[ele[0]]=suiji(ele[5].strip(),ps3[ele[0]])
		else:continue
	p=sorted(percent)
	p2=sorted(percent2)
	p3=sorted(percent3)
	if not any(p):
		result[c][b]['largest_species']='Null'
	        result[c][b]['species_num']='0'
	else:
		result[c][b]['largest_species']=ps[cv[p[-1]]]+'('+cv[p[-1]].strip()+')'
	        result[c][b]['species_num']=str(len(species))
	if not any(p2):
		result[c][b]['largest_genus']='Null'
	        result[c][b]['genus_num']='0'
	else:
		result[c][b]['largest_genus']=ps2[cv2[p2[-1]]]+'('+cv2[p2[-1]].strip()+')'
	        result[c][b]['genus_num']=str(len(species2))
	if not any(p3):
		result[c][b]['largest_phylum']='Null'
	        result[c][b]['phylum_num']='0'
	else:
		result[c][b]['largest_phylum']=ps3[cv3[p3[-1]]]+'('+cv3[p3[-1]].strip()+')'
	        result[c][b]['phylum_num']=str(len(species3))

####	quast	#############################	
f5=open('list/quast.list','r')
while True:
	line=f5.readline()
	if not line:break
	line=line.strip()
	ele=line.split('/')
	b=ele[-3]
	c=motify(line,'quast')
	f5_s=open(line+'/report.txt','r')	
	while True:
		line=f5_s.readline()
		line2=line
		#line=line.strip()
		if not line:break	
		if re.search('N50',line):
			line=line.split()
			n50=line[1]
			result[c][b]['N50']=n50
		if re.search('Total',line2) and not re.search('\(',line2):
			line4=line2.split()
			tl=line4[2]
			result[c][b]['total_length']=tl

####	coverage	####################
if  not cl=='':
	f6=open(cl,'r')
	#line=f6.readline()
	while True:
		line=f6.readline().strip()
		if not line:break
		f_cov=open(line,'r')
		if fq_type=='PE':
			line_cov=f_cov.readline()
			while True:
				line_cov=f_cov.readline().strip()
				if not line_cov:break
				ele=line_cov.split('\t')
				c=re.split('\.',ele[0])[0]
				#c=re.sub('\.fasta','',c)
				b=re.sub('\.fasta','',ele[0])
				cov=float(ele[1])+float(ele[2])
				if c not in result:
					print c+' of '+line+' is not in result dict! Please check.'
					exit()
				result[c][b]['coverage']=cov
		else:
			#SE fastq file ->  cal the bin coverage by formula
			contig={}
			bin_dir=re.split('/',line)[:-1]
			bin_dir='/'.join(bin_dir)
			while True:
				line_cov=f_cov.readline().strip()
				if not line_cov:break
				ele=re.split('\t',line_cov)	
				if ele[0] not in contig:
					contig[ele[0]]=float(ele[1])
			for filename in os.listdir(bin_dir):
				if not re.search('\.fasta',filename):continue
				c=re.split('\.',filename)[0]
				b=re.sub('\.fasta','',filename)
				fb=open(bin_dir+'/'+filename,'r')
				counter=1
				len_count=0
				contig_length_cov=0
				bin_length=0
				while True:
					line=fb.readline.strip()
					if not line:break	
					if re.search('>',line):
						if counter==1:
							contig_name=re.sub('>','',line)	
							counter+=1
						else:
							contig_length_cov+=contig[contig_name]*len_count
							contig_name=re.sub('>','',line)
							len_count=0	
							
					else:
						len_count+=len(line)
						bin_length+=len(line)
				contig_length_cov+=contig[contig_name]*len_count
				bin_length+=len(line)
				if not bin_length==float(result[c][b]['total_length']):
					print 'Warning : The Bin : '+b+' length has 2 values : '+str(bin_length)+'(Manual)  '+str(result[c][b]['total_length'])+'(quast) '
				bin_cov=float(contig_length)/float(bin_length)
				result[c][b]['coverage']=bin_cov
								
						
######## Bin relative abundance ##########
if not map_file=='':
	f7=open(map_file,'r')
	num_len={}
	while True:
		line=f7.readline().strip()
		if not line:break
		ele=line.split('\t')
		if ele[0] not in result:
			print ele[0]+' not in your dict ! Please check!'
			exit()
		num_len[ele[0]]={}
		num_len[ele[0]]['reads_len']=ele[2]
		fm=open(ele[1],'r')	
		mline=fm.readline().strip()
		reads_num=re.split('\+',mline)
		reads_num=reads_num[0].strip()
		num_len[ele[0]]['reads_num']=reads_num
	for key1 in result:
		for key2 in result[key1]:
			if key1 not in num_len:
				result[key1][key2]['abundance']='NA'
			else:
				result[key1][key2]['abundance']=(float(result[key1][key2]['coverage'])*float(result[key1][key2]['total_length']))/(float(num_len[key1]['reads_len'])*float(num_len[key1]['reads_num']))
	
####	quality #########
for key1 in result:
	for key2 in result[key1]:
		if float(result[key1][key2]['completeness'])>90.0 and float(result[key1][key2]['contamination'])<5.0 and (int(result[key1][key2]['tRNA_Num'])>18 or int(result[key1][key2]['tRNA_Num'])==18 ) and result[key1][key2]['5S_RNA']=='yes' and result[key1][key2]['16S_RNA']=='yes' and result[key1][key2]['23S_RNA']=='yes':
			result[key1][key2]['quality']='High'
		elif (float(result[key1][key2]['completeness'])>50.0 or float(result[key1][key2]['completeness'])==50.0) and float(result[key1][key2]['contamination'])<10.0:
			result[key1][key2]['quality']='Medium'
		elif float(result[key1][key2]['completeness'])<50.0 and float(result[key1][key2]['contamination'])<10.0:
                        result[key1][key2]['quality']='Low'
		else:	result[key1][key2]['quality']='Other'

####	Make Output File	#########
if not cl=='':
	feature=['completeness','contamination','tRNA_Num','5S_RNA','16S_RNA','23S_RNA','largest_species','species_num','largest_genus','genus_num','largest_phylum','phylum_num','N50','total_length','coverage','abundance','quality']
else:
	feature=['completeness','contamination','tRNA_Num','5S_RNA','16S_RNA','23S_RNA','largest_species','species_num','largest_genus','genus_num','largest_phylum','phylum_num','N50','total_length','quality']
if not os.path.exists('Tabel_'+tname):
                os.makedirs('Tabel_'+tname,0755)
for key1 in result:
	#if not os.path.exists('Result/'+key1):
		#os.makedirs('Result/'+key1,0755)
	o=open('Tabel_'+tname+'/'+key1+'.bins.result','w+')
	if not cl=='':
		o.write('bins_name\tcompleteness\tcontamination\ttRNA_Num\t5S_RNA\t16S_RNA\t23S_RNA\tlargest_species\tspecies_num\tlargest_genus\tgenus_num\tlargest_phylum\tphylum_num\tN50\ttotal_length\tcoverage\tabundance\tquality\n')
	else:
		o.write('bins_name\tcompleteness\tcontamination\ttRNA_Num\t5S_RNA\t16S_RNA\t23S_RNA\tlargest_species\tspecies_num\tlargest_genus\tgenus_num\tlargest_phylum\tphylum_num\tN50\ttotal_length\tquality\n')
	for key2 in result[key1]:
		o.write(key2)
		for fe in feature:
			if fe not in result[key1][key2]:
				print key1+' '+key2+' is wrong\n'+fe+' is missing!'
				exit()
			o.write('\t'+str(result[key1][key2][fe]))
		o.write('\n')
