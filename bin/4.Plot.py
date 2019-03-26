#coding=utf-8
from __future__ import unicode_literals
##### Author : Liao Herui ######
##### E-mail : liaoherui@mail.dlut.edu.cn #######
import os
import re
import getopt
import sys
import rpy2.robjects as robjects
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as colors
import matplotlib.cm as cmx
plt.switch_backend('agg')
from pyecharts import Scatter,Overlap,Grid,Page,Line,Bar
#############  User Info #############################
######################################################
## How to use this pipeline	            	    ##		
##					    	    ##	 
## python 4.Plot.py -v [overall] -s [sample_name_list]##
#############  Get  Option ###########################
opts,args=getopt.getopt(sys.argv[1:],"v:s:o:")
#fa_list=''
now=os.getcwd()
output=now+'/overall'  #[default]
overall='Y'
sl=''
ero=''   #Refer to table name here.ER output -> ero

for opt,arg in opts:
	if opt=="-v":
		overall=arg
	elif opt=='-s':
		sl=arg
	elif opt=='-o':
		ero=arg
if re.search('/',ero):
	ero=re.split('/',ero)[-1]
	

if sl=='':
	print 'Please give the conf file sample_name list! This parameter is required!'
	exit()
print 'Plot Start ... ... ...'
#### Sample Name  List ######
if re.search('/',sl):
	fsn=open(sl,'r')
	snl=[]
	while True:
		line=fsn.readline().strip()
		if not line:break
		snl.append(line)
else:
	snl=[]
	snl.append(sl)

#############	Make  the output dir and the shell dir#############
for s in snl:
	if not os.path.exists('Plot/'+s):
		os.makedirs('Plot/'+s,0755)
if not os.path.exists('Log'):
	os.makedirs('Log',0755)
###### Over All N50 Value ##############################
if overall=='Y':
	#n50_all=[]
	#contig_name=[]
	if len(os.listdir('overall/quast'))<1:
		print 'Something wrong in overall quast...Please check...'
		exit()
	for filename1 in os.listdir('overall/quast'):
		n50_all=[]
	        contig_name=[]
		for filename in os.listdir('overall/quast/'+filename1):
			contig_name.append(filename)
			fo=open('overall/quast/'+filename1+'/'+filename+'/'+'report.txt','r')
			while True:
				line=fo.readline()
				if not line:break
				if re.search('N50',line):
		                        line=line.split()
	        	                n50=float(line[1])/float(1000)
					n50_all.append(n50)
		x=range(1,len(contig_name)+1)
		plt.figure(figsize=(12,6))
		plt.bar(x,n50_all,width=0.35,align='center',color='c',alpha=0.8)
		plt.xticks(x,contig_name,size='small',rotation=20)
		for a,b in zip(x,n50_all):
			plt.text(a, b+0.05, '%.0f' % b, ha='center', va= 'bottom',fontsize=7)
		plt.savefig('Plot/'+filename1+'/0.n50_overall.png')

####### Major Hash Generate #######
##### Start... ##########
'''
####for coverage, this hash requires input list#########
co=open('Log/coverage.hash','r')
a=co.read()
cov=eval(a)
'''


bq={}   #means bin quality ,eg: bq['zxy_pacbio_HGA']['High']=10 /bq['zxy'][zxy_pacbio_HGA]['completeness']['zxy_pacbio_HGA.001']=90.2
for s in snl:
	bq[s]={}
for filename in os.listdir('Table_'+ero):
	fb=open('Table_'+ero+'/'+filename,'r')
	line=fb.readline()
	pre=re.split('\.',filename)[0]
	for s in snl:
		if pre not in bq[s]:
			if re.search(s,pre):
				bq[s][pre]={}
				bq[s][pre]['completeness']={}
				bq[s][pre]['contamination']={}
				bq[s][pre]['quality']={}
				bq[s][pre]['species']={}
				bq[s][pre]['species_per']={}
				bq[s][pre]['genus']={}
				bq[s][pre]['genus_per']={}
				bq[s][pre]['N50']={}
				bq[s][pre]['coverage']={}
				bq[s][pre]['abundance']={}
	while True:
		s=''
		line=fb.readline().strip()
		if not line:break
		bin_name=line.split('\t')[0]
		quality=line.split('\t')[-1]
		completeness=line.split('\t')[1]
		contamination=line.split('\t')[2]
		species=line.split('\t')[7]
		N50=line.split('\t')[-5]
		coverage=line.split('\t')[-3]
		species=re.sub('\(.*','',species)
		#sp_per
		species_per=line.split('\t')[7]
		species_per=re.sub('.*\(','',species_per)
		species_per=re.sub('\)','',species_per)
		abundance=line.split('\t')[-2]
		genus=line.split('\t')[9]
		genus=re.sub('\(.*','',genus)
		#genu_per
		genus_per=line.split('\t')[9]
                genus_per=re.sub('.*\(','',genus_per)
                genus_per=re.sub('\)','',genus_per)
		for sn in snl:
			if re.search(sn,line):s=sn
		#bq[s][pre][bin_name]={}
		'''
		if cov[pre]=='Null':
			bq[s][pre]['coverage']='Null'
		else:
			if 'coverage' not in bq[s][pre]:
				bq[s][pre]['coverage']={}
				bq[s][pre]['coverage'][bin_name]=cov[pre][bin_name]
			else:
				bq[s][pre]['coverage'][bin_name]=cov[pre][bin_name]
		'''
		bq[s][pre]['completeness'][bin_name]=float(completeness)
		bq[s][pre]['contamination'][bin_name]=float(contamination)
		bq[s][pre]['N50'][bin_name]=(float(N50)/1000000)
		bq[s][pre]['quality'][bin_name]=quality
		if float(species_per)>=60:
			bq[s][pre]['species'][bin_name]=species
		bq[s][pre]['coverage'][bin_name]=coverage
		bq[s][pre]['abundance'][bin_name]=abundance
		if float(genus_per)>=60:
			bq[s][pre]['genus'][bin_name]=genus
		if quality not in bq[s][pre]:
			bq[s][pre][quality]=1
		else:
			bq[s][pre][quality]+=1
##### Upset csv Generate ########
#### Both species and genus ####
### Species
sp={}
for key1 in bq:
	if key1 not in sp:
		sp[key1]={}
	for key2 in bq[key1]:
		if key2 not in sp[key1]:
			sp[key1][key2]=[]
		for key3 in bq[key1][key2]['species']:
			#if bq[key1][key2]['quality'][key3]=='Low' or bq[key1][key2]['quality'][key3]=='Other':continue
			if bq[key1][key2]['species'][key3] not in sp[key1][key2]:
				sp[key1][key2].append(bq[key1][key2]['species'][key3])
			else:continue
osp=open('Log/sp.hash','w+')
osp.write(str(sp))
osp.close()
### Genus
ge={}
for key1 in bq:
        if key1 not in ge:
                ge[key1]={}
        for key2 in bq[key1]:
                if key2 not in ge[key1]:
                        ge[key1][key2]=[]
                for key3 in bq[key1][key2]['genus']:
                        if bq[key1][key2]['genus'][key3] not in ge[key1][key2]:
                                ge[key1][key2].append(bq[key1][key2]['genus'][key3])
                        else:continue
oge=open('Log/ge.hash','w+')
oge.write(str(ge))
oge.close()

##### Bin Figure ########
#### 1 . Bar Figure and Level Bar Figure ####
def bar(bq,outdir):
	#bar=Bar("Bin质量分布柱形图","可点击分别查看other/low/medium/high Bin数量")
	bar=Bar("Bin Quality Distribution","other/low/medium/high Bin Num can be filtered")
	tp=[]	
	num_other=[]
	num_low=[]
	num_medium=[]
	num_high=[]
	for key1 in bq:
		if 'Other' not in bq[key1]:
			bq[key1]['Other']=0
		if 'Low' not in bq[key1]:
	                bq[key1]['Low']=0
		if 'Medium' not in bq[key1]:
	                bq[key1]['Medium']=0	
		if 'High' not in bq[key1]:
	        	bq[key1]['High']=0
		tp.append(key1)		
		num_other.append(bq[key1]['Other'])
		num_low.append(bq[key1]['Low'])
		num_medium.append(bq[key1]['Medium'])
		num_high.append(bq[key1]['High'])
	bar.add('Other',tp,num_other)
	bar.add('Low',tp,num_low)
	bar.add('Medium',tp,num_medium)
	bar.add('High',tp,num_high,is_label_show=True,xaxis_rotate=20)
	bar.width=1200
	bar.height=600
	bar.render(outdir+'/1_bar.html')
	#bar.render('Pyechart/bin_bar.html')
	
	#Literature Fig/Barplot
	plt.figure(figsize=(12,6))
	n_groups=len(bq.keys())
	index=np.arange(n_groups)
	barwidth=0.2
	qua_arr=['Low','Medium','High']
	#class_name=sorted(bq.keys())
	class_name=[]
	for key1 in bq:
                class_name.append(key1)
	count_qua=-1
	color_qua={}
	#color_bar=['#d62728','#9467bd','#1f77b4','#ff7f0e']
	#color_bar=['blue','cyan','green','olive']
	#color_bar=['red','pink','purple','orange']
	color_bar=['red','#04477C','#1291A9','#DA891E']
	cn_bar=0
	for q in qua_arr:
		color_qua[q]=color_bar[cn_bar]
		cn_bar+=1
	for q in qua_arr:
		y=[]
		count_qua+=1
		for c in class_name:
			y.append(bq[c][q])	
		plt.bar(index+(count_qua*barwidth),y,width=0.2,facecolor=color_qua[q],edgecolor='white',label=q)
	plt.xticks(index+barwidth/count_qua,class_name,rotation=11,fontsize=18)
	plt.ylabel('number  of  bins',fontsize=18)
	#plt.legend(loc='upper right')
	plt.legend(loc='center left', bbox_to_anchor=(0.2, 1.08),ncol=4,fontsize=18)
	plt.savefig(outdir+'/1.1_bar.png')
			
				
	
	###Stack Bar plot ######
	#bar2=Bar("Bin数据堆叠柱状图","显示不同质量Bin占比组成\n\n")
	bar=Bar("Bin Quality Distribution (stack)","Show the proportion of different quality bins\n\n")
	attr=[]
	other=[]
	low=[]
	medium=[]
	high=[]
	for key1 in bq:
		attr.append(key1)
		if 'Other' not in bq[key1]:
                        bq[key1]['Other']=0
                if 'Low' not in bq[key1]:
                        bq[key1]['Low']=0
                if 'Medium' not in bq[key1]:
                        bq[key1]['Medium']=0
                if 'High' not in bq[key1]:
                        bq[key1]['High']=0			
		total=bq[key1]['Other']+bq[key1]['Low']+bq[key1]['Medium']+bq[key1]['High']
		other.append(float(float(bq[key1]['Other'])/float(total)))	
		low.append(float(float(bq[key1]['Low'])/float(total)))
		medium.append(float(float(bq[key1]['Medium'])/float(total)))
		high.append(float(float(bq[key1]['High'])/float(total)))
	bar2.add('Other',attr,other,is_stack=True)	
	bar2.add('Low',attr,low,is_stack=True)
	bar2.add('Medium',attr,medium,is_stack=True)
	bar2.add('High',attr,high,is_stack=True,xaxis_rotate=20)
	bar2.width=1200
	bar2.height=600
	bar2.render(outdir+'/2_lb.html')
	#Literature fig /Stack bar plot
	plt.figure(figsize=(12,6))
        n_groups=len(bq.keys())
        index=np.arange(n_groups)
	width=0.35
	'''
	mul_qua_arr=[]
	#mul_qua_arr.append(other)
	mul_qua_arr.append(low)
	mul_qua_arr.append(medium)
	mul_qua_arr.append(high)
	count_bar_stack=0
	for m in mul_qua_arr:
		if count_bar_stack==0:
			plt.bar(index,m,width,color=color_qua[qua_arr[count_bar_stack]],label=qua_arr[count_bar_stack])
			tem=m
			
		else:
			if count_bar_stack==1:
				plt.bar(index,m,width,color=color_qua[qua_arr[count_bar_stack]],label=qua_arr[count_bar_stack],bottom=tem)
				s=np.array(m)+np.array(tem)
			else:
				plt.bar(index,m,width,color=color_qua[qua_arr[count_bar_stack]],label=qua_arr[count_bar_stack],bottom=s)
				s+=np.array(m)
				
		count_bar_stack+=1
	'''
	plt.bar(index,low,label='Low',fc=color_qua['Low'])
	plt.bar(index,medium,bottom=low,label='Medium',fc=color_qua['Medium'])
	base_medium=np.array(low)+np.array(medium)
	plt.bar(index,high,bottom=base_medium,label='High',fc=color_qua['High'])
	plt.xticks(index,attr,rotation=11,fontsize=18)
	plt.ylabel(' proportion  of  bins',fontsize=18)
	#plt.legend(loc='upper right')
	plt.legend(loc='center left', bbox_to_anchor=(0.2, 1.08),ncol=4,fontsize=18)
	plt.savefig(outdir+'/1.2_stack_bar.png')
	
			

##### 2.Black(completeness) and Red(contamination) Scatter Plot ######
sc={}
#name=[]
def scatter(bq,output):
	page=Page()
	page.height=900
	c=0
	fc=18   #fontsize
	if len(bq.keys())>4:
		fig=plt.figure(figsize=(24,10))
	else:
		fig=plt.figure(figsize=(20,5))
	#fig.subplots_adjust(wspace=0.5,hspace=0.5)
	#iden the fig shape
	#2 col situation
	if len(bq.keys())%2==0 and len(bq.keys())%3!=0 and len(bq.keys())%4!=0:
		all_col=2
	elif len(bq.keys())%3==0:
		all_col=3
	elif len(bq.keys())%4==0:
		all_col=4
	#iden shape end
	if int(len(bq.keys())%all_col)==0:
		row=int(len(bq.keys())/all_col)
	else:
		row=int(len(bq.keys())/all_col)+1
	for key1 in bq:
		cp=[] #completeness
		ct=[] #contamination
		a=sorted(bq[key1]['completeness'].items(),key=lambda x:x[1],reverse=True) #sort keys  according to their values
		#c=0
		for e in a:
			cp.append(bq[key1]['completeness'][e[0]])
			ct.append(bq[key1]['contamination'][e[0]])
		c+=1
		x=range(1,len(cp)+1)
		##### literature  figure part #####
		#fig.subplots_adjust(wspace=0.3,hspace=0.6)
		if max(ct)>100:
			ax1=fig.add_subplot(row,all_col,c)
                        ax1.set_xlabel(key1,color='black',fontsize=fc)
                        ax1.set_ylabel('Completeness',fontsize=fc)
                        ax1.scatter(x,cp,color='black',label='Completeness')
                        ax2=ax1.twinx()
                        ax2.scatter(x,ct,color='red',label='Contamination')
                        ax2.set_ylabel('Contamination',color='black',fontsize=fc)
                        ax2.tick_params(axis='y',colors='red')
                        ax2.spines['right'].set_color('red')
			if max(ct)%10==0:
				max_10_ct=max(ct)
			else:
				max_10_ct=(10-max(ct)%10)+max(ct)
			ax2.axis([min(x),max(x),0,max_10_ct])
			ax1.axis([min(x),max(x),0,max_10_ct])
			#ax2.axis([min(x),max(x),0,max(cp)])
		else:
			ax1=fig.add_subplot(row,all_col,c)
			ax1.set_xlabel(key1,color='black',fontsize=fc)
			ax1.set_ylabel('Completeness',fontsize=fc)
			ax1.scatter(x,cp,color='black',label='Completeness')
			ax2=ax1.twinx()
			ax2.scatter(x,ct,color='red',label='Contamination')
			ax2.set_ylabel('Contamination',color='black',fontsize=fc)
			ax2.tick_params(axis='y',colors='red')
			ax2.spines['right'].set_color('red')
			ax2.axis([min(x),max(x),0,100])
		#plt.savefig(output+'/3.scatter.png')				
		##### literature  figure end  #####

		#x=range(1,len(cp)+1)
		#grid=Grid()
		#overlap = Overlap()
		#scatter = Scatter()
		#scatter2 = Scatter(key1)
		if c%2==0:
			#2
			if c%4!=0:
				overlap2 = Overlap()
		                scatter3 = Scatter()
				scatter3.add("Completeness",x,cp,yaxis_label_textcolor='red',yaxis_line_color='red',xaxis_type='category',xaxis_max=max(x),legend_pos="40%",legend_orient="vertical")
				scatter4 = Scatter(key1,title_pos="25%")
				scatter4.add('Contamination',x,ct,xaxis_type='category',xaxis_max=max(x),legend_pos="40%",legend_orient="vertical")
				overlap2.add(scatter4)
		                overlap2.add(scatter3)
			#4
			else:
				overlap4 = Overlap()
                                scatter7 = Scatter()
                                scatter7.add("Completeness",x,cp,yaxis_label_textcolor='red',yaxis_line_color='red',xaxis_type='category',xaxis_max=max(x),legend_pos="90%",legend_orient="vertical")
                                scatter8 = Scatter(key1,title_pos="75%")
                                scatter8.add('Contamination',x,ct,xaxis_type='category',xaxis_max=max(x),legend_pos="90%",legend_orient="vertical")
                                overlap4.add(scatter8)
                                overlap4.add(scatter7)
			#overlap2.render(output+'/3.scatter_o2.html')
		else:
			#1
			if (c+1)%4!=0:
				overlap = Overlap()
		                scatter = Scatter()
				scatter.add("Completeness",x,cp,yaxis_label_textcolor='red',yaxis_line_color='red',xaxis_type='category',xaxis_max=max(x),legend_pos="15%",legend_orient="vertical")
	                        scatter2 = Scatter(key1)
                        	scatter2.add('Contamination',x,ct,xaxis_type='category',xaxis_max=max(x),legend_pos="15%",legend_orient="vertical")
				overlap.add(scatter2)
				overlap.add(scatter)
			#3
			else:
				overlap3 = Overlap()
                                scatter5= Scatter()
                                scatter5.add("Completeness",x,cp,yaxis_label_textcolor='red',yaxis_line_color='red',xaxis_type='category',xaxis_max=max(x),legend_pos="65%",legend_orient="vertical")
                                scatter6= Scatter(key1,title_pos="50%")
                                scatter6.add('Contamination',x,ct,xaxis_type='category',xaxis_max=max(x),legend_pos="65%",legend_orient="vertical")
                                overlap3.add(scatter6)
                                overlap3.add(scatter5)
			#overlap.render(output+'/3.scatter_o1.html')
		'''
		if c==len(bq.keys()) and c%2!=0:
			grid=Grid()
                        grid.add(overlap,grid_left="60%")
			page.add(grid)
			continue
		'''
		if c%2==0:
			#2
			if c%4!=0:
				grid.add(overlap2,grid_left="30%",grid_width=250)
			#4
			else:
				grid.add(overlap4,grid_left="80%",grid_width=200)
                                page.add(grid)
			#grid.render(output+'/3.scatter.html')
			#page.render(output+'/3.scatter.html')
		else:
			#3
			if (c+1)%4==0:
				grid.add(overlap3, grid_left="55%",grid_width=250)
			#1
			else:
				grid=Grid()
				grid.width=1500
				grid.height=300
				grid.add(overlap,grid_left="5%",grid_width=250)
			#grid.render(output+'/3.scatter.html')
                        #page.render(output+'/3.scatter.html')
                        #exit()
		if c==len(bq.keys()) and c%4!=0:
			page.add(grid)
	page.render(output+'/3_scatter.html')
	#fig.subplots_adjust(wspace=0.3,hspace=0.3)
	#ax1.legend(loc='center', bbox_to_anchor=(-1.3, 3.8),ncol=5,fontsize=20)
	fig.subplots_adjust(wspace=0.5)
	plt.savefig(output+'/3.scatter.png')
###### 3.Upset Plot ########################################
def upset(sam,output): #sam -> sample_name
	pwd=os.getcwd()
	fu=open('Log/sp.hash','r')
	fg=open('Log/ge.hash','r')
	fr=open('Log/upset.R','w+')
	a=fu.read()   #species
	b=fg.read()   #genus
	d=eval(a)     #species
	f=eval(b)     #genus
	tem={}
	all_species={}
	all_genus={}
	m=0
	#genus
	for key1 in f:
		af=[]
		tem=f[key1]
		for key2 in f[key1]:
			for e in f[key1][key2]:
				if e not in af:
					af.append(e)
		all_genus[key1]=af
	for key in f:
		inhash={}
		inhash['name']=[]
		ap=all_genus[key]
		for a in ap:
			inhash['name'].append(a)
		for key2 in f[key]:
			inhash[key2]=[]
			for a in ap:
				if a in f[key1][key2]:
					inhash[key2].append(1)
				else:
					inhash[key2].append(0)
		frame=pd.DataFrame(inhash)
	        frame.to_csv('Log/'+key+'_genus.csv', index=False, header=True)
	fr.write('library(UpSetR)\npng("'+pwd+'/'+output+'/4_upset_genus.png",width=1200,height=800,res=72*2)\n')
        fr.write('frame<-read.csv("'+pwd+'/Log/'+sam+'_genus.csv",header = TRUE, sep=",")\n')
        fr.write('upset(frame, nsets = 12, nintersects = 30,order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))\n')
	fr.write('dev.off()\n')
	#species
	for key1 in d:
		at=[]
		tem=d[key1]
		for key2 in d[key1]:
                	for e in d[key1][key2]:
	                        if e not in at:
        	                        at.append(e)
	        all_species[key1]=at
	for key in d:
		inhash={}
		inhash['name']=[]
		al=all_species[key]
		for a in al:
			inhash['name'].append(a)
		for key2 in d[key]:
			inhash[key2]=[]
			for a in al:
				if a in d[key][key2]:
					inhash[key2].append(1)
				else:
					inhash[key2].append(0)
		frame=pd.DataFrame(inhash)
	        frame.to_csv('Log/'+key+'_species.csv', index=False, header=True)
	fr.write('library(UpSetR)\npng("'+pwd+'/'+output+'/4_upset_species.png",width=1200,height=800,res=72*2)\n')
	fr.write('frame<-read.csv("'+pwd+'/Log/'+sam+'_species.csv",header = TRUE, sep=",")\n')
	fr.write('upset(frame, nsets = 12, nintersects = 30,order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))\n')
	fr.write('dev.off()')
	fr.close()
	robjects.r.source('Log/upset.R')	
###### 4. Bin N50 and Coverage ########
def bin_boxplot(bq,output):
	cpd={} #eg {'Low':'hlj_pacbio_HGA':'coverage':[12.3,12,3,1.5,102.5]}
	cpd['High']={}
        cpd['Medium']={}
        cpd['Low']={}
	overall={}  #eg overall['hlj_pacbio_HGA']['coverage']=[12.3,12,3,1.5,102.5,...]	
	plot_dict={}
	for key1 in bq:
		cpd['High'][key1]={}
		cpd['Medium'][key1]={}
		cpd['Low'][key1]={}
		if not bq[key1]['coverage']=='Null':
                        cpd['High'][key1]['coverage']=[]
                        cpd['Medium'][key1]['coverage']=[]
                        cpd['Low'][key1]['coverage']=[]
		cpd['High'][key1]['N50']=[]
                cpd['Medium'][key1]['N50']=[]
                cpd['Low'][key1]['N50']=[]
		#overall part
                overall[key1]={}
                overall[key1]['coverage']=[]
                overall[key1]['N50']=[]
		for key2 in bq[key1]['quality']:
                        #overall part
                        overall[key1]['coverage'].append(float(bq[key1]['coverage'][key2]))
                        overall[key1]['N50'].append(float(bq[key1]['N50'][key2]))
                        #split part
                        if bq[key1]['quality'][key2] in cpd:
                                if not bq[key1]['coverage']=='Null':
                                        cpd[bq[key1]['quality'][key2]][key1]['coverage'].append(float(bq[key1]['coverage'][key2]))
                                cpd[bq[key1]['quality'][key2]][key1]['N50'].append(float(bq[key1]['N50'][key2]))
	labels=sorted(overall.keys())
	plot_dict['overall']={}
	for key in cpd:
		plot_dict[key]={}
	for key in plot_dict:
		plot_dict[key]['coverage']=[]
		plot_dict[key]['N50']=[]
	for key in labels:
		plot_dict['overall']['coverage'].append(overall[key]['coverage'])
		plot_dict['overall']['N50'].append(overall[key]['N50'])
	for key in cpd:
		for l in labels:
			plot_dict[key]['coverage'].append(cpd[key][l]['coverage'])
			plot_dict[key]['N50'].append(cpd[key][l]['N50'])
	qtype=['overall','High','Medium','Low']
	feature=['N50','coverage']
	for q in qtype:
		fig=plt.figure(figsize=(12,6))
		#fig.subplots_adjust(wspace=0.3)
		c=1
		for f in feature:
			ax=fig.add_subplot(1,2,c)
			if f=='N50':
				ax.set_ylabel(f+'(Mb)',fontsize=18)
			else:
				ax.set_ylabel(f,fontsize=18)
			ax.boxplot(plot_dict[q][f],labels=labels)
			plt.xticks(rotation=15,fontsize=15)
			ax.set_title(q,fontsize=18)
			c+=1
		fig.subplots_adjust(wspace=0.2)
		fig.savefig(output+'/5_boxplot_cov_n50_'+q+'.png')

	
		
		
		
		
def bin_per(bq,output):	
	cpd={} #eg {'Low':'hlj_pacbio_HGA':'coverage':[12.3,12,3,1.5,102.5]}
	cpd['High']={}
	cpd['Medium']={}
	cpd['Low']={}
	overall={}  #eg overall['hlj_pacbio_HGA']['coverage']=[12.3,12,3,1.5,102.5,...]
	for key1 in bq:
		cpd['High'][key1]={}
		cpd['Medium'][key1]={}
		cpd['Low'][key1]={}
		if not bq[key1]['coverage']=='Null':
			cpd['High'][key1]['coverage']=[]
			cpd['Medium'][key1]['coverage']=[]
			cpd['Low'][key1]['coverage']=[]
		cpd['High'][key1]['N50']=[]
                cpd['Medium'][key1]['N50']=[]
                cpd['Low'][key1]['N50']=[]
		#overall part
		overall[key1]={}
		overall[key1]['coverage']=[]
		overall[key1]['N50']=[]
		for key2 in bq[key1]['quality']:
			#overall part
			overall[key1]['coverage'].append(bq[key1]['coverage'][key2])
			overall[key1]['N50'].append(bq[key1]['N50'][key2])
			#split part
			if bq[key1]['quality'][key2] in cpd:
				if not bq[key1]['coverage']=='Null':
					cpd[bq[key1]['quality'][key2]][key1]['coverage'].append(bq[key1]['coverage'][key2])
				cpd[bq[key1]['quality'][key2]][key1]['N50'].append(bq[key1]['N50'][key2])
	#add color dict
	color_list=['red','blue','green','saddlebrown','black','deepskyblue','blueviolet','magenta','cyan','#0cff0c']   #10 color now
	colordict={}
	for key1 in cpd:
		cnum=0
		for key2 in cpd[key1]:	
			if key2 not in colordict:
				colordict[key2]=color_list[cnum]
				cnum+=1
	fig=plt.figure(figsize=(20,5))	
	ax=plt.figure().add_subplot(111)
	if True:
		with plt.style.context('ggplot'):
			fig.subplots_adjust(wspace=0.3)	
			ax1=fig.add_subplot(1,4,1)
                        ax2=fig.add_subplot(1,4,2)
                        ax3=fig.add_subplot(1,4,3)
                        ax4=fig.add_subplot(1,4,4)
			ax1.set_xlabel('N50(Mb)',fontsize=15,color='black')
                        ax1.set_ylabel('number  of bins',fontsize=15,color='black')
                        #ax1.text(-0.2,18,"A",fontsize=22)
                        ax2.set_xlabel('N50(Mb)',fontsize=15,color='black')
                        ax2.set_ylabel('proportion of bins',fontsize=15,color='black')
                        ax3.set_xlabel('Coverage',fontsize=15,color='black')
                        ax3.set_ylabel('number of bins',fontsize=15,color='black')
			ax4.set_xlabel('Coverage',fontsize=15,color='black')
                        ax4.set_ylabel('proportion of bins',fontsize=15,color='black')
			for key in overall:
				le=key
				n,bins1,patches=ax.hist(overall[key]['N50'],bins=100,cumulative=-1,label=None,color=None,histtype='step')
                                n1,bins,patches1=ax.hist(overall[key]['N50'],bins=99,cumulative=-1,label=None,color=None,histtype='step')
                                ax1.plot(bins,n,linewidth=2,label=le,color=colordict[key])
				n,bins1,patches=ax.hist(overall[key]['N50'],bins=100,cumulative=1,normed=1,label=None,color=None,histtype='step')
                                n1,bins,patches1=ax.hist(overall[key]['N50'],bins=99,cumulative=1,normed=1,label=None,color=None,histtype='step')
                                ax2.plot(bins,n,linewidth=2,label=le,color=colordict[key])
				n,bins1,patches=ax.hist(overall[key]['coverage'],bins=100,cumulative=True,label=None,color=None,histtype='step')
                                n1,bins,patches1=ax.hist(overall[key]['coverage'],bins=99,cumulative=True,label=None,color=None,histtype='step')
                                ax3.plot(bins,n,linewidth=2,label=le,color=colordict[key])
				n,bins1,patches=ax.hist(overall[key]['coverage'],bins=100,cumulative=True,normed=1,label=None,color=None,histtype='step')
                                n1,bins,patches1=ax.hist(overall[key]['coverage'],bins=99,cumulative=True,normed=1,label=None,color=None,histtype='step')
                                ax4.plot(bins,n,linewidth=2,label=le,color=colordict[key])
	#ax4.legend(loc='center',bbox_to_anchor=(1.4,0.9),ncol=1,fontsize=10)
	#ax4.legend(loc='center ', bbox_to_anchor=(-0.2, 1.12),ncol=3,fontsize=10)
	#ax4.legend(loc='upper center')
	#ax4.legend(loc='upper left',bbox_to_anchor=(1.0,0.7),fontsize=9)
	ax4.legend(loc='center', bbox_to_anchor=(-1.4,1.1 ),ncol=5,fontsize=11)
        plt.tight_layout()
        fig.savefig(output+'/5_bin_cov_n50_overall.png')				
	#plot parameters
        #fig=plt.figure(figsize=(35,14),dpi=1000)
	#fig=plt.figure(figsize=(40,20))
	#fig=plt.figure(figsize=(35,18),constrained_layout=True)
	fig=plt.figure(figsize=(36,14))
        ax=plt.figure().add_subplot(111)
        count=1
        linew=2
        fs=20
	arr=['High','Medium','Low']
	for a in arr:
		with plt.style.context('ggplot'):
			fig.subplots_adjust(wspace=0.2,hspace=0.2)		
			#color map related
			'''
			legend_num=len(cpd[a])
			values=range(legend_num)
			jet=cm=plt.get_cmap('jet')
			cNorm=colors.Normalize(vmin=0, vmax=values[-1])
			scalarMap=cmx.ScalarMappable(norm=cNorm, cmap=jet)
			'''
			#color map down
			if count==1:
				max_1_x=0
				max_1_y=0
				max_2_x=0
				max_3_x=0
				max_3_y=0
				max_4_x=0
				ax1=fig.add_subplot(3,4,count)	
				ax2=fig.add_subplot(3,4,count+1)
	                        ax3=fig.add_subplot(3,4,count+2)
				ax4=fig.add_subplot(3,4,count+3)
				ax1.set_xlabel('N50(Mb)',fontsize=fs,color='black')
				ax1.set_ylabel('number  of bins',fontsize=fs,color='black')
				ax1.text(-0.2,18,"A",fontsize=fs+5)
				plt.xticks(fontsize=10)
				ax2.set_xlabel('N50(Mb)',fontsize=fs,color='black')
				ax2.set_ylabel('proportion of bins',fontsize=fs,color='black')
				plt.xticks(fontsize=10)
				ax1.text(2.3,18,"B",fontsize=fs+5)
				ax3.set_xlabel('Coverage',fontsize=fs,color='black')
        	                ax3.set_ylabel('number of bins',fontsize=fs,color='black')
				plt.xticks(fontsize=10)
				ax1.text(4.8,18,"C",fontsize=fs+5)
				ax4.set_xlabel('Coverage',fontsize=fs,color='black')
                        	ax4.set_ylabel('proportion of bins',fontsize=fs,color='black')
				plt.xticks(fontsize=10)
				ax1.text(7.3,18,"D",fontsize=fs+5)
				for key2 in cpd[a]:
					le=key2   #le means legend
					
					#N50
					n,bins1,patches=ax.hist(cpd[a][key2]['N50'],bins=100,cumulative=-1,label=None,color=None,histtype='step')
					n1,bins,patches1=ax.hist(cpd[a][key2]['N50'],bins=99,cumulative=-1,label=None,color=None,histtype='step')
					ax1.plot(bins,n,linewidth=linew,label=le,color=colordict[key2])
					if max(bins)>max_1_x:
						max_1_x=max(bins)
					if max(n)>max_1_y:
						max_1_y=max(n)
					ax1.axis([0,max_1_x,0,max_1_y])
					plt.xticks(fontsize=18)
					n,bins1,patches=ax.hist(cpd[a][key2]['N50'],bins=100,cumulative=1,normed=1,label=None,color=None,histtype='step')
                                        n1,bins,patches1=ax.hist(cpd[a][key2]['N50'],bins=99,cumulative=1,label=None,normed=1,color=None,histtype='step')
					ax2.plot(bins,n,linewidth=linew,label=le,color=colordict[key2])
					if max(bins)>max_2_x:
						max_2_x=max(bins)
                                        ax2.axis([0,max_2_x,0,1])
                                        plt.xticks(fontsize=18)	
					#coverage
					n,bins3,patches=ax.hist(cpd[a][key2]['coverage'],bins=100,cumulative=True,label=None,color=None,histtype='step')
					n3,bins,patches3=ax.hist(cpd[a][key2]['coverage'],bins=99,cumulative=True,label=None,color=None,histtype='step')
                                        ax3.plot(bins,n,linewidth=linew,label=le,color=colordict[key2])
					if max(bins)>max_3_x:
						max_3_x=max(bins)
					if max(n)>max_3_y:
						max_3_y=max(n)
                                        ax3.axis([0,max_3_x,0,max_3_y])
                                        plt.xticks(fontsize=18)
                                        n,bins,patches=ax.hist(cpd[a][key2]['coverage'],bins=100,normed=1,cumulative=True,label=None,color=None,histtype='step')
					n3,bins,patches3=ax.hist(cpd[a][key2]['coverage'],bins=99,normed=1,cumulative=True,label=None,color=None,histtype='step')
                                        ax4.plot(bins,n,linewidth=linew,label=le,color=colordict[key2])
					if max(bins)>max_4_x:
						max_4_x=max(bins)
                                        ax4.axis([0,max_4_x,0,1])
                                        plt.xticks(fontsize=18)
				plt.xticks(fontsize=18)
                                ax1.axis(fontsize=12)				
			if count==5:
				max_5_x=0
                                max_5_y=0
                                max_6_x=0
                                max_7_x=0
                                max_7_y=0
                                max_8_x=0
				ax5=fig.add_subplot(3,4,count)
                                ax6=fig.add_subplot(3,4,count+1)
                                ax7=fig.add_subplot(3,4,count+2)
                                ax8=fig.add_subplot(3,4,count+3)
                                ax5.set_xlabel('N50(Mb)',fontsize=fs,color='black')
                                ax5.set_ylabel('number of bins',fontsize=fs,color='black')
                                ax1.text(-0.2,-3,"E",fontsize=fs+5)
                                ax6.set_xlabel('N50(Mb)',fontsize=fs,color='black')
                                ax6.set_ylabel('proportion of bins',fontsize=fs,color='black')
                                ax1.text(2.3,-3,"F",fontsize=fs+5)
                                ax7.set_xlabel('Coverage',fontsize=fs,color='black')
                                ax7.set_ylabel('number of bins',fontsize=fs,color='black')
                                ax1.text(4.8,-3,"G",fontsize=fs+5)
                                ax8.set_xlabel('Coverage',fontsize=fs,color='black')
                                ax8.set_ylabel('proportion of bins',fontsize=fs,color='black')
                                ax1.text(7.3,-3,"H",fontsize=fs+5)
				for key2 in cpd[a]:
                                        le=key2   #le means legend
                                        #N50
                                        n,bins1,patches=ax.hist(cpd[a][key2]['N50'],bins=100,cumulative=-1,label=None,color=None,histtype='step')
                                        n1,bins,patches1=ax.hist(cpd[a][key2]['N50'],bins=99,cumulative=-1,label=None,color=None,histtype='step')
                                        ax5.plot(bins,n,linewidth=linew,label=le,color=colordict[key2])
					if max(bins)>max_5_x:
						max_5_x=max(bins)
					if max(n)>max_5_y:
						max_5_y=max(n)
                                        ax5.axis([0,max_5_x,0,max_5_y])
                                        plt.xticks(fontsize=10)
                                        n,bins1,patches=ax.hist(cpd[a][key2]['N50'],bins=100,cumulative=1,normed=1,label=None,color=None,histtype='step')
                                        n1,bins,patches1=ax.hist(cpd[a][key2]['N50'],bins=99,cumulative=1,label=None,normed=1,color=None,histtype='step')
                                        ax6.plot(bins,n,linewidth=linew,label=le,color=colordict[key2])
					if max(bins)>max_6_x:
						max_6_x=max(bins)
                                        ax6.axis([0,max_6_x,0,1])
                                        plt.xticks(fontsize=10)
                                        #coverage
                                        n,bins3,patches=ax.hist(cpd[a][key2]['coverage'],bins=100,cumulative=True,label=None,color=None,histtype='step')
                                        n3,bins,patches3=ax.hist(cpd[a][key2]['coverage'],bins=99,cumulative=True,label=None,color=None,histtype='step')
                                        ax7.plot(bins,n,linewidth=linew,label=le,color=colordict[key2])
					if max(bins)>max_7_x:
						max_7_x=max(bins)
					if max(n)>max_7_y:
						max_7_y=max(n)
                                        ax7.axis([0,max_7_x,0,max_7_y])
                                        plt.xticks(fontsize=10)
                                        n,bins,patches=ax.hist(cpd[a][key2]['coverage'],bins=100,normed=1,cumulative=True,label=None,color=None,histtype='step')
                                        n3,bins,patches3=ax.hist(cpd[a][key2]['coverage'],bins=99,normed=1,cumulative=True,label=None,color=None,histtype='step')
                                        ax8.plot(bins,n,linewidth=linew,label=le,color=colordict[key2])
					if max(bins)>max_8_x:
						max_8_x=max(bins)
                                        ax8.axis([0,max_8_x,0,1])
                                        plt.xticks(fontsize=10)
			if count==9:
				max_9_x=0
                                max_9_y=0
                                max_10_x=0
                                max_11_x=0
                                max_11_y=0
                                max_12_x=0
				ax9=fig.add_subplot(3,4,count)
                                ax10=fig.add_subplot(3,4,count+1)
                                ax11=fig.add_subplot(3,4,count+2)
                                ax12=fig.add_subplot(3,4,count+3)
                                ax9.set_xlabel('N50(Mb)',fontsize=fs,color='black')
                                ax9.set_ylabel('number of bins',fontsize=fs,color='black')
                                ax1.text(-0.2,-24,"I",fontsize=fs+5)
                                ax10.set_xlabel('N50(Mb)',fontsize=fs,color='black')
                                ax10.set_ylabel('proportion of bins',fontsize=fs,color='black')
                                ax1.text(2.3,-24,"J",fontsize=fs+5)
                                ax11.set_xlabel('Coverage',fontsize=fs,color='black')
                                ax11.set_ylabel('number of bins',fontsize=fs,color='black')
                                ax1.text(4.8,-24,"K",fontsize=fs+5)
                                ax12.set_xlabel('Coverage',fontsize=fs,color='black')
                                ax12.set_ylabel('proportion of bins',fontsize=fs,color='black')
                                ax1.text(7.3,-24,"L",fontsize=fs+5)
				for key2 in cpd[a]:
                                        le=key2   #le means legend
                                        #N50
                                        n,bins1,patches=ax.hist(cpd[a][key2]['N50'],bins=100,cumulative=-1,label=None,color=None,histtype='step')
                                        n1,bins,patches1=ax.hist(cpd[a][key2]['N50'],bins=99,cumulative=-1,label=None,color=None,histtype='step')
                                        ax9.plot(bins,n,linewidth=linew,label=le,color=colordict[key2])
					if max(bins)>max_9_x:
                                                max_9_x=max(bins)
                                        if max(n)>max_9_y:
                                                max_9_y=max(n)
                                        ax9.axis([0,max_9_x,0,max_9_y])
                                        plt.xticks(fontsize=10)
                                        n,bins1,patches=ax.hist(cpd[a][key2]['N50'],bins=100,cumulative=1,normed=1,label=None,color=None,histtype='step')
                                        n1,bins,patches1=ax.hist(cpd[a][key2]['N50'],bins=99,cumulative=1,label=None,normed=1,color=None,histtype='step')
                                        ax10.plot(bins,n,linewidth=linew,label=le,color=colordict[key2])
					if max(bins)>max_10_x:
						max_10_x=max(bins)
                                        ax10.axis([0,max_10_x,0,1])
                                        plt.xticks(fontsize=10)
                                        #coverage
                                        n,bins3,patches=ax.hist(cpd[a][key2]['coverage'],bins=100,cumulative=True,label=None,color=None,histtype='step')
                                        n3,bins,patches3=ax.hist(cpd[a][key2]['coverage'],bins=99,cumulative=True,label=None,color=None,histtype='step')
                                        ax11.plot(bins,n,linewidth=linew,label=le,color=colordict[key2])
					if max(bins)>max_11_x:
                                                max_11_x=max(bins)
                                        if max(n)>max_11_y:
                                                max_11_y=max(n)
                                        ax11.axis([0,max_11_x,0,max_11_y])
                                        #plt.xticks(fontsize=10)
                                        n,bins,patches=ax.hist(cpd[a][key2]['coverage'],bins=100,normed=1,cumulative=True,label=None,color=None,histtype='step')
                                        n3,bins,patches3=ax.hist(cpd[a][key2]['coverage'],bins=99,normed=1,cumulative=True,label=None,color=None,histtype='step')
                                        ax12.plot(bins,n,linewidth=linew,label=le,color=colordict[key2])
					if max(bins)>max_12_x:
                                                max_12_x=max(bins)
                                        ax12.axis([0,max_12_x,0,1])
                                        #plt.xticks(fontsize=10)
			count+=4
										
	#ax12.legend(loc='center',bbox_to_anchor=(1.3,2.9),ncol=1,fontsize=18)
	ax12.legend(loc='center', bbox_to_anchor=(-1.3, 3.8),ncol=5,fontsize=20)
	plt.tight_layout()
	fig.savefig(output+'/5_bin_cov_n50.png')	
	
###### 5. Species/Genus Relative Abundance  ########		
def tax_abundance(bq,output):
	genus={}
	species={}
	for key1 in bq:		
		genus[key1]={}
		species[key1]={}
		for key2 in bq[key1]['genus']:
			if bq[key1]['genus'][key2] not in genus[key1]:
				genus[key1][bq[key1]['genus'][key2]]=float(bq[key1]['abundance'][key2])
			else:
				genus[key1][bq[key1]['genus'][key2]]+=float(bq[key1]['abundance'][key2])
		for key2 in bq[key1]['species']:
			if bq[key1]['species'][key2] not in species[key1]:
				species[key1][bq[key1]['species'][key2]]=float(bq[key1]['abundance'][key2])
			else:
				species[key1][bq[key1]['species'][key2]]+=float(bq[key1]['abundance'][key2])
	uni_spe={}
	uni_genus={}
	for key1 in species:
		for key2 in species[key1]:
			if key2 not in uni_spe:
				uni_spe[key2]=1
			else:
				uni_spe[key2]+=1
	for key1 in genus:
		for key2 in genus[key1]:
			if key2 not in uni_genus:
				uni_genus[key2]=1
			else:
				uni_genus[key2]+=1
	#color related
        color_bar=['#fc5c65','#fd9644','#fed330','#26de81','#2bcbba','#45aaf2','#4b7bec','#a55eea','#d1d8e0','#778ca3','#2C3A47','#B33771','#485460','#0be881','#3742fa']
	#Genus stack bar
	plt.figure(figsize=(12,11))
	class_name=sorted(genus.keys())
        n_groups=len(genus.keys())
        index=np.arange(n_groups)
        all_genus=sorted(uni_genus.keys())
        gen_arr=sorted(genus.keys())
        width=0.35
	mul_gen_arr=[]
	total_genus={}	
	#color start
	color_genus={}
        cn_bar=0
        for s in all_genus:
                color_genus[s]=color_bar[cn_bar]
                cn_bar+=1
	#color end
	for key1 in sorted(genus.keys()):
		total=0
		for g in all_genus:
			if g not in genus[key1]:continue
			else:
				total+=genus[key1][g]
		total_genus[key1]=total
	for s in all_genus:
                tem=[]
                for key1 in sorted(genus.keys()):
                        if s not in genus[key1]:
                                tem.append(0)
                        else:
                                tem.append(float(genus[key1][s])/float(total_genus[key1]))
                mul_gen_arr.append(tem)
	count_bar_stack=0
	#Pyechart Report Part
	#bar_genus=Bar("物种丰度数据堆叠柱状图","Genus层面\n\n")
	bar_genus=Bar("Relative abundance","Genus Level\n\n")
	attr=sorted(bq.keys())
	c=0
	for k in all_genus:
		bar_genus.add(k,attr,mul_gen_arr[c],is_stack=True,xaxis_rotate=20,legend_top='10%',legend_pos='10%',yaxis_max=1)
		c+=1
	gg=Grid(width=1200,height=600)
        gg.add(bar_genus,grid_top='20%')
	gg.render(output+'/6_stack_genus.html')
	#literature figure part
	for m in mul_gen_arr:
                if count_bar_stack==0:
                        plt.bar(index,m,width,color=color_genus[all_genus[count_bar_stack]],label=all_genus[count_bar_stack])
                        tem=m
                else:
                        if count_bar_stack==1:
                                plt.bar(index,m,width,color=color_genus[all_genus[count_bar_stack]],label=all_genus[count_bar_stack],bottom=tem)
                                s=np.array(m)+np.array(tem)
                        else:
                                plt.bar(index,m,width,color=color_genus[all_genus[count_bar_stack]],label=all_genus[count_bar_stack],bottom=s)
                                s+=np.array(m)
                count_bar_stack+=1
	plt.xticks(index,class_name,rotation=13,fontsize=18)
        plt.ylabel(' proportion  of  genus',fontsize=18)
        plt.legend(loc='center left', bbox_to_anchor=(0.0, 1.08),ncol=4,fontsize=15)
        plt.savefig(output+'/6_1_stack_bar_genus.png')
	#Species stack bar
	plt.figure(figsize=(12,11))
	class_name=sorted(species.keys())
	n_groups=len(species.keys())
	index=np.arange(n_groups)
	all_species=sorted(uni_spe.keys())
	spe_arr=sorted(species.keys())
	width=0.35
	mul_spe_arr=[]	
	total_spe={}
	#color start
	color_spe={}
        cn_bar=0
	#print len(all_species)
        for s in all_species:
                color_spe[s]=color_bar[cn_bar]
                cn_bar+=1
	#color end
	for key1 in sorted(species.keys()):
		total=0
		for s in all_species:
			if s not in species[key1]:continue
			else:
				total+=species[key1][s]
		total_spe[key1]=total
	for s in all_species:
		tem=[]
		for key1 in sorted(species.keys()):
			if s not in species[key1]:
				tem.append(0)
			else:
				tem.append(float(species[key1][s])/float(total_spe[key1]))
		mul_spe_arr.append(tem)
	#Pyechart Report Part
        #bar_species=Bar("物种丰度数据堆叠柱状图","Species层面\n\n")
	bar_species=Bar("Relative abundance","Species Level\n\n"
        attr=sorted(bq.keys())
        c=0
        for k in all_species:
                bar_species.add(k,attr,mul_spe_arr[c],is_stack=True,xaxis_rotate=20,yaxis_max=1,legend_top='10%',legend_pos='10%')
                c+=1
        #bar_species.width=1200
        #bar_species.height=600
	g=Grid(width=1200,height=600)
	g.add(bar_species,grid_top='20%')
        g.render(output+'/6_stack_species.html')
	#literature figure part
	count_bar_stack=0
	for m in mul_spe_arr:
		if count_bar_stack==0:
			plt.bar(index,m,width,color=color_spe[all_species[count_bar_stack]],label=all_species[count_bar_stack])	
			tem=m
		else:
			if count_bar_stack==1:
				plt.bar(index,m,width,color=color_spe[all_species[count_bar_stack]],label=all_species[count_bar_stack],bottom=tem)
				s=np.array(m)+np.array(tem)
			else:
				plt.bar(index,m,width,color=color_spe[all_species[count_bar_stack]],label=all_species[count_bar_stack],bottom=s)
                                s+=np.array(m)
		count_bar_stack+=1
	plt.xticks(index,class_name,rotation=13,fontsize=18)
	plt.ylabel(' proportion  of  species',fontsize=18)
	plt.legend(loc='center left', bbox_to_anchor=(-0.1, 1.08),ncol=3,fontsize=15)
        plt.savefig(output+'/6_2_stack_bar_species.png')
	
	

tem={}
for key1 in bq:
        tem=bq[key1]
	#bar(tem,'Plot/'+key1)
	scatter(tem,'Plot/'+key1)
	#upset(key1,'Plot/'+key1)
	#bin_boxplot(tem,'Plot/'+key1)
	#bin_per(tem,'Plot/'+key1)		
	#tax_abundance(tem,'Plot/'+key1)

exit()
ob1=open('Log/br.html','w+')
ob1.write('<br/>\n<br/>\n<br/>\n<br/>\n<br/>\n')
oe=open('Merge_Plot.sh','w+')
ob3=open('Log/scatter_p.html','w+')
ob3.write('<p><a style="font-weight:bold;font-size:18px;">Bin\'s Completeness and Contamination(Scatter)</a></p>\n<br\>\n')
for s in snl:
	if overall=='Y':
		ob2=open('Log/'+s+'_img.html','w+')
		ob2.write('<p><a style="font-weight:bold;font-size:18px;">Overall N50(kb) Bar</a></p>\n<br\>\n')
		ob2.write('<img src="'+s+'/0.n50_overall.png" alt="Overall N50 "/>\n')
		ob2.write('<br/>\n<br/>\n<br/>\n<br/>\n<br/>\n')
	ob5=open('Log/'+s+'_upset.html','w+')
	ob5.write('<p><a style="font-weight:bold;font-size:18px;">Genus Upset Comparison</a></p>\n<br\>\n')
        ob5.write('<img src="'+s+'/4_upset_genus.png" alt="Genus Upset"/>\n')
        ob5.write('<br/>\n<br/>\n<br/>\n<br/>\n<br/>\n')
	ob5.write('<p><a style="font-weight:bold;font-size:18px;">Species Upset Comparison</a></p>\n<br\>\n')
	ob5.write('<img src="'+s+'/4_upset_species.png" alt="Species Upset"/>\n')
	ob5.write('<br/>\n<br/>\n<br/>\n<br/>\n<br/>\n')
	ob6=open('Log/'+s+'_bin_cov_n50.html','w+')
	ob6.write('<style>#bin_coverage_n50{width:1200px;height:600px;}\n#bin_coverage_n50 img{width:1200px; height:600px}\n#overall_bin_coverage_n50{width:1200px;height:300px;}\n#overall_bin_coverage_n50 img{width:1200px; height:300px}\n</style>')
	ob6.write('<p><a style="font-weight:bold;font-size:18px;">Overall Bins Coverage and N50</a></p>\n<br\>\n')
	ob6.write('<div id="overall_bin_coverage_n50">\n<img src="'+s+'/5_bin_cov_n50_overall.png" alt="Overall Bins Coverage and N50"/>\n</div> \n')
        ob6.write('<p><a style="font-weight:bold;font-size:18px;">Bins Coverage and N50</a></p>\n<br\>\n')
	ob6.write('<p><a style="font-weight:bold;font-size:15px;">A-D means High Quality Bins</a></p>\n<br\>\n')
	ob6.write('<p><a style="font-weight:bold;font-size:15px;">E-H means Medium Quality Bins</a></p>\n<br\>\n')
	ob6.write('<p><a style="font-weight:bold;font-size:15px;">I-L means Low Quality Bins</a></p>\n<br\>\n')
        ob6.write('<div id="bin_coverage_n50">\n<img src="'+s+'/5_bin_cov_n50.png" alt="Bins Coverage and N50"/>\n</div> \n')
        ob6.write('<br/>\n<br/>\n<br/>\n<br/>\n<br/>\n')
	if not os.path.exists('Report_'+ero+'/'+s):
		os.makedirs('Report_'+ero+'/'+s,0755)
	if overall=='Y':
		os.system('cp Plot/'+s+'/0.n50_overall.png Report_'+ero+'/'+s+'/0.n50_overall.png')
	os.system('cp Plot/'+s+'/5_bin_cov_n50_overall.png Report_'+ero+'/'+s+'/5_bin_cov_n50_overall.png')
	os.system('cp Plot/'+s+'/5_bin_cov_n50.png Report_'+ero+'/'+s+'/5_bin_cov_n50.png')
	os.popen('cp Plot/'+s+'/4_upset_genus.png  Report_'+ero+'/'+s+'/4_upset_genus.png')
	os.system('cp Plot/'+s+'/4_upset_species.png  Report_'+ero+'/'+s+'/4_upset_species.png')
	if overall=='Y':
		oe.write('cat Log/'+s+'_img.html '+' Plot/'+s+'/1_bar.html Log/br.html '+' Plot/'+s+'/2_lb.html Log/br.html Log/scatter_p.html  Plot/'+s+'/3_scatter.html   Log/'+s+'_upset.html Log/'+s+'_bin_cov_n50.html Plot/'+s+'/6_stack_genus.html Log/br.html Plot/'+s+'/6_stack_species.html  >Report_'+ero+'/'+s+'.html\n')
	else:
		oe.write('cat  Plot/'+s+'/1_bar.html Log/br.html '+' Plot/'+s+'/2_lb.html Log/br.html Log/scatter_p.html  Plot/'+s+'/3_scatter.html   Log/'+s+'_upset.html Log/'+s+'_bin_cov_n50.html Plot/'+s+'/6_stack_genus.html Log/br.html Plot/'+s+'/6_stack_species.html >Report_'+ero+'/'+s+'.html\n')
