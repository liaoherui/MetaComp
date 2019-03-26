import os
import re
import getopt
import sys
import time
opts,args=getopt.getopt(sys.argv[1:],"o:g:")
pwd=os.getcwd()
output=pwd+'/Result'  #[default]
crun='N'   #[default without coverage]
gl='Y'   #[defalut]
for opt,arg in opts:
        if opt=="-o":
                output=arg
        elif opt=='-g':
                gl=arg
if not os.path.exists('list'):
	os.makedirs('list',0755)
if not gl=='N':
        #### Checkm list ####
        os.system('find '+output+' -name \'checkm\' >checkm.list')
        fc=open('checkm.list','r')
        oc=open('list/checkm.list','w+')
        while True:
                line=fc.readline().strip()
                if not line:break
                oc.write(line+'/result'+'\n')
        fc.close()
        os.system('rm checkm.list')
        #### Aragorn list ####
        a=os.popen('find '+output+' -name \'aragorn\' > list/aragorn.list | echo \'aragorn_list done\'').read().strip()
        #### Barrnap list ####
        b=os.popen('find '+output+' -name \'barrnap\' > list/barrnap.list| echo \'barrnap_list done\'').read().strip()
        #### Quast list ####
        q=os.popen('find '+output+' -name \'quast\' > list/quast.list| echo \'quast_list done\'').read().strip()
        #### Kraken list ####
        k=os.popen('find '+output+' -name \'kraken\' > list/kraken.list| echo \'kraken_list done\'').read().strip()
