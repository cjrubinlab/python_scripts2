#This script processes output files from samtools mpileup run on one single chromosome for one single bam and calculates nucleotide diversity in 500kb windows

#usage:
# first run samtools mpileup, e.g.: "samtools mpileup  [process flags] --reference ref.fasta --region chr1  sample.bam > out.chr1.sample.mpileup"

# then to calculate nucleotide diversity estimates:

# "cat out.chr1.sample.mpileup | python WinNuclDivmPileup.py > OUT_nucleotide_diversity.chr1.txt"


## COPYRIGHT 2020 CARL-JOHAN RUBIN AND UPPSALA UNIVERSITY. INQUERIES: carl(dot)rubin(at)gmail.com

import sys
import re
import fileinput
import numpy

# NC_048943.1     2375    t       19      ...................     FJJJFJJJFsJsJJJJJJG
# NC_048943.1     2376    t       19      ...................     JJJJJJJ<sJsJJJJFJJG
# NC_048943.1     2377    t       20      ....................    JJJJAFJJJsJsnJJJJJJI
# NC_048943.1     2378    a       19      ...................     JJJJJJJJFsJssJJJJJI
# NC_048943.1     2379    a       20      ....................    JJJJJAJJ<sJssJJJFJJI
# NC_048943.1     2380    a       20      ....................    JJJJJJJJJsJssJJJFJJI
# NC_048943.1     2381    t       20      ....................    JJJJJJJJJsJssJJJJJJI
# NC_048943.1     2382    a       20      ....................    JJJJJJJJJsJosJJJFJJI
# NC_048943.1     2383    c       20      ....................    JJJJJJJJJsJsoJJJJJJI
# NC_048943.1     2384    t       20      ....................    JJJJJJJJJsJssJJJFJJI
# NC_048943.1     2385    t       20      ....................    JJJJJJJJJsJsoJJJJJJI
# NC_048943.1     2386    t       20      ....................    JJJJFJJJJsJsoJJJJAJI
# NC_048943.1     2387    t       20      ....................    JJJJJJJJJsJjsJJJJFJG
# NC_048943.1     2388    t       20      ....................    JJJJJJJJFsJssJJJJJJI
# NC_048943.1     2389    c       20      ....................    JJJJJJJJ<sJssJJJJJJI
# NC_048943.1     2390    a       20      ....................    JJJJJJJJFsJssJJJJJJI
# NC_048943.1     2391    a       20      ....................    JJJJJFJJJsJssJJJJJJI
# NC_048943.1     2392    a       20      ....................    JJJJJJJJJsJssJJJJJJI
# NC_048943.1     2393    a       20      ....................    JJJJJJJJJsJssJJJFJJI
# NC_048943.1     2394    a       20      ....................    JJJJJJJJJsJssJJJJJJI
# NC_048943.1     2395    t       20      ....................    JJJJJJJJJsJssJJJJJJI
# NC_048943.1     2396    a       20      ....................    JJJJJJJJJsJssJJJJJJG
# NC_048943.1     2397    t       20      ....................    JJJJFJJJJsJssJJJJJJI



def NuclDivFrommPileup():
	nVar=0
	nRef=0
	nSamples=0
	printRef='0'
	printVar='0'
	cols=''
	counterPopulation=0
	counts_array=numpy.zeros([199000000,6],int)
	linecounter=-1
	linecounter2=0
	Arraycounter=-1
	for line in fileinput.input():
		linecounter=linecounter+1	
		linecounter2=linecounter2+1

		sum2_A=0
		sum2_C=0
		sum2_G=0
		sum2_T=0
		sum2_REF=0
		sum2_MostAbundantVar=0
		singleVarCount2=0
		
		
		line=line.replace('\n', '')
		line=line.replace('chr', '')
		line=line.replace(' ', '\t')
		bols = line.split('\t')
		collengths=len(bols)
		chromosome= bols[0]
		bp= int(bols[1])
		if linecounter2==100000:
			linecounter2=0
		nReads=int(bols[3])
		AlleleString=bols[4]
		lenAlleleString=len(AlleleString)
		refAllele=str(bols[2])
		refAllele=refAllele.replace('A','1')
		refAllele=refAllele.replace('a','1')
		refAllele=refAllele.replace('c','2')
		refAllele=refAllele.replace('C','2')
		refAllele=refAllele.replace('G','3')
		refAllele=refAllele.replace('g','3')
		refAllele=refAllele.replace('t','4')
		refAllele=refAllele.replace('T','4')
		refAllele=refAllele.replace('n','5')
		refAllele=refAllele.replace('N','5')

		counterINDEL=0
		counter=-1
		counterA=-1
		counterC=-1
		counterG=-1
		counterT=-1
		counterREF=-1
		Extra=0
		cc=''




		out_str = ''
		totSum=0
		sumA=0
		sumC=0
		sumG=0
		sumT=0
		sumREF=0
		counter=-1
		counterA=0
		counterC=0
		counterG=0
		counterT=0
		counterREF=0
		Extra=0
		out_str=''
		removeNT=''

		for iter in range(lenAlleleString):
			if AlleleString[iter]=='+' or AlleleString[iter]=='-':
				counterINDEL=counterINDEL+1



		for iter in range(lenAlleleString):
			if iter+Extra < lenAlleleString:
				data = AlleleString[iter+Extra]
				if data != '-' and data != '+' and data !='^' and data !='$':
					counter = counter+1
					out_str += data
					if data == 'a' or data == 'A':
						counterA=counterA+1

					if data == 'c' or data == 'C':
						counterC=counterC+1

					if data == 'g' or data == 'G':
						counterG=counterG+1

					if data == 't' or data == 'T':
						counterT=counterT+1

					if data == ',' or data == '.':
						counterREF=counterREF+1

				if data == '^':
					Extra=Extra+1				
				if data == '-' or data == '+':
					if iter <= lenAlleleString - 2:
						removeNT = AlleleString[iter+Extra +1]
						for yy in range(10):
							if AlleleString[iter+Extra +2] == str(yy):
								r=0
								r=str(yy)
								removeNT+= r
						if int(removeNT) < (lenAlleleString-iter+Extra):
							Extra = Extra + int(removeNT) + len(removeNT)				
							


		


		sum2_A=counterA
		sum2_C=counterC
		sum2_G=counterG
		sum2_T=counterT
		sum2_REF=counterREF
		sumall=counterA+counterT+counterC+counterG+counterREF
		sumallvar=counterA+counterT+counterC+counterG
		
		if sumall>=10:
			Arraycounter=Arraycounter+1
			counts_array[Arraycounter,1]=sum2_A
			counts_array[Arraycounter,2]=sum2_C
			counts_array[Arraycounter,3]=sum2_G
			counts_array[Arraycounter,4]=sum2_T
			counts_array[Arraycounter,5]=sum2_REF
			counts_array[Arraycounter,0]=bp

	winsize=500000
	Wins=numpy.zeros([80,5],float)
	CurrWinArray=numpy.zeros([5000000],float)
	q=-1.000
	qq=-1.000
	Het=-100.00
	CurrWinCounter=-1
	CurrWin=0
	for iter in range(80):
		Wins[iter,0]=(iter*(1*winsize))+1
		Wins[iter,1]=(iter+1)*(1*winsize)

	for itr in range(Arraycounter+1):
		data=counts_array[itr,:]

		if  data[0]> Wins[CurrWin,1]:
			Wins[CurrWin,2]=numpy.mean(CurrWinArray[0:CurrWinCounter])
			Wins[CurrWin,3]=CurrWinCounter

			CurrWinCounter=-1
			CurrWinArray=numpy.zeros([5000000],float)
			CurrWin=CurrWin+1
			
		if data[0]>= Wins[CurrWin,0] and data[0]<= Wins[CurrWin,1]:
			CurrWinCounter=CurrWinCounter+1
			if numpy.sum(data[1:5])>0:
				VarCount=numpy.max(data[1:5])
			elif numpy.sum(data[1:5])==0:
				VarCount=0
			
			RefCount=data[5]
			
			q=float(VarCount)/float(RefCount+VarCount)
			qq=float(RefCount)/float(RefCount+VarCount)
			Het=1.00-float(float(q*q)+float(qq*qq))
			CurrWinArray[CurrWinCounter]=float(Het)
		
		if itr==Arraycounter:
			Wins[CurrWin,2]=numpy.mean(CurrWinArray[0:CurrWinCounter])
			Wins[CurrWin,1]=data[0]
			Wins[CurrWin,3]=CurrWinCounter
			
	for ii in range(CurrWin+1):
		print str(chromosome)+'\t'+ str(Wins[ii,0])+'\t'+ str(Wins[ii,1])+'\t'+ str(Wins[ii,2])+'\t'+ str(Wins[ii,3])+'\n',



			

			



if __name__ == '__main__':
    if len(sys.argv) != 1:
        print 'Wrong number of arguments!\n'
        sys.exit(2)
    NuclDivFrommPileup()
