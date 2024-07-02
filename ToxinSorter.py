import sys

targets = open(sys.argv[1], 'r+')
CDSseqs = open(sys.argv[2], 'r+')
FullPepSeqs = open(sys.argv[3], 'r+')
MaturePepSeqs = open(sys.argv[4], 'r+')

fp1 = open(sys.argv[2]+"ToxinSorterout.fasta", 'w+')
fp2 = open(sys.argv[3]+"ToxinSorterout.fasta", 'w+')
fp3 = open(sys.argv[4]+"ToxinSorterout.fasta", 'w+')

#print(targets,"\n",CDSseqs,"\n",FullPepSeqs,"\n",MaturePepSeqs)

targetnames = targets.readlines()
CDSlines = CDSseqs.readlines()
FullPeplines = FullPepSeqs.readlines()
MaturePeplines = MaturePepSeqs.readlines()

for name in targetnames:
    for CDSline in CDSlines:
        if CDSline.startswith('>'):  
            if name in CDSline:
                name_found = True
            else:
                name_found = False
        elif name_found:
            #print(">",name,"\n",CDSline.strip())
            fp1.write(">"+name+CDSline+"\n")
    for FullPepline in FullPeplines:
        if FullPepline.startswith('>'):  
            if name in FullPepline:
                name_found = True
            else:
                name_found = False
        elif name_found:
            #print(">",name,"\n",FullPepline.strip())
            fp2.write(">"+name+FullPepline+"\n")
    for MaturePepline in MaturePeplines:
        if MaturePepline.startswith('>'):  
            if name in MaturePepline:
                name_found = True
            else:
                name_found = False
        elif name_found:
            #print(">",name,"\n",MaturePepline.strip())
            fp3.write(">"+name+MaturePepline+"\n")

fp1.close()
fp2.close()
fp3.close()

#with open('PnTx3.fasta.rtf', 'r') as file:
#    lines = file.readlines()
#
#fileName = 'Ctenus_corniger_sra_clean_trinity_Trinity_TRINITY_DN27263_c0_g2_i1'
#
#name_found = False
#
#for line in lines:
#    if line.startswith('>'):  
#        if fileName in line:
#            name_found = True
#        else:
#            name_found = False
#    elif name_found:
#        print(line.strip())
#        break  
