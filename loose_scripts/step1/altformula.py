import re
#this script gives a list of base editing gRNAs with an score
file = open("NewData including antisense.txt")
lines =[]
for l in file.readlines():
    lines.append(l)

dic={}
for i in lines:
    if re.match("^>",i):
        title = i.strip()
        dic[title]=""
    else:
        dic[title] = dic[title] + i.strip("\n")

fif = []
output=open("Result heading scoring with distance.txt","w")
for k in dic:
    res=re.findall("(?=(.{20}.GG))",dic[k])
    if res:
        output.write("\n" + k+"\n")
        for i in res:
            PD1 = 0
            PA1 = 0
            PB1 = 0
            PC1 = 0
            if "Antisense" in k:
                In = len(dic[k]) - dic[k].index(i)
                POS = (In - 6)/len(dic[k])
                if "CCA" in i[0:11]:
                    PD = In - i[0:11].index("CCA")
                    if PD % 3 == 0:
                        PD1 = PD
                    else:
                        PD1 = 0
                else:
                    PD1 = 0
            else:
                In = dic[k].index(i)
                POS = (In + 17)/len(dic[k])
                if "CGA" in i[1:11]:
                    PA = In + i[1:11].index("CGA") + 4
                    if PA % 3 == 0:
                        PA1 = PA
                    else:
                        PA1 = 0
                else:
                    PA1 = 0
                if "CAA" in i[1:11]:
                    PB = In + i[1:11].index("CAA") + 4
                    if PB % 3 == 0:
                        PB1 = PB
                    else:
                        PB1 = 0
                else:
                    PB1 = 0
                if "CAG" in i[1:11]:
                    PC = In + i[1:11].index("CAG") + 4
                    if PC % 3 == 0:
                        PC1 = PC
                    else:
                        PC1 = 0
                else:
                    PC1 = 0
            if POS <= 0.50:  #position formula
                Sd = 2 * POS
            else:
                Sd= 5.1535*(POS*POS)-(9.6508*POS)+4.5061 #Sd = 2-2 * POS
            filters = ["AAAA","TTTT","CCCC","GGGG"] #homopolymer penalty component
            if any(name in i for name in filters) == False:
                shp = 1
            else:
                shp = 0
            Sub = i[:20]
            GC = (Sub.count('G') + Sub.count('C'))/ 20
            if GC < 0.60: #GC content formula
                Sgc = GC * 2 - 0.20
            else:
                Sgc = 2.53 - GC * 2.55
            lj = format((Sd * 0.75 + Sgc * 0.13 + shp * 0.12), '.2f')
            if PD1 !=0 or PA1 != 0 or PB1 != 0 or PC1 != 0 or PD1 != 0:
                fif = str(lj) + "," +  i + ","+ format((POS),'.2f') + "\n"
                output.write(str(fif))
output.close()
