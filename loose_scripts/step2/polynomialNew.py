import re
with open("Final result all NGG and distance 1018.txt") as m: #this is a list with all the NGG 20 nt sequences in each gene from the the ORF catalogue NRRL3
    list1 = m.readlines()

with open("Result heading scoring with distance_second formula.txt") as f: #this is a list of base editor candidates in the CDS catalogue NRRL3. Obtained from altformula.py project
    list2 = f.readlines()

output=open("Both existed sgRNA_polynomial_secondFormula.txt","w") #the overlap between the two files described aboved
Outputlist = list1
for n in list2:
    if re.match("^>", n):
        title = n.strip()
        # print(title)
        output.write("\n" + title + "\n")
    else:
        for j in list1:
            if j[0:21] == n[5:26]:
                #print("5'-TGAGGACGAAACGAGTAAGCTCGTC" + n[-24:-4] + "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG-3" )
                output.write(n)
                # print(n + "\n")
output.close()
