#the final script from Jun

1.
import re
with open("Both existed sgRNA_polynomial_FirstFormula.txt") as f: #output from the ORF vs CDS overlap
    list1 = f.readlines()

output = open("2.Candidates above 0.6.txt", "w")

for i in list1:
    if re.match("^>", i):
        title = i.strip()
        output.write('\n' + title + "\n")
    else:
        if i != '\n':
            if float(i[0:4]) >= 0.7:  #cut-off
                output.write(i)
output.close()


2.
import re
file = open("2.Candidates above 0.6.txt")
lines =[]
for l in file.readlines():
    lines.append(l)

list1 = []
list2 = []
for i in lines:
    list1.append(i)

output=open("3.Sorted sgRNA list.txt","w")
n = 0
for k in list1:
    list3 = []
    if re.match("^>", k):
        title = k
        # print(k)
        output.write("\n" + k)
    else:
        if k != "\n":
            list2.append(k.strip())
        else:
            list3 = sorted(list2, key=lambda x:x[0:4], reverse=True)
            list2 = []
            # print("\n".join(list3))
            # print(str(list3))
            output.write("\n".join(list3))
output.close()

3.
import re
file = open("3.Sorted sgRNA list.txt")

lines =[]
for l in file.readlines():
    lines.append(l)

output=open("4.Combine sense and antisense.txt","w")
Outputlist = lines
for i in lines:
    if "Antisense" not in i:
        output.write(i)
output.close()

3.
import re
file = open("4.Combine sense and antisense.txt")

lines =[]
for l in file.readlines():
    lines.append(l)

output=open("5.line Combine sense and antisense.txt","w")
for i in lines:
    if re.match("^>", i):
        title = i.strip()
        output.write("\n" + title + "\n")
    else:
        output.write(i)
output.close()

4.
import re
file = open("5.line Combine sense and antisense.txt")
lines =[]
for l in file.readlines():
    lines.append(l)

list1 = []
list2 = []
for i in lines:
    list1.append(i)

output=open("6.Sorted sgRNA list.txt","w")
n = 0
for k in list1:
    list3 = []
    if re.match("^>", k):
        title = k
        # print(k)
        output.write("\n" + "\n"+ k)
    else:
        if k != "\n":
            list2.append(k.strip())
        else:
            list3 = sorted(list2, key=lambda x:x[0:4], reverse=True)
            list2 = []
            # print("\n".join(list3))
            # print(str(list3))
            output.write("\n".join(list3))
output.close()

5.
import re

with open("6.Sorted sgRNA list.txt") as f:
    list1 = f.readlines()

list2 = []
Str = []
output=open("7.keep highest sgRNA library.txt","w")

for i in list1:
    if re.match("^>", i):
        list2 = []
        title = i.strip()
        # print(title)
        output.write("\n" + title + "\n")
    else:
        if i != "\n":
            list2.append(i.strip())
        else:
            Str = ",".join(list2)
            Str1 = Str[0:33]
            output.write(Str1 + "\n")
output.close()

6.
import re

with open("7.keep highest sgRNA library.txt") as f:
    list1 = f.readlines()

output=open("8.final sgRNA library for ordering.txt","w")

for i in list1:
    if re.match("^>", i):
        title = i.strip()
        # print(title)
    else:
        if i != "\n":
            output.write(i)
output.close()

7.
with open("8.final sgRNA library for ordering.txt") as f:
    list1 = f.readlines()

output = open("9.distributionATG_distance.txt","w")
for i in list1:
    output.write(i[-5:]) #to get the the final list for ordering fix the slicing window and use add or maybe append to add the cloning flanks
