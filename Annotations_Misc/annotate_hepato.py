import sys
import pandas as pd

inmark = sys.argv[1]
outmark=sys.argv[2]

### save markers gene list as markers.txt  ###
markers1 = pd.read_csv("/Users/akhaliq/Desktop/MarkerGenesListCode/markers/hepato_markers.txt",sep="\t",encoding= 'unicode_escape')
markers1.to_excel(outmark, sheet_name = "markers")


d={}
for i in open("/Users/akhaliq/Desktop/MarkerGenesListCode/markers/hepato_markers.txt",encoding= 'unicode_escape'):
	line1=i.strip().split("\t")
	if len(line1) == 2:
		#w,code=line.split('\t')
		w=line1[0].upper()
		code=line1[1]
		d.setdefault(w,[]).append(code)


e = {i:str(j) for i,j in d.items() }
print(e)
markers = pd.Series(e).to_frame('subtype')
markers.reset_index(inplace=True)
markers = markers.rename(columns = {'index':'gene'})
print(markers)

### save output from FindMarkers as marker_genes_list.csv ###
infile=pd.read_csv(inmark,sep =",",header=0)
rnum = int(infile.iloc[-1,-2]) + 1
print(rnum)
### set number of clusters ###
for i in range(1,rnum):
	j=str(i)
	a = infile[infile["cluster"] == i]
	b = a.merge(markers, on = "gene" , how = "left")
	b.sort_values('avg_log2FC', ascending=False, inplace = True)
	print(b)
	with pd.ExcelWriter(outmark, engine="openpyxl" ,mode="a") as f:
		b.to_excel(f, sheet_name="cluster_"+str(i) )
    
"""
xl = pd.ExcelFile('marker_genes_all.xlsx')
list1=xl.sheet_names[1:]

markers1= pd.read_excel(open('marker_genes_all.xlsx', 'rb'),sheet_name='markers')


markers1.to_excel("marker_genes_annotated.xlsx", sheet_name = "markers")

for i in list1:
	a = pd.read_excel(open('marker_genes_all.xlsx', 'rb'),sheet_name=i)
	print(a)
	b = a.merge(markers, on = "gene" , how = "left")
	b.sort_values('avg_log2FC', ascending=False, inplace = True)

	with pd.ExcelWriter("marker_genes_annotated.xlsx", engine="openpyxl", mode="a") as f:
		b.to_excel(f, sheet_name=i )

"""
