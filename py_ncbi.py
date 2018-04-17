#!/usr/bin/python2.7
##A script to get data from NCBI, parse out taxon numbers and convert to kraken compatible fasta files.

##to install Bio do 'conda install biopython'


from Bio import Entrez
Entrez.email='rharbert@amnh.org'
handle = Entrez.esearch(db="nucleotide", retmax=1000000, term="chloroplast[All Fields] NOT pseudogene[All Fields] NOT unverified[All Fields] AND chloroplast[filter]", idtype="acc")
#handle = Entrez.esearch(db="nucleotide",  term="chloroplast", field="filter", idtype = "acc")
record = Entrez.read(handle)
handle.close()


zz=200
cuts = len(record['IdList'])/zz
taxid = [];
seqs = [];
for n in range(0, cuts):
	z = (n*zz)
	fetch = Entrez.efetch(db='nucleotide', id=",".join(record['IdList'][z:z+zz-1]), rettype='gb', retmode='text')
	fast = fetch.read()
	arfast = fast.split("//\n");
	for pi in arfast:
		if len(pi)>20:
			seqs.append(pi);
			arpi = pi.split("taxon:");
			nexpi = arpi[1].split("\"");
			taxid.append(nexpi[0]);
	

fast = "//\n",join(seqs)
target = open('seqs.gb', 'w')
target.write(fast)
target.close()

from Bio import SeqIO
count = SeqIO.convert("seqs.gb", "genbank", "seqs.fasta", "fasta")
print("Converted %i records" % count)

infile = open('seqs.fasta', 'r')
ffast = infile.read()
infile.close()
arffast = ffast.split("\n>")
n = 0;
arrout = [];
for line in arffast:
	str2 = "|kraken:taxid|%s" % taxid[n]
	str2 = str2 + " "
	ss= line.replace(" ", str2, 1)
	ss = ">" + ss
	arrout.append(ss);
	n=n+1

outstr = "\n\n".join(arrout)
out = open('seqs.fa', 'w')
out.write(outstr)
out.close()
exit()


