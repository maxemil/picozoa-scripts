from Bio import SeqIO
import os
import sys

hmmin = sys.argv[1]
fastain = sys.argv[2]
fastaout = sys.argv[3]

best_score = 0
hits = []
for line in open(hmmin):
  if not line.startswith('#'):
      line = line.split()
      if float(line[5]) > best_score:
          best_score = float(line[5])
      if float(line[5]) >= best_score * .7:
        hits.append(line[0])
        

with open(fastaout, 'w') as out:
  for rec in SeqIO.parse(fastain, 'fasta'):
      if rec.id in hits:
          SeqIO.write(rec, out, 'fasta')
