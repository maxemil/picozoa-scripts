params.faa = ""
params.outdir = "results"

faa_input = Channel.fromPath(params.faa)

process alignSequences {
  input:
  file faa from faa_input

  output:
  file "${faa.simpleName}.mafft" into mafft_alignments
  file "${faa.simpleName}.aln" into filtered_alignments

  tag {"${faa.simpleName}"}
  publishDir "$params.outdir/alignments", mode: 'copy'

  script:
  """
  sed 's/*//g' $faa > ${faa.simpleName}.clean
  sed -i 's/@/../g' ${faa.simpleName}.clean
  mafft-einsi --thread ${task.cpus} ${faa.simpleName}.clean > ${faa.simpleName}.mafft
  trimal -in ${faa.simpleName}.mafft -out ${faa.simpleName}.aln -gt 0.01
  """
}

process removeShortSeqs {
  input:
  file aln from filtered_alignments

  output:
  file "${aln.simpleName}.noshort.aln" into alignments_noshort
  
  tag {"${aln.simpleName}"}
  publishDir "$params.outdir/alignments", mode: 'copy'

  script:
  """
  #! /usr/bin/env python3
  from Bio import SeqIO
  from numpy import median
  
  threshold = 0.1
  lengths = []
  
  for rec in SeqIO.parse("$aln", 'fasta'):
      lengths.append(len(rec.seq.ungap()))
      
  with open("${aln.simpleName}.noshort.aln", 'w') as out:
    for rec in SeqIO.parse("$aln", 'fasta'):
        if len(rec.seq.ungap()) > threshold * median(lengths):
            SeqIO.write(rec, out, 'fasta')
  """
}


process makeTrees {
  input:
  file aln from alignments_noshort

  output:
  file "${aln.simpleName}.*" into trees
  file "${aln.simpleName}.treefile" into treefiles

  tag {"${aln.simpleName}"}
  publishDir "$params.outdir/trees", mode: 'copy'

  script:
  """
  iqtree -s $aln -m LG+G+F -nt ${task.cpus} -fast -pre ${aln.simpleName} -keep-ident
  """
}
