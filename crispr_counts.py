
import sys
import os
from collections import defaultdict

#key is sequence
def create_library_dict(filename=None,length_dict=None,dup_seqs=None):
  d={}
  seq_array=[]
  with open(filename) as f:
    for line in f:
      sgrna = line.split(",")
      seq_array.append(sgrna[1])
      seq_length=len(sgrna[1])
      length_dict.update({seq_length:0})
      d[sgrna[1]]= sgrna
  #dup_seqs=dup_seq(seq_array)
  #flagged_seq=check_seq(seq_array)
  return d


def check_seq(sequence_array=None):
  seq_score=[]
  for i in sequence_array:
    #print(i)
    max_sim_score=calculate_seq_score(i,[x for x in sequence_array if x not in i])
    seq_score.append(max_sim_score)
  return seq_score

#key is read name
def process_fastq(filename=None,input_dict=None,n=0,seq_length_counter=None,dup_seq=None,all_length_counter=None,correct_len_seq_counter=None):
  unmappedreads=0
  numreads=0
  nummulti=0
  num_wronglength=0
  num_mapped=0
  with open(filename, 'r') as fh:
    lines = []
    for line in fh:
      lines.append(line.rstrip())
      if len(lines) == n:
        numreads+=1
        sequence = lines[1]
        if sequence in dup_seq:
          nummulti+=1
        elif sequence in input_dict:
          input_dict[sequence]+=1
          num_mapped+=1
        else:
          unmappedreads+=1
        if len(sequence) in seq_length_counter:
          seq_length_counter[len(sequence)]+=1
          if sequence in correct_len_seq_counter.keys():
            if sequence not in input_dict:
              correct_len_seq_counter[sequence]+=1
          else:
            if sequence not in input_dict:
              correct_len_seq_counter.update({sequence:1})
          if len(sequence) in all_length_counter.keys(): 
            all_length_counter[len(sequence)]+=1
          else:
            all_length_counter.update({len(sequence):1})
        else:
          num_wronglength+=1
          if len(sequence) in all_length_counter.keys():
            all_length_counter[len(sequence)]+=1
          else:
            all_length_counter.update({len(sequence):1})
        lines=[]
  return [unmappedreads,numreads,nummulti,num_wronglength,num_mapped]

try:
    fastq = sys.argv[1]
except IndexError as ie:
    raise SystemError("Error: Specify fastq file name\n")

if not os.path.exists(fastq):
    raise SystemError("Error: fastq file does not exist\n")

try:
    lib_file = sys.argv[2]
except IndexError as ie:
    raise SystemError("Error: Specify library file name\n")

if not os.path.exists(lib_file):
    raise SystemError("Error: lib_file file does not exist\n")

try:
    statistics_outfile = sys.argv[3]
except IndexError as ie:
    raise SystemError("Error: Specify statistics out file name\n")


try:
    count_outfile = sys.argv[4]
except IndexError as ie:
    raise SystemError("Error: Specify count out file name\n")

try:
    correct_len_seq_outfile = sys.argv[5]
except IndexError as ie:
    raise SystemError("Error: Specify correct length sequence out file name\n")

try:
    all_read_length_outfile = sys.argv[6]
except IndexError as ie:
    raise SystemError("Error: Specify all length out file name\n")


#number of lines for a fastq file
#consider if fastq has funny business
record_lines = 4
#sample_reads =  process_fastq(fastq)
#keeps track of sgRNA count
reads_length_dict={}
#keeps track of duplicated sgRNA counts
dupe_sgnra_seq=()
#histogram for other reads lengths
all_reads_length_dict={}
#keeps track of sequence for reads of correct length
correct_length_seq_dict={}
sgrna_dict=create_library_dict(lib_file,reads_length_dict,dupe_sgnra_seq)
sgrna_seq = list(sgrna_dict.keys())
sample_dict = {key: 0 for key in sgrna_seq}

reads_summary=process_fastq(fastq,sample_dict,record_lines,reads_length_dict,dupe_sgnra_seq,all_reads_length_dict,correct_length_seq_dict)

#output_filename=fastq[0:fastq.index('.')]

with open(statistics_outfile, 'w') as f:
    f.write('Number_of_processed_reads\tNum_mapped_reads\tNumber_of_unmapped_reads\tNumber_of_multimapped_reads')
    seq_lens=reads_length_dict.keys()
    for i in seq_lens:
        f.write('\tNumber of reads with length ' + str(i))
    f.write('\tNumber reads with different length\n')
    f.write(str(reads_summary[1]) + '\t' + str(reads_summary[4]) + '\t' + str(reads_summary[0]) + '\t' + str(reads_summary[2]))
    for i in seq_lens:
        f.write('\t'+str(reads_length_dict[i]))
    f.write('\t'+str(reads_summary[3])+'\n')
f.close()


with open(correct_len_seq_outfile, 'w') as f:
    f.write('Sequence\tCounts\n')
    seqs=correct_length_seq_dict.keys()
    for i in seqs:
        f.write(str(i) + '\t' + str(correct_length_seq_dict[i]) + '\n')
f.close()

with open(all_read_length_outfile, 'w') as f:
    f.write('Length\tCounts\n')
    len=all_reads_length_dict.keys()
    for i in len:
        f.write(str(i) + '\t' + str(all_reads_length_dict[i]) + '\n')
f.close()

ds = [sgrna_dict, sample_dict]
sample_counter = defaultdict(list)

for d in (sgrna_dict,sample_dict): 
    for key, value in d.items():
        sample_counter[key].append(value)

df=open(count_outfile,"w")
for i in sgrna_seq:
  df.write(sample_counter[i][0][0] + "\t" + str(sample_counter[i][1]) + "\n")
df.close()


