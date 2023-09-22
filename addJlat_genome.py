from Bio import SeqIO
from Bio.Seq import Seq
from os import path

###############################################################################
# constants
###############################################################################

# fasta file from MN989412 which is 10.6 provirus
MN989412_FA_FN = "fasta/jlat10-6.orig.fa"

JLAT_SEQUENCE_FN = "fasta/jlat.cleaned.fa"
HG38_GENOME_FN = "fasta/hg38_only_genome.fa"

HG38_DENYLIST_FN = "denylist/hg38_boyleLab_denylist.v2.bed"

# take`n from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5237276/
# position is the hg38 nucleotide position following 3' LTR end.
# note that the above paper has incorrect numbers in their table
JLAT_META = {
  "10.6": ["chr9", "+", 136468579],
  "8.4": ["chr1", "-", 77946384],
  "9.2": ["chr19", "+", 46381104],
  "15.4": ["chr19", "+", 34441293],
  "6.3": ["chr19", "+", 46381104]
}


###############################################################################
# basic functions 
###############################################################################
def get_jlat_only_sequence(fn: str) -> Seq:
  """Save and return JLat proviral sequence.
  
  This function relies on sequencing performed in Chung CH, Mele AR, Allen AG,
  Costello R, Dampier W, Nonnemacher MR, Wigdahl B. Integrated Human 
  Immunodeficiency Virus Type 1 Sequence in J-Lat 10.6. Microbiol Resour 
  Announc. 2020 Apr 30;9(18):e00179-20. doi: 10.1128/MRA.00179-20. PMID: 
  32354973; PMCID: PMC7193928. Fasta file can be downloaded from GenBank under
  accession number: MN989412. Since this sequence include the integration
  site junction, the host sequence is trimmed off and only the proviral
  sequence is returned.
  
  Args:
    fn: filename for the fasta record from the above publication.

  Returns:
    A Seq object containing only the entire proviral sequence.
  """

  # read jlat fasta and clean
  with open(fn) as handle:
    for record in SeqIO.parse(handle, "fasta"):
      jlat_cleaned = record
      
  jlat_cleaned.seq = jlat_cleaned.seq[145:10345]
  jlat_cleaned.id = "jlat_proviral_sequence"
  jlat_cleaned.description = ""

  with open(JLAT_SEQUENCE_FN, "w") as handle:
    SeqIO.write(jlat_cleaned, handle, "fasta")

  return(jlat_cleaned.seq)

def generate_integrated_jlat(jlat_seq: Seq, target_dir: str, jlat_name: str, jlat_type: list) -> None:
  """todo will add this stuff
  
  Args:
    jlat_seq:
    tardet_dir:
    jlat_name:
    jlat_type:

  Returns:
    None

  """
  chrName = jlat_type[0]
  orient = jlat_type[1]
  pos = jlat_type[2]

  proviralSize = len(jlat_seq)
  targetFn = path.join(target_dir, jlat_name + ".fa")
  
  # some notes...
  # for 10.6, the recorded integration site (host nucleotide that is proximal to 3' LTR end)
  # is 136468579. However, that's because there's a 5bp duplication at the site of integration
  # so for 10.6, which is a positive oriented provirus, the proximal point at the 5' site is
  # actually 136468583 since the -4bp as well +0 (136468583; which means a total of 5bp) are 
  # duplicated at the 3' end...hence why the 3' end integration site is marked as 136468579
  
  # For all positive oriented proviruses, the integration site in the meta constant
  # needs to be subtracted by 4 for the 5' insertion site.
  
  # For all negative oriented proviruses, the integration site in the meta constant
  # is the same for the host genome oriented 5' insertion site. the reverse complement of
  # the provirus is then added. Next, the 5bp duplication occurs (5bp adjacent to where
  # reverse complement provirus was just entered)
  #insertionPos = 136468583


  with open(HG38_GENOME_FN, "r") as in_handle:
    with open(targetFn, "w") as out_handle:
      for record in SeqIO.parse(in_handle, "fasta"):
        if record.id == chrName and orient == "+":
          insertionPos = pos + 4

          newChr = record.seq[:insertionPos] # reminder that seq is 0 based
          newChr += jlat_seq
          newChr += record.seq[insertionPos - 5:]
          
          record.seq = newChr
          
          print(record.seq[insertionPos:insertionPos+50])
          print(record.seq[(insertionPos + proviralSize): (insertionPos + proviralSize + 50)])

        elif record.id == chrName and orient == "-":
          insertionPos = pos + 1 # because we are working in 0-based system

          newChr = record.seq[:insertionPos]
          newChr += jlat_seq.reverse_complement()
          newChr += record.seq[insertionPos - 5:]
          
          record.seq = newChr

          print(record.seq[insertionPos:insertionPos+50])
          print(record.seq[(insertionPos + proviralSize): (insertionPos + proviralSize + 50)])

        # write out chr
        SeqIO.write(record, out_handle, "fasta")


def generate_shifted_denylist(jlat_name: str, jlat: list, jlat_len: int = 10200, duplication: int = 5) -> None: 
  """todo will add this stuff
  
  Args:
    jlat_name:
    jlat:
    jlat_len:
    duplication:

  Returns:
    None

  """
  chrom = jlat[0]
  orient = jlat[1]
  pos = jlat[2]
  
  # reminder that bed file is 0 indexed. start is inclusive. end is exclusive.
  if orient == "-":
    insertionPos = pos + 1
  elif orient == "+":
    insertionPos = pos + 4
  
  modified_fn = "denylist/" + jlat_name + ".bed"

  with open(HG38_DENYLIST_FN, "r") as f:
    with open(modified_fn, "w") as fo:
      for line in f:
        fields = line.split("\t")

        if fields[0] == chrom and int(fields[1]) > insertionPos:
          print("updating {}".format(str(fields)))

          fields[1] = int(fields[1]) + jlat_len + duplication
          fields[2] = int(fields[2]) + jlat_len + duplication

        elif fields[0] == chrom and int(fields[2]) > insertionPos:
          print("updating {}".format(str(fields)))
          
          fields[2] = int(fields[2]) + jlat_len + duplication

      fo.write("\t".join(map(lambda x: str(x), fields)))


def main():
  jlat_seq = get_jlat_only_sequence(MN989412_FA_FN)
  
  for jlat in JLAT_META:
    print("Working on jlat:", jlat)
    generate_integrated_jlat(jlat_seq, "fasta", jlat, JLAT_META[jlat])
  
    generate_shifted_denylist(jlat, JLAT_META[jlat])


if __name__ == "__main__":
  main()


def make_integrated_gtf():
  # integrate jlat to hg38 gtf
  with open("hg38_with_jlat10-6.unprocesesed.gtf", "r") as infile:
    with open("hg38_with_jlat10-6.procesesed.gtf", "w") as outfile:
      counter = 0
      for line in infile:
        if counter < 5:
          outfile.write(line)
          counter += 1
          continue

        fields = line.split("\t")
        fields[3] = int(fields[3])
        fields[4] = int(fields[4])

        if fields[0] == "chr9" and fields[3] > insertionPos:
          fields[3] += proviralSize + 5 # account for 5bp duplication

        if fields[0] == "chr9" and fields[4] > insertionPos:
          fields[4] += proviralSize + 5 # account for 5bp duplication

        if fields[0] == "chrHIVJLat10.6":
          fields[0] = "chr9"
          fields[3] += insertionPos
          fields[4] += insertionPos

          print(str(fields))

        outfile.write("\t".join(map(lambda x: str(x), fields)))

        counter += 1
