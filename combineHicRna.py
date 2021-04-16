
import sys
import re
import pathlib as pl
from distutils.version import LooseVersion

# ####### #
# classes #
# ####### #
class DefaultSortOrderKey:
	counter = 0

	def key_func(*args):
		DefaultSortOrderKey.counter += 1
		return DefaultSortOrderKey.counter

class Agp2ScaffException(Exception):
	pass

class Gap:
	def __init__(self, length=0, evidence="na"):
		self.length = length
		self.evidence = evidence

	def setLength(self, length):
		self.length = length

	def getLength(self):
		return self.length

	def setEvidence(self, evidence):
		self.evidence = evidence

	def getEvidence(self):
		return self.evidence

	def getSeq(self):
		return 'N' * self.length

	def getOrientedSeq(self):
		return self.getSeq()
	
	def flipOrientation(self):
		pass

	def getAgpRecordTail(self):
		raise NotImplementedError

	def __repr__(self):
		return f"length: {self.length}, evidence: {self.evidence}"

	def __str__(self):
		return f"Gap: {{{self.__repr__()}}}"

class GapKnownSize(Gap):
	def __parseAgpRecord__(self, record, fullAgpRecord=False):
		if fullAgpRecord:
			fields = record.rstrip('\n').split('\t')
			record_type = fields[4]
			length = int(fields[5])
			evidence= fields[8]
			if record_type != 'N':
				raise Agp2ScaffException("AGP records not labelled as type N don't belong as gaps of known size\n")
		else:
			length = int(record[0])
			evidence = fields[3]
		return length, evidence

	def __init__(self, length=100, evidence="na", agp_record=None, fullAgpRecord=False): # agp_record overrides length
		self.length = length
		self.evidence = evidence
		if not agp_record is None:
			self.length, self.evidence = self.__parseAgpRecord__(agp_record, fullAgpRecord=fullAgpRecord)

	def deepCopy(self):
		return GapKnownSize(length=self.length, evidence=self.evidence)

	def getAgpRecordTail(self):
		return ("N", self.length, "scaffold", "yes", self.evidence)

	def __str__(self):
		return f"GapKnownSize: {{{self.__repr__()}}}"

class GapUnknownSize(Gap):
	def __parseAgpRecord__(self, record, fullAgpRecord=False):
		length = 0
		evidence = "na"
		fields = record.rstrip('\n').split('\t')
		if fullAgpRecord:
			record_type = fields[4]
			length = int(fields[5])
			evidence= fields[8]
			if record_type != 'U':
				raise Agp2ScaffException("AGP records not labelled as type U don't belong as gaps of unknown size\n")
		else:
			length = int(fields[0])
			evidence = fields[3]
		return length, evidence

	def __init__(self, length=100, evidence="na", agp_record=None, fullAgpRecord=False): # agp_record overrides length
		self.length = length
		self.evidence = evidence
		if not agp_record is None:
			self.length, self.evidence = self.__parseAgpRecord__(agp_record, fullAgpRecord=fullAgpRecord)

	def deepCopy(self):
		return GapUnknownSize(length=self.length, evidence=self.evidence)

	def getAgpRecordTail(self):
		return ("U", self.length, "scaffold", "yes", self.evidence)

	def __str__(self):
		return f"GapUnknownSize: {{{self.__repr__()}}}"

class Contig:
	__comp__ = {"A": "T", "C": "G", "G": "C", "T": "A",
					"a": "t", "c": "g", "g": "c", "t": "a",
					"N": "N", "n": "n", "R": "N", "r": "n",
					"Y": "N", "y": "n", "S": "N", "s": "n",
					"W": "N", "w": "n", "K": "N", "k": "n",
					"M": "N", "m": "n", "B": "N", "b": "n",
					"D": "N", "d": "n", "H": "N", "h": "n",
					"V": "N", "v": "n" }
	
	def __parseAgpRecord__(self, record, fullAgpRecord=False):
		header = None
		length = None
		orient = None
		fields = record.rstrip('\n').split('\t')
		if fullAgpRecord:
			record_type = fields[4]
			header = fields[5]
			if record_type != 'W':
				raise Agp2ScaffException("AGP records not labelled as type W don't belong as Contigs\n")
			#start = int(fields[6])
			length = int(fields[7])
			orient = fields[8]
		else:
			header = fields[0]
			#start = int(fields[1])
			length = int(fields[2])
			orient = fields[3]
		return header, length, orient


	def __init__(self, header='', length='0', orientation='+', seq='', agp_record=None, fasta=None, fullAgpRecord=False): # agp_record overrides header, length, and orientation. fasta overrides seq.
		if not agp_record is None:
			self.header, self.length, self.orientation = self.__parseAgpRecord__(agp_record, fullAgpRecord=fullAgpRecord)
		else:
			self.header = header
			self.length = length
			self.orientation = orientation

		if not fasta is None:
			self.seq = fasta[self.header]
		else:
			self.seq = seq

	def getHeader(self):
		return self.header
	
	def getName(self):
		return self.header

	def getLength(self):
		return self.length

	def getOrientation(self):
		return self.orientation
	
	def setOrientation(self, orientation):
		self.orientation = orientation

	def flipOrientation(self):
		self.orientation = '+' if self.orientation == '-' else '-'

	def deepCopy(self):
		return Contig(header=self.header, length=self.length, orientation=self.orientation, seq=self.seq)

	def getSeq(self):
		return self.seq

	def __revComp__(self):
		return ''.join(Contig.__comp__[nuc] for nuc in self.seq[::-1])

	def getOrientedSeq(self):
		return self.__revComp__() if self.orientation == '-' else self.seq

	def getAgpRecordTail(self):
		return ("W", self.header, 1, self.length, self.orientation)
	
	def isEqualByName(self, other):
		return self.header == other.header

	def __repr__(self):
		seq = f"{self.seq[:3]}..{self.seq[-3:]}" if self.length > 6 else self.seq
		return f"header: {self.header}, length: {self.length}, orientation: {self.orientation}, seq: {seq}"

	def __str__(self):
		return f"Contig: {{{self.__repr__()}}}\n"

class Scaffold:

	def complementOrientations(orientations):
		comp = []
		for orientation in orientations:
			c = '+' if orientation == '-' else '-'
			comp.append(c)
		return comp
	
	def __twoListsAreEqual__(l1, l2):
		return len(l1) == len(l2) and all(map(lambda x, y: x == y, l1, l2))

	def __init__(self, name=''):
		self.name = name
		self.num_components = 0
		self.num_gaps = 0
		self.num_contigs = 0
		self.num_bases = 0
		self.num_bases_known = 0
		self.num_bases_unknown = 0
		self.components = []
	
	def __init__(self, name='', components=[]):
		self.name = name
		self.num_components = 0
		self.num_gaps = 0
		self.num_contigs = 0
		self.num_bases = 0
		self.num_bases_known = 0
		self.num_bases_unknown = 0
		self.components = []
		for component in components:
			self.components.append(component.deepCopy())
			if isinstance(component, Gap):
				self.num_gaps += 1
				self.num_bases_unknown += component.getLength()
			elif isinstance(component, Contig):
				self.num_contigs += 1
				self.num_bases_known += component.getLength()
			self.num_components += 1
			self.num_bases += component.getLength()

	def getName(self):
		return self.name

	def setName(self, name):
		self.name = name
	
	def addContig(self, contig):
		self.components.append(contig)
		self.num_components += 1
		self.num_contigs += 1
		self.num_bases += contig.length
		self.num_bases_known += contig.length

	def addGap(self, gap):
		self.components.append(gap)
		self.num_components += 1
		self.num_gaps += 1
		self.num_bases_unknown += gap.length

	def addComponent(self, component_type, agp_record_tail, fasta_dict):
		if component_type == 'W':
			self.addContig( Contig(agp_record=agp_record_tail, fasta=fasta_dict) )
		elif component_type == 'N':
			self.addGap(GapKnownSize(agp_record=agp_record_tail))
		elif component_type == 'U':
			self.addGap(GapUnknownSize(agp_record=agp_record_tail))
		else:
			raise Agp2ScaffException("ERROR: Expected the component type ({component_type}) to be 'N', 'U', or 'W'.")

	def deepCopy(self):
		s = Scaffold(name=self.name)
		for component in self.components:
			c = component.deepCopy()
			if isinstance(c, Contig):
				s.addContig(c)
			else: #if isinstance(c, GapKnownSize) or isinstance(c, GapUnknownSize):
				s.addGap(c)
		return s
	
	def getContigHeadsAndOrients(self):
		heads = []
		orients = []

		for component in self.components:
			if isinstance(component, Contig):
				heads.append(component.header)
				orients.append(component.orientation)

		return heads, orients
	
	def getComponents(self):
		comps = []
		for comp in self.components:
			#comps.append(comp) # shallow copy
			comps.append(comp.deepCopy())
		return comps

	def getContigs(self):
		contigs = []
		for component in self.components:
			if isinstance(component, Contig):
				contigs.append(component.deepCopy())
		return contigs
	
	def getContigHeaders(self):
		return [component.header for component in self.components if isinstance(component, Contig)]
	
	def getContigOrientations(self):
		return [component.orientation for component in self.components if isinstance(component, Contig)]
	
	def isEqualByComponents(self, other):
		my_headers = self.getContigHeaders()
		my_orients = self.getContigOrientations()
		other_headers = other.getContigHeaders()
		other_orients = other.getContigOrientations()

		return ( ( Scaffold.__twoListsAreEqual__(my_headers, other_headers) and Scaffold.__twoListsAreEqual__(my_orients, other_orients) ) or 
			( Scaffold.__twoListsAreEqual__(my_headers, other_headers[::-1]) and Scaffold.__twoListsAreEqual__(my_orients, Scaffold.complementOrientations(other_orients[::-1])) ) )
	
	def containsContig(self, contig):
		for component in self.components:
			if isinstance(component, Contig):
				if component.isEqualByName(contig):
					return True
		return False

	def setFixedUnknownGapLength(self, length):
		for component in self.components:
			if isinstance(component, GapUnknownSize):
				component.setLength(length)
	
	def supplantFrom(self, other_scaffolds):
		replacement_components = []
		for component in self.components:
			if isinstance(component, Gap):
				replacement_components.append(component.deepCopy())
			elif isinstance(component, Contig):
				head = component.getHeader()
				scaff = None
				for other_scaff in other_scaffolds:
					if head == other_scaff.getName():
						scaff = other_scaff
						break
				assert not scaff is None
				comps = scaff.getComponents()
				if component.getOrientation() == '-':
					for comp in comps:
						comp.flipOrientation()
					comps = comps[::-1]
				replacement_components.extend(comps)

		return Scaffold(name=self.name, components=replacement_components)
	
	def getSeq(self):
		return ''.join(component.getOrientedSeq() for component in self.components)

	def getFasta(self):
		return f">{self.name}\n{self.getSeq()}\n"

	def getAgp(self):
		agp_records = []
		counter = 0
		object_coord = 1
		if len(self.components) == 1:
			self.components[0].setOrientation('+')
		for component in self.components:
			counter += 1
			agp_record = [self.name, object_coord]
			object_coord += component.getLength()
			agp_record.extend((object_coord - 1,counter))
			agp_record.extend(component.getAgpRecordTail())
			agp_record = '\t'.join(map(str, agp_record))
			agp_records.append(agp_record)
		return '\n'.join(agp_records) + '\n'

	def __repr__(self):
		s = f"name: {self.name}, num_components: {self.num_components}, num_gaps: {self.num_gaps}, num_contigs: {self.num_contigs}, num_bases: {self.num_bases}, num_bases_known: {self.num_bases_known}, num_bases_unknown: {self.num_bases_unknown}, components: "
		s += ' '.join(map(str, self.components))
		return s

	def __str__(self):
		return f"Scaffold: {{{self.__repr__()}}}\n"
	

# ######### #
# functions #
# ######### #
def parseFasta(ifn):
	fasta = {}
	with open(ifn, 'r') as ifd:
		line = ifd.readline()
		while line != '':
			header = line.rstrip('\n')[1:]
			seq = ''
			line = ifd.readline()
			while line != '' and line[0] != '>':
				seq += line.rstrip('\n')
				line = ifd.readline()
			fasta[header] = seq
	return fasta

def generateScaffoldsFromAgpAndContigFasta(agp, fasta):
	scaffolds = []
	visited_scaffolds = set()

	heads_x_seqs = parseFasta(fasta)

	with open(agp, 'r') as ifd:
		line = ifd.readline() # grab the first line
		# skip the header lines
		while line != '' and line[0] == '#':
			line = ifd.readline()
		# process each data line
		while line != '':
			fields = line.rstrip('\n').split('\t')
			scaffold_name = fields[0]
			object_start = int(fields[1])
			object_end = int(fields[2])
			part_num = int(fields[3])
			component_type = fields[4]

			if not scaffold_name in visited_scaffolds:
				visited_scaffolds.add(scaffold_name)
				scaffolds.append(Scaffold(name=scaffold_name, components=[]))

			if part_num != (scaffolds[-1].num_components + 1):
				raise Agp2ScaffException(f"ERROR: AGP file likely malformed. Part number ({part_num}) has already been added to the scaffold from a previous record -OR- the part numbers are out-of-order / a part number has been skipped. Currently, {scaffolds[-1].num_components} parts have been added to the scaffold.")

			scaffolds[-1].addComponent(component_type, '\t'.join(fields[5:]), heads_x_seqs)
			line = ifd.readline()

	return scaffolds

def scaffoldsSortedIterate(scaffolds, sort_method, reverse):
	sorted_key = DefaultSortOrderKey.key_func

	if sort_method != "unordered":
		if sort_method == "alpha":
			sorted_key = lambda x: x.name
		elif sort_method == "length":
			sorted_key = lambda x: x.num_bases
		elif sort_method == "version":
			sorted_key = lambda x: LooseVersion(x.name)
		else:
			raise Agp2ScaffException("ERROR: sort method ({sort_method}) was not 'agp', 'alpha', 'length', 'unordered', or 'version'.")

	for scaffold in sorted(scaffolds, key=sorted_key, reverse=reverse):
		yield scaffold

def replaceRnaScaffsWithHicComponents(hics, rnas):
	out = []
	for rna in rnas:
		out.append(rna.supplantFrom(hics))
	return out


# #### #
# main #
# #### #
if __name__ == "__main__":

	# parse command line arguments
	contigs_fa_fn = "contigs.fa"
	hic_fa_fn = "hic.fa"
	hic_agp_fn = "hic.agp"
	rna_agp_fn = "rna.agp"
	out_agp_fn = "merged.agp"
	out_fa_fn = "merged.fa"
	out2_agp_fn = "merged_fixedGaps.agp"
	out2_fa_fn = "merged_fixedGaps.fa"
	sort_method = "version"
	reverse = False

	contigs_fa = pl.Path(contigs_fa_fn)
	hic_fa = pl.Path(hic_fa_fn)
	hic_agp = pl.Path(hic_agp_fn)
	rna_agp = pl.Path(rna_agp_fn)
	out_agp = pl.Path(out_agp_fn)
	out_fa = pl.Path(out_fa_fn)
	out2_agp = pl.Path(out2_agp_fn)
	out2_fa = pl.Path(out2_fa_fn)

	# read the input files and create the scaffolds in memory
	hic_scaffs = generateScaffoldsFromAgpAndContigFasta(hic_agp, contigs_fa)
	rna_scaffs = generateScaffoldsFromAgpAndContigFasta(rna_agp, hic_fa)

	scaffs = replaceRnaScaffsWithHicComponents(hic_scaffs, rna_scaffs)

	# write the output
	with open(out_agp, 'w') as ofd:
		with open(out_fa, 'w') as ofd_fa:
			with open(out2_agp, 'w') as ofd2:
				with open(out2_fa, 'w') as ofd2_fa:
					for scaffold in scaffoldsSortedIterate(scaffs, sort_method, reverse):
						ofd.write(scaffold.getAgp())
						ofd_fa.write(scaffold.getFasta())
						scaffold.setFixedUnknownGapLength(100)
						ofd2.write(scaffold.getAgp())
						ofd2_fa.write(scaffold.getFasta())

