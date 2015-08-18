"""
Classes for parsing MEME results out of XML files.
"""

__version__ = '0.1'


### IMPORTS

from xml.etree import ElementTree
from collections import namedtuple


### CONSTANTS & DEFINES

### CODE ###

def find (needle, haystack):
	"""
	Return the index of the needle in the haystack
	
	Parameters:
		needle: any iterable
		haystack: any other iterable
		
	Returns:
		the index of the start of needle or -1 if it is not found.
	
	Looking for a sub-list of a list is actually a tricky thing. This
	approach uses the Boyer-Moore-Horspool algorithm. Needle and haystack
	should be any iterable, as long as their elements are hashable.

	Example:
	
		>>> find ([1, 2], [1, 1, 2])
		1
		>>> find ((1, 2, 3), range (10))
		1
		>>> find ('gh', 'abcdefghi')
		6
		>>> find ([2, 3], [7, 8, 9])
		-1

	"""
	h = len (haystack)
	n = len (needle)
	skip = {needle[i]: n - i - 1 for i in range(n - 1)}
	i = n - 1
	while i < h:
		for j in range(n):
			if haystack[i - j] != needle[-j - 1]:
				i += skip.get(haystack[i], n)
				break
		else:
			return i - n + 1
	return -1

def iter_find (needle, haystack):
	"""
	Iterates over the index of the needle in the haystack
	
	Parameters:
		needle: as per `find`
		haystack: as per `find`
		
	Returns:
		A generator function over the indices of needle in haystack
	
	This behaves as per `find` except it iterates over occurences of
	needle and so may be used to find more than 1 instance, unlike
	the original algorithm.

	Example:
	
		>>> find ([1, 2], [1, 1, 2])
		1
		>>> find ((1, 2, 3), range (10))
		1
		>>> find ('gh', 'abcdefghi')
		6
		>>> find ([2, 3], [7, 8, 9])
		-1

	"""
	h = len (haystack)
	n = len (needle)
	skip = {needle[i]: n - i - 1 for i in range(n - 1)}
	i = n - 1
	while i < h:
		for j in range(n):
			if haystack[i - j] != needle[-j - 1]:
				i += skip.get(haystack[i], n)
				break
		else:
			yield i - n + 1
			i += 1
	

def motif_to_item (m, rev=False):
	"""
	Given motif information, produce an item (symbol)
	
	Parameters:
		m: a motif information tuple
		rev: if True, produce the symbol for the opposite strand, i.e. flip
			the sign
	"""
	# XXX: unsure if motifs are always nicely named but will abreviate for
	# time being
	strand = m.strand
	if rev:
		if strand == 'plus': strand = 'minus'
		elif strand == 'minus': strand = 'plus'
		
	item_name = '%s_%s' % (m.motif_id, strand)
	item_name = item_name.replace ('motif_', '')
	item_name = item_name.replace ('_plus', '+').replace ('_minus', '-')
	
	return item_name



def parse_meme_results (in_pth_or_hndl):
	## Subfunctions:
	def tree_to_obj():
		# first get overview info
		ts_node = tree.findall ('./training_set')[0]
		results_obj.datafile = ts_node.attrib['datafile']
		
		cl_node = tree.findall ('./model/command_line')
		results_obj.command_line = cl_node[0].text		
		
		# get sequence information
		for sn in tree.findall ('./training_set/sequence'):
			results_obj.seqs.append (SeqInfo (
				sn.attrib['id'],
				sn.attrib['name'],
				sn.attrib['length'],
				sn.attrib['weight'],
			)
		)
		
		# get summary - all motifs arranged over all seqs
		for seq in tree.findall ('./scanned_sites_summary/scanned_sites'):
			children = []
			for child in seq.iter ('scanned_site'):
				children.append (SeqSite (
						child.attrib['motif_id'],
						child.attrib['strand'],
						int (child.attrib['position']),
						float (child.attrib['pvalue']),
					)
				)
			results_obj.scanned_sites_summary.append (SiteSummary (
					seq.attrib['sequence_id'],
					float (seq.attrib['pvalue']),
					children,
				)
			)
	
		# extract the motif information
		motifs = tree.findall ('./motifs')[0]
		for child_motif in motifs:
			motif = Motif (
				child_motif.attrib['id'],
				child_motif.attrib['name'],
				int (child_motif.attrib['width']),
				int (child_motif.attrib['sites']),
				float (child_motif.attrib['ic']),
				float (child_motif.attrib['re']),
				int (child_motif.attrib['llr']),
				float (child_motif.attrib['e_value']),
				float (child_motif.attrib['bayes_threshold']),
				float (child_motif.attrib['elapsed_time']),
			)
			results_obj.motifs.append (motif)
	
	## Preconditions & preparations:
	if hasattr (in_pth_or_hndl, 'open'):
		# it's a raw handle
		in_hndl = in_pth_or_hndl
	else:
		# it's a raw string
		in_hndl = open (in_pth_or_hndl, 'rU')
		hndl_opened = True
	
	tree = ElementTree.parse (in_hndl)
	results_obj = MemeResults()
	tree_to_obj ()
	
	## Postconditions & return:
	return results_obj
	
	
Motif = namedtuple ('Motif',
	['id', 'name', 'width', 'sites', 'ic', 're',
	'llr', 'e_value', 'bayes_threshold', 'elapsed_time'])

SiteSummary = namedtuple ('SiteSummary',
	['sequence_id', 'pvalue', 'scanned_sites'])

SeqSite = namedtuple ('SeqSite',
	['motif_id', 'strand', 'position', 'pvalue'])
	
SeqInfo = namedtuple ('SeqInfo',
	['id', 'name', 'length', 'weight'])


class MemeResults (object):
	def __init__ (self):
		self.datafile = None
		self.seqs = []
		self.scanned_sites_summary = []
		self.motifs = []
		
	def seq_id_to_name (self, seq_id):
		for s in self.seqs:
			if s.id == seq_id:
				return s.name
		raise exceptions.ValueError ("no sequence with id '%s'" % seq_id)
		
	def seq_name_to_id (self, seq_name):
		for s in self.seqs:
			if s.name == seq_name:
				return s.id
		raise exceptions.ValueError ("no sequence with name '%s'" % seq_name)
		
	def get_motif_by_id (self, motif_id):
		for m in self.motifs:
			if motif_id == m.id:
				return m
		raise exceptions.ValueError ("no motif with id '%s'" % motif_id)
		
	def __repr__ (self):
		"""
		Return a nicely formatted representation string
		"""
		return '%(classname)s(datafile=%(datafile)s, seqs=%(seqs)s, ' \
			'scanned_sites_summary=%(summary)s), motifs=%(motifs)s' %  {
				'classname': self.__class__.__name__,
				'datafile': self.datafile,
				'seqs': self.seqs,
				'summary': self.scanned_sites_summary,
				'motifs': self.motifs,
		}
	
	
### TESTING

if __name__ == '__main__':
	print (parse_meme_results ('meme.xml'))


### END ###
