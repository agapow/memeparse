"""
Tests for memeparse.
"""
	
DATA_FILE = 'tests/meme.xml'


class test_memeparse (object):
	# overhead
	def setup (self):
		print ("settingup test_memeparse class")
		import memeparse
		self.meme_results = memeparse.parse_meme_results (DATA_FILE)
		
	def teardown (self):
		print ("tearing down test_memeparse class")
	
	# tests
	def test_num_motifs (self):
		assert len (self.meme_results.seqs) == 100, "wrong number of sequences"
		# just test the very first scanned site
		assert len (self.meme_results.scanned_sites_summary[0].scanned_sites) == 11, \
			"wrong number of scanned sites"
		assert len (self.meme_results.motifs) == 8, "wrong number of motifs"


