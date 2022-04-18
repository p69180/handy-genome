import importlib
top_package_name = __name__.split('.')[0]

common = importlib.import_module('.'.join([top_package_name, 'common']))
readhandler = importlib.import_module('.'.join([top_package_name, 'readplus', 'readhandler']))

'''
class Params(common.Params):
	pass
'''


class ReadPlusPair:
	'''
	Requires two readplus objects which are mates of each other and all of primary alignment.
	Non-primary reads with the same queryname is stored in 'rplist_nonprimary' attribute.
	'''
	def __init__(self, rplist_primary, rplist_nonprimary, chromdict):
		assert len(rplist_primary) == 2

		set_rp1_rp2(self, rplist_primary, chromdict)
		set_basic_attributes(self, rplist_nonprimary)
		set_is_proper_pair(self)
		set_SV_supporting(self)
		set_irrelevant(self)


def set_rp1_rp2(rpp, rplist_primary, chromdict):
	order = common.get_order(
			rplist_primary[0].read.reference_name, 
			rplist_primary[0].fiveprime_end,
			rplist_primary[1].read.reference_name, 
			rplist_primary[1].fiveprime_end,
			chromdict,
			)
	if order <= 0:
		rpp.rp1 = rplist_primary[0]
		rpp.rp2 = rplist_primary[1]
	else:
		rpp.rp1 = rplist_primary[1]
		rpp.rp2 = rplist_primary[0]


def set_basic_attributes(rpp, rplist_nonprimary):
	rpp.rplist_nonprimary = rplist_nonprimary
	rpp.query_name = rpp.rp1.read.query_name
	rpp.is_TRA = rpp.rp1.read.reference_name != rpp.rp2.read.reference_name
	rpp.pairorient = readhandler.get_pairorient(rpp.rp1.read) # may be None (mate unmapped, TLEN == 0)
	rpp.template_length = None if rpp.is_TRA else abs(rpp.rp1.read.template_length)


def set_is_proper_pair(rpp):
	if rpp.pairorient == None:
		rpp.is_proper_pair = False
	else:
		rpp.is_proper_pair = ( rpp.pairorient[0] == 'F' and rpp.pairorient[2] == 'R' )


def set_SV_supporting(rpp, threshold_template_length = common.THRESHOLD_TEMPLATE_LENGTH):
	rpp.SV_supporting = \
			(not rpp.is_proper_pair) or \
			rpp.is_TRA or \
			rpp.template_length > threshold_template_length


def set_irrelevant(rpp):
	rpp.irrelevant = rpp.rp1.irrelevant or rpp.rp2.irrelevant




