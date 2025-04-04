
# parse sdevice cmd and par files and extract parameter values.
# Note this parser is _not_ very robust... fails in some cases that the Sentaurus parser would allow.
# requires brackets on separate rows, section headings on separate rows, spaces for delimiters (not tabs) etc.

class sdevice_file:
	def __init__(self, file):
		# print 'Reading sdevice file', file
		lines = []
		with open(file, 'r') as f:
			for line in f:
				lines.append(line.strip().split())
		
		# create a dictionary to store various sections of data.
		par = dict()
		material_files = dict()
		region_files = dict()
		r = iter(lines)
		for row in r:
			if len(row)>0:
				#print('is true?',row[0]==b'Material')
				if row[0] == 'Bandgap':
					# Found Bandgap section
					par['Bandgap'] = dict()
					row = r.__next__()
				
					while row[0] != '}':
						if row[0] == 'Eg0':
							par['Bandgap']['Eg0'] = float(row[2])
						if row[0] == 'alpha':
							par['Bandgap']['alpha'] = float(row[2])
						if row[0] == 'beta':
							par['Bandgap']['beta'] = float(row[2])
						if row[0] == 'Tpar':
							par['Bandgap']['Tpar'] = float(row[2])
						row = r.__next__()
					# print 'Bandgap', par['Bandgap']

				elif row[0] == 'ComplexRefractiveIndex': #TableODB no longer supported changing to Complex Refractive Index table format pwils
					# Found TableODB section
					par['TableODB'] = dict()
					par['TableODB']['wl'] = []
					par['TableODB']['n'] = []
					par['TableODB']['k'] = []
					row = r.__next__()
					while row != [')']: # Changed } to ) as complex refractive index is now in a numeric table pwils
						if len(row) == 3:
							if row[0] != '*' and row[0] != '{' and row[0] != 'Formula' and row[0] != 'NumericalTable(' and row[0] != 'TableInterpolation':
								par['TableODB']['wl'].append(float(row[0]))
								par['TableODB']['n'].append(float(row[1]))
								par['TableODB']['k'].append(float(row[2].rstrip(';')))
						row = r.__next__()
				elif row[0] == 'RadiativeRecombination':
					# Found TableODB section
					par['RadiativeRecombination'] = dict()
					row = r.__next__()
					while row[0] != '}':
						if row[0] == 'C':
							par['RadiativeRecombination']['C'] = float(row[2])
						row = r.__next__()
				elif row [0] == 'Material':
					material_files[row[2].strip('"')] = row[6].strip('"')
				elif row [0] == 'Region':
					region_files[row[2].strip('"')] = row[6].strip('"')

		self.par = par
		self.material_files = material_files
		self.region_files = region_files


	# Calculate Band gap for a given temperature.
	def Eg_T(self, T):
		assert 'Bandgap' in self.par , "No bandgap data in layer {0}, material {1}".format(self.name, self.material)
		bg = self.par['Bandgap']
		assert 'Eg0' in bg
		assert 'alpha' in bg
		assert 'beta' in bg
		assert 'Tpar' in bg

		Eg = bg['Eg0'] + bg['alpha']*(bg['Tpar']**2)/(bg['Tpar']+bg['beta']) - bg['alpha']*(T**2)/(T+bg['beta'])
		return Eg
		
		
