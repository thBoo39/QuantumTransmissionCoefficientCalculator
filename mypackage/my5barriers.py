# list of effective mass from www.ioffe.ru/SVA/NSM/Semicond/index.html
m = {'GaAs': 0.063,
     'AlxGa1-xAs': lambda x: 0.063+0.083*x,
     'AlAs': 0.146,
     'InN': 0.11,
     'GaN': 0.20,
     'InxGa1-xN': lambda x: 0.2+x*(0.11-0.2),
     'InAs': 0.023}

# Set Up structure Here ---------------------------------------
# unit:
# potential (eV)
# position (nm)
# effective_mass (m_0 electron free mass)
# Multiple Barriers
potential = [0.0, 0.4, 0.0, 0.4, 0.0, 0.4, 0.0, 0.4, 0.0, 0.4, 0.0]
thickness = [4, 2, 6, 2, 6, 2, 6, 2, 6, 2, 4]
effective_mass = [m['GaAs'], m['GaAs'], m['GaAs'], m['GaAs'],
                  m['GaAs'], m['GaAs'], m['GaAs'], m['GaAs'], m['GaAs'],
                  m['GaAs'], m['GaAs']]
