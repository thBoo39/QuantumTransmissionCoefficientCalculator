from matplotlib import rc, rcParams
import sys
import os


def setup_fonts():
    # Fonts setup
    # activate latex text rendering
    # rc('text', usetex=True)
    rc('axes', linewidth=2)
    rc('axes', labelweight='bold')
    rc('axes', labelsize='large')
    rc('legend', fontsize='medium')
    rc('lines', linewidth=3)
    rc('font', weight='bold')
    rc('mathtext', default='regular')
    # rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
    rcParams.update({'font.size': 16})
    return


def redirect_output_to_file(target_file='test_hcsc_log.txt'):
    fname = os.path.join(os.getcwd(), target_file)
    sys.stdout = open(fname, 'w')
    return
