import pandas as pd
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


_banks_names = ['winkler', 'lawaetz']
def load_params(bankname):
    """Load material parameters from specified databank.

    output: pandas dataframe
    """
    if bankname not in _banks_names:
        msg = "Unkown bankname. Possible options are {}"
        raise ValueError(msg.format(_banks_names))
    fname = 'bank_' + bankname + '.csv'
    fpath = os.path.join(BASE_DIR, 'databank', fname)
    return pd.read_csv(fpath, index_col=0)
