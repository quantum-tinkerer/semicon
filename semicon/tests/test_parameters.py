import collections
import pytest

import pandas as pd

from semicon.parameters import DataBank


@pytest.mark.parametrize('databank_name', ['winkler', 'lawaetz'])
def test_databank_loading(databank_name):
    db = DataBank(databank_name)
    assert isinstance(db, collections.abc.Mapping)


@pytest.mark.parametrize('databank_name', ['winkler', 'lawaetz'])
def test_databank_dataframe_conversion(databank_name):
    db = DataBank(databank_name)
    df = db.to_dataframe()

    assert sorted(list(db)) == sorted(list(df.index))
    assert isinstance(df, pd.DataFrame)
