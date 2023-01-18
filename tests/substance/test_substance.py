import numpy as np

import reactord as rd


@pytest.mark.parametrize("name", compounds)
def test_create_substance_file(name):
    substance = rd.Substance(name)
    rd.Substance.to_pickle(substance, "name_file")
    file = rd.Substance.from_pickle("name_file")
    assert file.name == substance.name
