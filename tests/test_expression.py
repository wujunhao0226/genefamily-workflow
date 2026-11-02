import numpy as np

from genefamily.expression import summarise, tpm, zscore


def test_tpm_normalises_each_column_to_one_million():
    counts = np.array([[200.0, 100.0], [200.0, 100.0]])
    lengths = np.array([1000.0, 2000.0])
    matrix = tpm(counts, lengths)
    assert np.allclose(matrix.sum(axis=0), 1e6)
    # gene 1 is half as long, so it has twice the RPK of gene 2
    assert np.allclose(matrix[0] / matrix[1], 2.0)


def test_zscore_row_mean_zero_and_sd_one():
    matrix = np.array([[1.0, 2.0, 3.0], [5.0, 5.0, 5.0]])
    z = zscore(matrix)
    assert np.allclose(z[0].mean(), 0.0)
    assert np.allclose(z[0].std(ddof=1), 1.0)
    assert np.allclose(z[1], 0.0)          # constant row guarded


def test_summarise_reports_best_sample():
    tpm_matrix = np.array([[1.0, 9.0], [5.0, 5.0]])
    rows = summarise(["G1", "G2"], ["s1", "s2"], tpm_matrix)
    assert rows[0][3] == "s2"
    assert rows[0][4] == 2
