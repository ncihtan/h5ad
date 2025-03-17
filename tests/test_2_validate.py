import pytest
from validate import validate


def test_validate():
    """Test main validate function with incorrect path and
    non-h5ad format file """

    incorrect_path = "tests/somewhere.h5ad"
    non_h5ad = "tests/not_real_h5ad"

    with pytest.raises(SystemExit) as excinfo:
        validate(incorrect_path, "test_validate_out.txt")
    assert excinfo.value.code == 2

    with pytest.raises(SystemExit) as excinfo:
        validate(non_h5ad, "test_validate_out.txt")
    assert excinfo.value.code == 2
