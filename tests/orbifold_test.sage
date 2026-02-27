"""Unit tests for orbifold.sage.

Tests lambda-sector basis computation for zero, trivial, and
non-trivial sectors.

Run with: sage orbifold_test.sage
"""

load('orbifold.sage')


def assert_same_elements(left_list, right_list):
    """Asserts two lists contain the same elements, ignoring order.

    Converts each element to a tuple for hashing, then checks that
    the sets are equal and have the same cardinality.

    Args:
        left_list: First list of sequences to compare.
        right_list: Second list of sequences to compare.
    """
    left_set = set([tuple(element) for element in left_list])
    right_set = set([tuple(element) for element in right_list])
    assert left_set.intersection(right_set) == left_set
    assert len(left_set) == len(right_set)


def test_get_lambda_sector_basis_zero_sector():
    """Tests lambda-sector basis for the zero sector lambda = (0, 0).

    With a 2x2 chain matrix, the zero sector should yield two basis
    elements with all coordinates fixed.
    """
    A = matrix([[2, 1], [0, 3]])
    y_lambda = (0, 0)

    basis = get_lambda_sector_basis(A, {}, y_lambda, 3)
    assert_same_elements(
        [b.to_tuple() for b in basis],
        [((2, 1), (0, 0), (0, 1)), ((1, 2), (0, 0), (0, 1))]
    )

def test_get_lambda_sector_basis_trivial_sector():
    """Tests lambda-sector basis for the trivial sector lambda = (1, 1).

    With no fixed coordinates, the sector contributes a single element
    with zero Koszul part.
    """
    A = matrix([[2, 1], [0, 3]])
    y_lambda = (1, 1)

    basis = get_lambda_sector_basis(A, {}, y_lambda, 3)
    assert_same_elements(
        [b.to_tuple() for b in basis],
        [((0, 0), (1, 1), ())]
    )

def test_get_lambda_sector_basis_with_sector():
    """Tests lambda-sector basis for a non-trivial sector.

    With a 3x3 Fermat matrix and lambda = (0, 0, 1), two of three
    coordinates are fixed, yielding two basis elements.
    """
    A = matrix([[3, 0, 0], [0, 3, 0], [0, 0, 3]])
    y_lambda = (0, 0, 1)

    basis = get_lambda_sector_basis(A, {}, y_lambda, 3)
    assert_same_elements(
        [b.to_tuple() for b in basis],
        [((1, 2, 0), (0, 0, 1), (0, 1)), ((2, 1, 0), (0, 0, 1), (0, 1))]
    )


if __name__ == "__main__":
    print("Running orbifold.sage tests...")
    test_get_lambda_sector_basis_zero_sector()
    test_get_lambda_sector_basis_trivial_sector()
    test_get_lambda_sector_basis_with_sector()
    print("All passed.")
