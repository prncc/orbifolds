"""Unit tests for potentials.sage.

Tests matrix decomposition into basic potentials, matrix builders,
basis computation, and variable permutation invariance.

Run with: sage potentials_test.sage
"""

load("potentials.sage")


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


def test_decompose_into_basic_potentials():
    """Tests decomposition across 11 variable permutations.

    Constructs a matrix with 2 chains, 2 loops, and 1 Fermat block,
    then verifies that decompose_into_basic_potentials correctly
    identifies all 5 components regardless of variable ordering.
    """
    chain = build_chain_matrix((2, 2))
    loop = build_loop_matrix((3, 3, 3))
    fermat = build_fermat_matrix((4, 4, 4, 4))
    A = block_diagonal_matrix([chain, chain, loop, loop, fermat])
    for variable_ordering in [
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
        [9, 4, 8, 10, 1, 3, 11, 13, 0, 7, 12, 5, 6, 2],
        [10, 13, 4, 6, 7, 9, 0, 2, 11, 5, 3, 12, 8, 1],
        [10, 2, 8, 11, 3, 6, 1, 5, 12, 7, 0, 9, 4, 13],
        [6, 11, 5, 10, 7, 1, 4, 13, 8, 2, 9, 3, 12, 0],
        [10, 2, 8, 3, 6, 12, 9, 0, 13, 4, 11, 7, 5, 1],
        [6, 1, 2, 9, 3, 7, 11, 13, 8, 4, 5, 12, 10, 0],
        [4, 5, 7, 6, 11, 10, 12, 8, 2, 0, 1, 13, 3, 9],
        [11, 13, 2, 10, 12, 4, 7, 6, 3, 0, 1, 5, 9, 8],
        [1, 8, 6, 13, 9, 2, 7, 0, 5, 10, 11, 4, 12, 3],
        [5, 9, 4, 8, 3, 2, 10, 12, 11, 1, 13, 7, 0, 6],
    ]:
        A = A.matrix_from_rows_and_columns(
            variable_ordering,
            variable_ordering,
        )
        potentials = decompose_into_basic_potentials(A)
        assert len(potentials) == 5
        assert len([p for p in potentials if p.kind == "fermat"]) == 1
        assert len([p for p in potentials if p.kind == "chain"]) == 2
        assert len([p for p in potentials if p.kind == "loop"]) == 2


def test_build_chain_matrix():
    """Tests that build_chain_matrix produces the expected structure."""
    chain = build_chain_matrix((2, 3))
    assert chain == matrix([[2, 1], [0, 3]])


def test_build_fermat_matrix():
    """Tests that build_fermat_matrix produces a diagonal matrix."""
    chain = build_fermat_matrix((2, 3))
    assert chain == matrix([[2, 0], [0, 3]])


def test_build_loop_matrix():
    """Tests that build_loop_matrix produces the expected cyclic structure."""
    chain = build_loop_matrix((2, 3))
    assert chain == matrix([[2, 1], [1, 3]])


def _check_potential_basis(potential, expected_basis, expected_expanded_basis):
    """Validates a basic potential's basis and its expanded form.

    Args:
        potential: A BasicPotential object to check.
        expected_basis: Expected basis elements in sub-potential coordinates.
        expected_expanded_basis: Expected basis elements expanded to
            full matrix coordinates.
    """
    basis = potential.get_basis()
    expanded_basis = potential.get_basis(expand_to_full_matrix=True)
    assert len(basis) == get_milnor_number(potential.matrix)
    assert_same_elements(basis, expected_basis)
    assert_same_elements(expanded_basis, expected_expanded_basis)


def test_basic_potential_bases():
    """Tests basis computation for each basic potential type.

    Constructs a matrix with one chain, one loop, and one Fermat block,
    decomposes it, and checks each sub-potential's basis against
    known expected values.
    """
    chain = build_chain_matrix([2, 2, 3])
    loop = build_loop_matrix([2, 2, 2])
    fermat = build_fermat_matrix([3, 3])
    A = block_diagonal_matrix([chain, loop, fermat])
    potentials = decompose_into_basic_potentials(A)

    chain_potential = [p for p in potentials if p.kind == "chain"][0]
    _check_potential_basis(
        chain_potential,
        [
            [1, 1, 1],
            [1, 1, 2],
            [1, 1, 3],
            [1, 2, 1],
            [1, 2, 2],
            [1, 2, 3],
            [2, 1, 1],
            [2, 1, 2],
        ],
        [
            [1, 1, 1, 0, 0, 0, 0, 0],
            [1, 1, 2, 0, 0, 0, 0, 0],
            [1, 1, 3, 0, 0, 0, 0, 0],
            [1, 2, 1, 0, 0, 0, 0, 0],
            [1, 2, 2, 0, 0, 0, 0, 0],
            [1, 2, 3, 0, 0, 0, 0, 0],
            [2, 1, 1, 0, 0, 0, 0, 0],
            [2, 1, 2, 0, 0, 0, 0, 0],
        ],
    )

    loop_potential = [p for p in potentials if p.kind == "loop"][0]
    _check_potential_basis(
        loop_potential,
        [
            [1, 1, 1],
            [1, 1, 2],
            [1, 2, 1],
            [1, 2, 2],
            [2, 1, 1],
            [2, 1, 2],
            [2, 2, 1],
            [2, 2, 2],
        ],
        [
            [0, 0, 0, 1, 1, 1, 0, 0],
            [0, 0, 0, 2, 1, 1, 0, 0],
            [0, 0, 0, 1, 1, 2, 0, 0],
            [0, 0, 0, 2, 1, 2, 0, 0],
            [0, 0, 0, 1, 2, 1, 0, 0],
            [0, 0, 0, 2, 2, 1, 0, 0],
            [0, 0, 0, 1, 2, 2, 0, 0],
            [0, 0, 0, 2, 2, 2, 0, 0],
        ],
    )

    fermat_potential = [p for p in potentials if p.kind == "fermat"][0]
    _check_potential_basis(
        fermat_potential,
        [[1, 1], [1, 2], [2, 1], [2, 2]],
        [
            [0, 0, 0, 0, 0, 0, 1, 1],
            [0, 0, 0, 0, 0, 0, 1, 2],
            [0, 0, 0, 0, 0, 0, 2, 1],
            [0, 0, 0, 0, 0, 0, 2, 2],
        ],
    )


def test_potential_get_basis():
    """Tests get_gamma_basis with permuted variables.

    Constructs a matrix from chain and Fermat blocks, permutes the
    variables, and verifies the combined basis matches expectations
    and has the correct Milnor number.
    """
    A = block_diagonal_matrix(
        [
            build_chain_matrix([2, 3]),
            build_fermat_matrix([3, 2]),
        ]
    ).matrix_from_rows_and_columns([0, 2, 3, 1], [0, 2, 3, 1])
    basis = get_gamma_basis(A)
    expected_unshuffled_basis = [
        [1, 1, 1, 1],
        [1, 1, 2, 1],
        [1, 2, 1, 1],
        [1, 2, 2, 1],
        [1, 3, 1, 1],
        [1, 3, 2, 1],
        [2, 1, 1, 1],
        [2, 1, 2, 1],
    ]
    assert len(basis) == get_milnor_number(A)
    assert_same_elements(
        basis,
        [[b[0], b[2], b[3], b[1]] for b in expected_unshuffled_basis],
    )


if __name__ == "__main__":
    print("Running potentials.sage tests...")
    test_build_chain_matrix()
    test_build_fermat_matrix()
    test_build_loop_matrix()
    test_decompose_into_basic_potentials()
    test_basic_potential_bases()
    test_potential_get_basis()
    print("All passed.")
