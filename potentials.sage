"""Matrix potential decomposition and Koszul cohomology basis computation.

Decomposes invertible integer matrices into Fermat, chain, and loop basic
potentials, and computes cohomology bases for each type.
"""

import dataclasses
import enum
import itertools


def get_milnor_number(A):
    """Computes the size of the Koszul cohomology basis of A.

    Uses the formula prod(1/q_i - 1) where q = A^{-1} J.

    Args:
        A: Square integer matrix defining the potential.

    Returns:
        The Milnor number (basis size) as an integer.
    """
    epsilon = vector([1] * A.nrows())
    total = 1
    for qq in A.inverse() * epsilon:
        total *= 1 / (qq) - 1
    return total


def get_chain_basis(A, log_differentials=True):
    """Returns the Koszul cohomology basis for a chain potential.

    Recursively constructs basis elements for the chain potential defined
    by matrix A, using the structure of the chain's off-diagonal entries.

    Args:
        A: Square integer matrix defining the chain potential.
        log_differentials: If True, shift all basis indices by +1 to
            convert from polynomial to log-differential coordinates.

    Returns:
        List of integer lists, each representing a basis element as
        exponent vectors.
    """
    N = A.nrows()
    if N == 1:
        basis = []
        for i in range(0, A[0, 0] - 1):
            basis.append([i])
    else:
        basis = [[]]
        samples = []
        samples.append(range(A[0, 0] - 1))
        for i in range(1, N):
            samples.append(range(A[i, i]))
        for sample in samples:
            basis = [b + [s] for b in basis for s in sample]
        if N == 2:
            basis.append([A[0, 0] - 1, 0])
        else:
            sub_A = A.matrix_from_rows_and_columns(range(2, N), range(2, N))
            for b in get_chain_basis(sub_A, log_differentials=False):
                basis.append([A[0, 0] - 1, 0] + b)
    if log_differentials:
        return [
            [basis_element[i] + 1 for i in range(len(basis_element))]
            for basis_element in basis
        ]
    return basis


def get_loop_basis(A, log_differentials=True):
    r"""Returns the Koszul cohomology basis for a loop potential.

    Constructs all basis elements by taking the Cartesian product of
    ranges [0, a_ii) for each diagonal entry a_ii of A.

    Args:
        A: Square integer matrix defining the loop potential.
        log_differentials: If True, shift all basis indices by +1 to
            use log-differential coordinates.

    Returns:
        List of integer lists, each representing a basis element as
        exponent vectors.
    """
    N = A.nrows()
    if N == 1:
        basis = []
        for i in range(0, A[0, 0] - 1):
            basis.append([i])
    else:
        basis = [[]]
        samples = []
        for i in range(N):
            samples.append(range(A[i, i]))
        for sample in samples:
            basis = [b + [s] for b in basis for s in sample]
    if log_differentials:
        return [
            [basis_element[i] + 1 for i in range(len(basis_element))]
            for basis_element in basis
        ]
    return basis


def get_fermat_basis(A, log_differentials=True):
    r"""Returns the Koszul cohomology basis for a Fermat potential.

    Constructs all basis elements by taking the Cartesian product of
    ranges [0, a_ii - 1) for each diagonal entry a_ii of A.

    Args:
        A: Square integer matrix defining the Fermat potential (diagonal).
        log_differentials: If True, shift all basis indices by +1 to
            use log-differential coordinates.

    Returns:
        List of integer lists, each representing a basis element as
        exponent vectors.
    """
    N = A.nrows()
    if N == 1:
        basis = []
        for i in range(A[0, 0] - 1):
            basis.append([i])
    else:
        basis = [[]]
        samples = []
        for i in range(N):
            samples.append(range(A[i, i] - 1))
        for sample in samples:
            basis = [b + [s] for b in basis for s in sample]
    if log_differentials:
        return [
            [basis_element[i] + 1 for i in range(len(basis_element))]
            for basis_element in basis
        ]
    return basis


class BasicPotentialType(enum.StrEnum):
    """Enumeration of the three basic potential types."""

    fermat = "fermat"
    chain = "chain"
    loop = "loop"


@dataclasses.dataclass
class BasicPotential:
    """A basic (Fermat, chain, or loop) sub-potential of a matrix A.

    Represents one indecomposable component after decomposing A into
    basic potential blocks.

    Attributes:
        kind: The type of basic potential (fermat, chain, or loop).
        indices_in_matrix: Tuple of row/column indices in the full matrix
            that belong to this sub-potential.
        full_matrix: The original full matrix A before decomposition.
    """

    kind: BasicPotentialType
    indices_in_matrix: tuple
    full_matrix: object

    @property
    def matrix(self):
        """Returns the sub-matrix of A corresponding to this potential.

        Returns:
            The square sub-matrix extracted from full_matrix at the
            rows and columns given by indices_in_matrix.
        """
        return self.full_matrix.matrix_from_rows_and_columns(
            self.indices_in_matrix,
            self.indices_in_matrix,
        )

    def get_basis(self, expand_to_full_matrix=False):
        """Computes the Koszul cohomology basis for this basic potential.

        Args:
            expand_to_full_matrix: If True, embed each basis element into
                the full matrix dimension by zero-padding.

        Returns:
            List of integer lists representing basis elements. If
            expand_to_full_matrix is True, each list has length
            full_matrix.nrows(); otherwise, length matches this
            sub-potential's dimension.
        """
        A = self.matrix
        if self.kind == "fermat":
            basis = get_fermat_basis(A)
        elif self.kind == "chain":
            basis = get_chain_basis(A)
        elif self.kind == "loop":
            basis = get_loop_basis(A)
        else:
            raise NotImplementedError(
                "Can only generate a basis for basic potential types."
                " Decompose your potential first."
            )
        if expand_to_full_matrix:
            expanded_basis = []
            for gamma in basis:
                expanded_gamma = [0] * self.full_matrix.nrows()
                for basis_index, expanded_basis_index in enumerate(
                    self.indices_in_matrix
                ):
                    expanded_gamma[expanded_basis_index] = gamma[basis_index]
                expanded_basis.append(expanded_gamma)
            return expanded_basis
        else:
            return basis


def is_fermat_monomial(A, row_index, variable_index):
    """Determines whether a monomial of A is a Fermat monomial.

    A monomial is Fermat if the only nonzero entry in its row (aside
    from the diagonal) is the diagonal entry itself.

    Args:
        A: Square integer matrix defining the potential.
        row_index: Row index of the monomial to check.
        variable_index: Column index of the diagonal entry (must be > 1).

    Returns:
        True if the monomial is Fermat, False otherwise.

    Raises:
        ValueError: If A[row_index][variable_index] <= 1.
    """
    row = A[row_index]
    if row[variable_index] <= 1:
        raise ValueError(
            "variable_index should point to the variable"
            " with the largest exponent in monomial."
        )
    sum_before = sum(row[:variable_index])
    sum_after = sum(row[variable_index + 1 :])
    return (sum_before == 0) and (sum_after == 0)


def is_connecting_variable(A, row_index, variable_index):
    """Determines whether a variable connects multiple monomials.

    Checks if the variable at variable_index appears in any monomial
    other than the one at row_index.

    Args:
        A: Square integer matrix defining the potential.
        row_index: Row index of the monomial containing this variable.
        variable_index: Column index of the variable to check (must be > 1
            at A[row_index]).

    Returns:
        True if the variable appears in other monomials, False otherwise.

    Raises:
        ValueError: If A[row_index][variable_index] <= 1.
    """
    if A[row_index][variable_index] <= 1:
        raise ValueError(
            "variable_index should point to the variable"
            " with the largest exponent in monomial."
        )
    column = A.columns()[variable_index]
    sum_before = sum(column[:row_index])
    sum_after = sum(column[row_index + 1 :])
    return (sum_before > 0) or (sum_after > 0)


def get_diagonal_value_index(row):
    """Finds the index of the (presumably unique) non-unit positive entry.

    Scans the row for the first entry greater than 1, which corresponds
    to the diagonal entry of a basic potential.

    Args:
        row: A row vector or sequence of integers.

    Returns:
        The index of the first entry > 1, or None if no such entry exists.
    """
    for j, value in enumerate(row):
        if value > 1:
            return j
    return None


def mark_chain_from_fermat(A, fermat_row_index, fermat_variable_index,
                           monomial_slots, mark):
    """Traces a chain's monomials in A starting from its Fermat term.

    Follows the off-diagonal connections from a Fermat monomial to
    identify all monomials belonging to the same chain sub-potential.

    Args:
        A: Square integer matrix defining the potential.
        fermat_row_index: Row index of the Fermat monomial at the
            start of the chain.
        fermat_variable_index: Column index of the Fermat variable.
        monomial_slots: Mutable list tracking which sub-potential each
            monomial belongs to (modified in place).
        mark: Integer label to assign to this chain's monomials.

    Returns:
        A BasicPotential of kind "chain" with the discovered indices.
    """
    monomial_slots[fermat_row_index] = mark
    indices_in_matrix = [fermat_row_index]

    connecting_column = A.columns()[fermat_variable_index]
    while 1 in connecting_column:
        # Find the other monomial that contains this variable.
        row_index = tuple(connecting_column).index(1)
        monomial_slots[row_index] = mark
        indices_in_matrix.append(row_index)

        # Find the next "diagonal" entry and to continue marking.
        column_index = get_diagonal_value_index(A[row_index])
        connecting_column = A.columns()[column_index]

    return BasicPotential(
        kind="chain",
        indices_in_matrix=tuple(reversed(indices_in_matrix)),
        full_matrix=A,
    )


def mark_loop(A, starting_row_index, starting_column_index,
              monomial_slots, mark):
    """Traces a loop's monomials in A starting from one of its terms.

    Follows the off-diagonal connections until returning to the starting
    monomial, identifying all monomials in the loop sub-potential.

    Args:
        A: Square integer matrix defining the potential.
        starting_row_index: Row index of the starting monomial.
        starting_column_index: Column index of the starting variable.
        monomial_slots: Mutable list tracking which sub-potential each
            monomial belongs to (modified in place).
        mark: Integer label to assign to this loop's monomials.

    Returns:
        A BasicPotential of kind "loop" with the discovered indices.
    """
    monomial_slots[starting_row_index] = mark
    indices_in_matrix = [starting_row_index]

    connecting_column = A.columns()[starting_column_index]
    next_row_index = tuple(connecting_column).index(1)

    # Iterate until we complete the "loop".
    while next_row_index != starting_row_index:
        monomial_slots[next_row_index] = mark
        connecting_column = A.columns()[next_row_index]
        indices_in_matrix.append(next_row_index)
        next_row_index = tuple(connecting_column).index(1)

    return BasicPotential(
        kind="loop",
        indices_in_matrix=tuple(reversed(indices_in_matrix)),
        full_matrix=A,
    )


def reconstruct_combined_potential(potentials):
    """Reconstructs matrix A from its basic potential decomposition.

    Assembles a block-diagonal matrix from the sub-potentials and then
    permutes rows and columns to recover the original variable ordering.

    Args:
        potentials: Sequence of BasicPotential objects comprising the
            decomposition.

    Returns:
        The reconstructed square integer matrix.
    """
    n = 0
    for one_potential in potentials:
        n += len(one_potential.indices_in_matrix)

    # We will add each sub-potential matrix into a block diagonal matrix.
    all_blocks = []

    # And then permute the variables to recover A in the original order.
    permutation_indices = [0] * n
    start_from_index = 0

    for potential in potentials:
        for from_index, to_index in enumerate(
            potential.indices_in_matrix, start=start_from_index
        ):
            permutation_indices[to_index] = from_index
        start_from_index = from_index + 1
        all_blocks.append(potential.matrix)

    A = block_diagonal_matrix(all_blocks)
    A = A.matrix_from_rows_and_columns(permutation_indices, permutation_indices)

    return A


def decompose_into_basic_potentials(A):
    """Decomposes a potential matrix into basic (Fermat, chain, loop) potentials.

    Identifies Fermat monomials first, distinguishes chains from isolated
    Fermats by checking for connecting variables, then assigns all
    remaining monomials to loops.

    Args:
        A: Square integer matrix defining the potential W_A.

    Returns:
        List of BasicPotential objects whose sub-matrices reconstruct A.

    Raises:
        ValueError: If A is not square or cannot be decomposed.
    """
    if A.nrows() != A.ncols():
        raise ValueError("A should be a square matrix.")

    n = A.nrows()
    # Each slot corresponds to a monomial of W_A. The value of
    # each slot indicates which sub-potential it belongs to.
    monomial_slots = [0] * n
    # 1 is reserved for all Fermat terms.
    fermat_mark = 1
    mark = 2

    basic_potentials = []
    fermat_indices = []

    # Only Fermats and chains contain Fermat monomials,
    # and are easy to distinguish.
    for row_index in range(n):
        if monomial_slots[row_index] != 0:
            continue
        row = A[row_index]
        variable_index = get_diagonal_value_index(row)
        if is_fermat_monomial(A, row_index, variable_index):
            if is_connecting_variable(A, row_index, variable_index):
                chain = mark_chain_from_fermat(
                    A, row_index, variable_index, monomial_slots, mark
                )
                basic_potentials.append(chain)
                mark += 1
            else:
                monomial_slots[row_index] = fermat_mark
                fermat_indices.append(row_index)

    if fermat_indices:
        basic_potentials.append(
            BasicPotential(
                kind="fermat",
                indices_in_matrix=tuple(fermat_indices),
                full_matrix=A,
            )
        )

    # All that remain are loops.
    for row_index in range(n):
        if monomial_slots[row_index] != 0:
            continue
        variable_index = get_diagonal_value_index(A[row_index])
        loop = mark_loop(A, row_index, variable_index, monomial_slots, mark)
        basic_potentials.append(loop)
        mark += 1

    # Check decomposition (should only raise an error if A was invalid).
    if A != reconstruct_combined_potential(basic_potentials):
        raise ValueError(
            "Could not decompose A into basic potentials. Check that"
            " A is a combination of loop, fermat, and chain blocks."
        )

    return basic_potentials


def get_gamma_basis(A):
    """Computes the full Koszul cohomology basis for an arbitrary potential.

    Decomposes A into basic potentials, computes each sub-basis, and
    takes their Cartesian product to form the combined basis with
    variables reordered to match A's original ordering.

    Args:
        A: Square integer matrix defining the potential W_A.

    Returns:
        List of integer lists, each representing a cohomology basis
        element as exponent vectors in the original variable ordering.
    """
    basic_potentials = decompose_into_basic_potentials(A)

    basic_bases = []
    permutation_indices = []
    for basic_potential in basic_potentials:
        basic_bases.append(basic_potential.get_basis())
        permutation_indices.extend(basic_potential.indices_in_matrix)

    basis = []
    for gamma_parts in itertools.product(*basic_bases):
        # Merge product of basic bases into a single element.
        unordered_gamma = []
        for part in gamma_parts:
            unordered_gamma.extend(part)

        # Order variables.
        gamma = [0] * len(unordered_gamma)
        for from_index, to_index in enumerate(permutation_indices):
            gamma[to_index] = unordered_gamma[from_index]

        basis.append(gamma)

    return basis


def build_chain_matrix(diagonal):
    """Builds a chain potential matrix from diagonal entries.

    Constructs a square matrix with the given diagonal and 1s on the
    first superdiagonal, representing a chain potential.

    Args:
        diagonal: Sequence of positive integers for the diagonal entries.

    Returns:
        A square integer matrix with the chain structure.
    """
    n = len(diagonal)
    rows = [[0] * n for _ in range(n)]
    for i, a_i in enumerate(diagonal):
        rows[i][i] = a_i
        if i < (n - 1):
            rows[i][i + 1] = 1
    return matrix(rows)


def build_loop_matrix(diagonal):
    """Builds a loop potential matrix from diagonal entries.

    Constructs a square matrix with the given diagonal and 1s on the
    first superdiagonal plus position (n-1, 0), forming a cyclic loop.

    Args:
        diagonal: Sequence of positive integers for the diagonal entries.

    Returns:
        A square integer matrix with the loop structure.
    """
    n = len(diagonal)
    rows = [[0] * n for _ in range(n)]
    for i, a_i in enumerate(diagonal):
        rows[i][i] = a_i
        rows[i][(i + 1) % n] = 1
    return matrix(rows)


def build_fermat_matrix(diagonal):
    """Builds a Fermat potential matrix from diagonal entries.

    Constructs a diagonal matrix, representing an uncoupled Fermat
    potential.

    Args:
        diagonal: Sequence of positive integers for the diagonal entries.

    Returns:
        A diagonal integer matrix.
    """
    return diagonal_matrix(diagonal)


def get_d(A):
    """Computes the minimal degree d of the potential W_A.

    Since A * weights = d * epsilon, we have weights = d * A^{-1} * epsilon,
    and d is the minimal positive integer that produces integer weights.
    Equivalently, d = lcm of denominators of A^{-1} * epsilon.

    This d also clears denominators of epsilon * A^{-T}, making it useful
    for writing symmetry group elements as xi_d^(d * lambda * A^{-T}).

    Args:
        A: Square invertible integer matrix defining the potential.

    Returns:
        The minimal degree as a positive integer.
    """
    epsilon = vector([1] * A.nrows())
    weights = A.inverse() * epsilon
    return lcm([denominator(wt) for wt in weights])


def get_weights(A):
    """Computes the variable weights for the potential defined by A.

    The weight vector w satisfies A * w = d * epsilon, so
    w = d * A^{-1} * epsilon.

    Args:
        A: Square invertible integer matrix defining the potential.

    Returns:
        Integer vector of variable weights.
    """
    return get_d(A) * A.inverse() * vector([1] * A.nrows())
