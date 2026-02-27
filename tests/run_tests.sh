#!/bin/bash

DIR="$(cd "$(dirname "$0")" && pwd)"

sage "$DIR/potentials_test.sage"
sage "$DIR/orbifold_test.sage"
sage "$DIR/diagonal_k3_test.sage"
sage "$DIR/deformed_diagonal_k3_test.sage"
