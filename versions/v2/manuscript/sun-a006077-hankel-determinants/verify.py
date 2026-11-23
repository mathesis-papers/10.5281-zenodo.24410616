import json
from fractions import Fraction
from math import comb


def determinant(matrix):
    size = len(matrix)
    if size == 0:
        return 1
    entries = [list(row) for row in matrix]
    denominator = 1
    sign = 1
    for column in range(size - 1):
        pivot_row = next((row for row in range(column, size) if entries[row][column]), None)
        if pivot_row is None:
            return 0
        if pivot_row != column:
            entries[column], entries[pivot_row] = entries[pivot_row], entries[column]
            sign = -sign
        pivot = entries[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                numerator = entries[row][target] * pivot - entries[row][column] * entries[column][target]
                quotient, remainder = divmod(numerator, denominator)
                if remainder:
                    raise ArithmeticError('Nonintegral fraction-free elimination')
                entries[row][target] = quotient
            entries[row][column] = 0
        denominator = pivot
    return sign * entries[-1][-1]


def multinomial(index):
    return comb(2 * index, index) * comb(3 * index, index)


def coefficient(index, translation=3):
    return sum((-1) ** term * translation ** (index - 3 * term)
               * comb(index, 3 * term) * multinomial(term)
               for term in range(index // 3 + 1))


def normalized_coefficient(index):
    return comb(2 * index - 1, index - 1) * comb(3 * index - 1, index - 1)


def hankel(sequence, size, shift=0):
    return determinant([[sequence(row + column + shift) for column in range(size)]
                        for row in range(size)])


def rational_determinant(matrix):
    entries = [[Fraction(entry) for entry in row] for row in matrix]
    product = Fraction(1)
    for column in range(len(entries)):
        pivot_row = next((row for row in range(column, len(entries))
                          if entries[row][column]), None)
        if pivot_row is None:
            return Fraction(0)
        if pivot_row != column:
            entries[column], entries[pivot_row] = entries[pivot_row], entries[column]
            product = -product
        pivot = entries[column][column]
        product *= pivot
        for row in range(column + 1, len(entries)):
            multiplier = entries[row][column] / pivot
            for target in range(column + 1, len(entries)):
                entries[row][target] -= multiplier * entries[column][target]
    return product


def inertia(matrix):
    entries = [[Fraction(entry) for entry in row] for row in matrix]
    positive = negative = 0
    while entries:
        size = len(entries)
        diagonal = next((index for index in range(size) if entries[index][index]), None)
        if diagonal is not None:
            order = [diagonal] + [index for index in range(size) if index != diagonal]
            entries = [[entries[row][column] for column in order] for row in order]
            pivot = entries[0][0]
            positive += int(pivot > 0)
            negative += int(pivot < 0)
            entries = [[entries[row][column] - entries[row][0] * entries[0][column] / pivot
                        for column in range(1, size)] for row in range(1, size)]
        else:
            pair = next(((row, column) for row in range(size)
                         for column in range(row + 1, size) if entries[row][column]), None)
            if pair is None:
                return [positive, negative, size]
            order = list(pair) + [index for index in range(size) if index not in pair]
            entries = [[entries[row][column] for column in order] for row in order]
            pivot = entries[0][1]
            positive += 1
            negative += 1
            entries = [[entries[row][column] -
                        (entries[row][0] * entries[1][column] +
                         entries[row][1] * entries[0][column]) / pivot
                        for column in range(2, size)] for row in range(2, size)]
    return [positive, negative, 0]


def deformation(index, translation, scale):
    return sum(scale ** term * translation ** (index - 3 * term)
               * comb(index, 3 * term) * multinomial(term)
               for term in range(index // 3 + 1))


def sibling_checks():
    for index in range(17):
        assert hankel(lambda entry: comb(2 * entry, entry), index + 1) == 2 ** index
    for order in range(1, 7):
        harmonics = [Fraction(0)]
        for index in range(1, 21):
            harmonics.append(harmonics[-1] + Fraction(1, index ** order))
        for index in range(1, 11):
            matrix = [[harmonics[row + column] for column in range(index + 1)]
                      for row in range(index + 1)]
            assert (-1) ** index * rational_determinant(matrix) > 0
    bernoulli = [Fraction(1)]
    euler = [1]
    for index in range(1, 21):
        bernoulli.append(-sum(Fraction(comb(index + 1, term)) * bernoulli[term]
                              for term in range(index)) / (index + 1))
        euler.append(0 if index % 2 else
                     -sum(comb(index, term) * euler[term]
                          for term in range(0, index, 2)))
    for index in range(1, 11):
        bernoulli_matrix = [[bernoulli[row + column] ** 2
                             for column in range(index + 1)] for row in range(index + 1)]
        euler_matrix = [[euler[row + column] ** 2
                        for column in range(index + 1)] for row in range(index + 1)]
        assert rational_determinant(bernoulli_matrix) < 0
        assert rational_determinant(euler_matrix) > 0
    return {'order_two_franel': {'n_min': 0, 'n_max': 16},
            'harmonic': {'m_min': 1, 'm_max': 6, 'n_min': 1, 'n_max': 10},
            'bernoulli_and_euler_squares': {'n_min': 1, 'n_max': 10}}


def check():
    limit = 32
    rows = []
    for index in range(limit + 1):
        actual = hankel(coefficient, index + 1)
        block_size = index // 3
        unshifted = hankel(multinomial, block_size + 1)
        shifted_size = block_size if index % 3 == 0 else block_size + 1
        shifted = hankel(multinomial, shifted_size, 1)
        predicted = (0 if index % 3 == 1 else
                     (-1) ** ((index + 1) // 3) * unshifted * shifted ** 2)
        assert actual == predicted, ('residue factorization', index)
        row = {'n': index, 'W': str(actual)}
        if index % 3 != 1:
            quotient, remainder = divmod((-1) ** ((index + 1) // 3) * actual, 6 ** index)
            assert remainder == 0 and quotient > 0 and quotient % 2 == 1
            row['positive_odd_quotient'] = str(quotient)
        for translation in [-2, 0, 1, 3, 5]:
            transformed = hankel(lambda entry: coefficient(entry, translation), index + 1)
            assert transformed == actual, ('translation', index, translation)
        rows.append(row)
    deformations = []
    for index in range(21):
        block_size, residue = divmod(index, 3)
        expected = ([2 * block_size + 1, block_size, 0] if residue == 0 else
                    [2 * block_size + 1, block_size, 1] if residue == 1 else
                    [2 * block_size + 2, block_size + 1, 0])
        for translation in [-2, 0, 3]:
            for scale in [-2, -1, 0, 1, 2, 3]:
                matrix = [[deformation(row + column, translation, scale)
                           for column in range(index + 1)] for row in range(index + 1)]
                actual = determinant(matrix)
                predicted = (0 if residue == 1 else
                             scale ** (index * (index + 1) // 3) * int(rows[index]['W']))
                assert actual == predicted, ('deformation', index, translation, scale)
                expected_inertia = expected if scale else [1, 0, index]
                assert inertia(matrix) == expected_inertia, ('inertia', index, translation, scale)
        deformations.append({'n': index, 'inertia_for_nonzero_scale': expected})
    for index in range(1, 513):
        normalized = normalized_coefficient(index)
        assert multinomial(index) == 6 * normalized
        assert normalized % 2 == int(index & (index - 1) == 0)
    epsilon = lambda index: int(index > 0 and index & (index - 1) == 0)
    for size in range(65):
        assert hankel(epsilon, size, 1) % 2 == 1
        assert hankel(epsilon, size, 2) % 2 == 1
    for size in range(17):
        integral_unshifted = determinant([
            [multinomial(column) if row == 0 else normalized_coefficient(row + column)
             for column in range(size + 1)] for row in range(size + 1)])
        integral_shifted = hankel(normalized_coefficient, size, 1)
        assert hankel(multinomial, size + 1) == 6 ** size * integral_unshifted
        assert hankel(multinomial, size, 1) == 6 ** size * integral_shifted
        assert integral_unshifted > 0 and integral_unshifted % 2 == 1
        assert integral_shifted > 0 and integral_shifted % 2 == 1
    return {
        'arithmetic': 'exact integers; fraction-free determinant elimination',
        'coverage': {
            'literal_sum_hankel_and_residue_factorization': {'n_min': 0, 'n_max': limit},
            'translation_parameters': [-2, 0, 1, 3, 5],
            'coefficient_integrality_and_parity': {'r_min': 1, 'r_max': 512},
            'binary_shifted_hankel_determinants': {'size_min': 0, 'size_max': 64},
            'integral_normalizations': {'m_min': 0, 'm_max': 16},
            'deformation_and_exact_inertia': {'n_min': 0, 'n_max': 20,
                'translation_parameters': [-2, 0, 3], 'scale_parameters': [-2, -1, 0, 1, 2, 3]},
            'sibling_signs': sibling_checks(),
        },
        'status': 'all exact finite checks agree with the stated identities',
        'values': rows,
        'inertia': deformations,
    }


if __name__ == '__main__':
    print(json.dumps(check(), indent=2, sort_keys=True))
