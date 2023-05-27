"""Parent class of Substance.

Class to redefine algebraic operators.
"""

from typing import Union

import numpy as np
from numpy.typing import NDArray

from sympy import Eq, Equality, Symbol, symbols


class Symbolic:
    """Parent class of Substance. Redefines algebraic operators.

    This class is created as a parent class of Substance and redefines the
    __add__, __mul__, __rmul__, __gt__ and __truediv__ methods to allow the
    construction of chemical equations from algebraic operations with Substance
    objects.

    Parameters
    ----------
    names : Union[str, np.ndarray]
        Array that stores the name attribute of the substances that
        participates in a chemical equation.
    expression : Symbol, optional
        Full chemical expression as a SymPy symbol, by default None.
    chem_equality : Equality, optional
        Full chemical expression as a SymPy Equality, by default None.

    Raises
    ------
    TypeError
        "Substance's name must be a string."
    """

    def __init__(
        self,
        names: Union[str, NDArray],
        expression: Symbol = None,
        chem_equality: Equality = None,
    ) -> None:
        # =====================================================================
        # Name storing:
        # =====================================================================
        if isinstance(names, str):
            self._names = np.array([names])
        elif isinstance(names, np.ndarray):
            self._names = names
        else:
            raise TypeError("Substance's name must be a string.")

        # =====================================================================
        # Expression:
        # =====================================================================
        if expression is None:
            self._expression = symbols(names)
        else:
            self._expression = expression

        # =====================================================================
        # Chem equality:
        # =====================================================================
        if chem_equality is None:
            self._chem_equality = None
        else:
            self._chem_equality = chem_equality

    def __rmul__(self, coeff: Union[int, float]) -> "Symbolic":
        """Redefine the reverse '*' operand.

        Redefines the reverse '*' operand to specify the stoichiometric of each
        substance.

        Parameters
        ----------
        coeff : Union[int, float]
            Stoichiometric coefficient.

        Returns
        -------
        Symbolic
            New Symbolic object with the respective stoichiometric coefficient.

        Raises
        ------
        TypeError
            "'*' operand defined only for 'int' or 'float' types."
        """
        if not isinstance(coeff, (int, float)):
            raise TypeError(
                "'*' operand defined only for 'int' or 'float' types."
            )

        new_sym = Symbolic(
            names=self._names,
            expression=coeff * self._expression,
        )
        return new_sym

    def __mul__(self, coeff: Union[int, float]) -> "Symbolic":
        """Redefine the '*' operand.

        Redefines the '*' operand to specify the stoichiometric of each
        substance.

        Parameters
        ----------
        coeff : Union[int, float]
            Stoichiometric coefficient.

        Returns
        -------
        Symbolic
            New Symbolic object with the respective stoichiometric coefficient.
        """
        return self.__rmul__(coeff)

    def __truediv__(self, coeff: Union[int, float]) -> "Symbolic":
        """Redefine the '/' operand.

        Redefines the '/' operand to specify the stoichiometric of each
        substance.

        Parameters
        ----------
        coeff : Union[int, float]
            Stoichiometric coefficient.

        Returns
        -------
        Symbolic
            New Symbolic object with the respective stoichiometric coefficient.

        Raises
        ------
        TypeError
            "'/' operand defined only for 'int' or 'float' types."
        """
        if not isinstance(coeff, (int, float)):
            raise TypeError(
                "'/' operand defined only for 'int' or 'float' types."
            )

        new_sym = Symbolic(
            names=self._names,
            expression=self._expression / coeff,
        )
        return new_sym

    def __add__(self, other_sym: "Symbolic") -> "Symbolic":
        """Redefine the '+' operand.

        Add two Symbolic objects to build a chemical expression.

        Parameters
        ----------
        other_sym : Symbolic
            Symbolic object to add.

        Returns
        -------
        Symbolic
            New Symbolic object that contains the sum of the two Symbolics, and
            stores the names and individual symbols of the addends.

        Raises
        ------
        TypeError
            "'+' defined only for Symbolic type."
        """
        if not isinstance(other_sym, Symbolic):
            raise TypeError("'+' defined only for Symbolic type.")

        new_sym = Symbolic(
            names=np.unique(np.append(self._names, other_sym._names)),
            expression=self._expression + other_sym._expression,
        )
        return new_sym

    def __gt__(self, other_sym: "Symbolic") -> "Symbolic":
        """Redefine the '>' operand.

        Construct a chemical expression.

        Parameters
        ----------
        other_sym : Symbolic
            The right-hand side of the chemical expression.

        Returns
        -------
        Symbolic
            New Symbolic object that contains the sum of the two Symbolics, and
            stores the names and individual symbols of the addends. Also
            creates the full chemical equation to print.

        Raises
        ------
        TypeError
            "'>' defined only for Symbolic type."
        """
        if not isinstance(other_sym, Symbolic):
            raise TypeError("'>' defined only for Symbolic type.")

        new_sym = Symbolic(
            names=np.append(self._names, other_sym._names),
            expression=-self._expression + other_sym._expression,
            chem_equality=Eq(self._expression, other_sym._expression),
        )
        return new_sym
