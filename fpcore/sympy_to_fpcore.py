# transform output of Sympy to FPCores, also does a few basic optimizations that
# Herbie does not do for some reason, including:
# - (pow x 2) -> (* x x)
# - (pow x 4) -> (* (* x x) (* x x))
# - (pow x (/ 1 2)) -> (sqrt x)
# - (pow x (/ 1 4)) -> (sqrt (sqrt x))


import ast
from dataclasses import dataclass
from typing import Any
import sys


def codegen(tree: Any, title: str | None = None, names: set | None = None) -> str:
    if names is None:
        names = set()
    match tree:
        case ast.Expression(body):
            s = codegen(body, names=names)
            args = " ".join(names)
            return f'(FPCore ({args}) :precision binary64 :name "{title}" {s})'
        case ast.BinOp(left, ast.Pow(), ast.Constant(value=2)):
            # (pow x 2) -> (* x x)
            left_s = codegen(left, names=names)
            return f"(* {left_s} {left_s})"
        case ast.BinOp(left, ast.Pow(), ast.Constant(value=4)):
            # (pow x 4) -> (* (* x x) (* x x))
            left_s = codegen(left, names=names)
            return f"(* (* {left_s} {left_s}) (* {left_s} {left_s}))"
        case ast.BinOp(
            left,
            ast.Pow(),
            ast.BinOp(ast.Constant(value=1), ast.Div(), ast.Constant(value=2)),
        ):
            # (pow x (/ 1 2)) -> (sqrt x)
            return f"(sqrt {codegen(left, names=names)})"
        case ast.BinOp(
            left,
            ast.Pow(),
            ast.BinOp(ast.Constant(value=1), ast.Div(), ast.Constant(value=4)),
        ):
            # (pow x (/ 1 4)) -> (sqrt (sqrt (x))
            return f"(sqrt (sqrt {codegen(left, names=names)}))"
        case ast.BinOp(left, op, right):
            return f"({codegen(op, names=names)} {codegen(left, names=names)} {codegen(right, names=names)})"
        case ast.UnaryOp(op, operand):
            return f"({codegen(op, names=names)} {codegen(operand, names=names)})"
        case ast.Name(id, _):
            names.add(id)
            return id
        case ast.Call(ast.Name("sin"), [operand]):
            return f"(sin {codegen(operand, names=names)})"
        case ast.Call(ast.Name("cos"), [operand]):
            return f"(cos {codegen(operand, names=names)})"
        case ast.Call(ast.Name("atan2"), [op1, op2]):
            return f"(atan2 {codegen(op1, names=names)} {codegen(op2, names=names)})"
        case ast.Add():
            return "+"
        case ast.Sub():
            return "-"
        case ast.USub():
            return "-"
        case ast.Mult():
            return "*"
        case ast.Div():
            return "/"
        case ast.Pow():
            return "pow"
        case ast.Constant(value):
            return str(value)
        case _:
            raise Exception(ast.dump(tree))


def parse(s: str) -> ast.AST:
    return ast.parse(s, mode="eval")


if __name__ == "__main__":
    arg = sys.argv[1]
    tree = ast.parse(arg, mode="eval")
    print(codegen(tree))
