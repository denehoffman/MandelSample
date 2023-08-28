# Code from Boris Grube (he got it from somewhere else I think)
# https://github.com/bgrube/Moments/blob/main/typings/ROOT/__init__.pyi
"""
type stub file
see https://github.com/microsoft/pyright/blob/main/docs/type-stubs.md
"""

# tell type checker that it should allow access to any attribute within the ROOT module
# this silences the `is not a known member of module "ROOT"` errors
# see https://github.com/microsoft/pylance-release/issues/3216#issuecomment-1222512796
def __getattr__(name: str) -> Any: ...
