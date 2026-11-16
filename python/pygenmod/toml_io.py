"""
TOML experiment-record writer — Python mirror of `mod_toml`
(src/mod_toml.f90), same canonical-format contract as `json_io`.

Produces byte-identical documents to the Fortran writer under the contract
fixed by tests/fixtures/canonical_record.toml and validated by tomllib/tomli:

  * ``key = value`` pairs, root-table keys first;
  * sections ``[table]`` / ``[table.sub]``; arrays of tables ``[[name]]``;
  * arrays as ``key = [`` + one element per line (2-space indent) + ``]``;
  * empty arrays inline (``key = []``);
  * numbers as ``d.dddddddddddddddddE±ddd`` (17 fractional digits, three
    exponent digits; matches Fortran ``ES26.17E3`` trimmed);
  * basic-string escaping per TOML 1.0 (same escape set as JSON).
"""
from .json_io import format_number, _escape

__all__ = ["TomlWriter"]


class TomlWriter:
    """Streaming TOML writer with the canonical formatting contract."""

    def __init__(self):
        self._parts = []

    # ---------- low level ----------
    def _put(self, txt):
        self._parts.append(txt)

    def _newline(self):
        if self._parts and not self._parts[-1].endswith("\n"):
            self._put("\n")

    def _key(self, key):
        self._newline()
        self._put(f"{key} = ")

    # ---------- tables ----------
    def begin_table(self, name):
        """Declare [name] (dotted names allowed); empty name = root table."""
        if not name:
            return
        self._newline()
        self._put(f"[{name}]\n")

    def begin_array_table(self, name):
        """Declare [[name]] (one element of an array of tables)."""
        self._newline()
        self._put(f"[[{name}]]\n")

    # ---------- values ----------
    def add_string(self, key, val):
        self._key(key)
        self._put('"' + _escape(str(val).strip()) + '"')

    def add_int(self, key, val):
        self._key(key)
        self._put(str(int(val)))

    def add_real(self, key, val):
        self._key(key)
        self._put(format_number(val))

    def add_bool(self, key, val):
        self._key(key)
        self._put("true" if val else "false")

    def _array(self, key, vals, fmt):
        self._key(key)
        vals = [fmt(v) for v in vals]
        if not vals:
            self._put("[]")
            return
        self._put("[")
        for i, v in enumerate(vals):
            if i:
                self._put(",")
            self._put("\n  " + v)
        self._put("\n]")

    def add_int_array(self, key, vals):
        self._array(key, vals, lambda v: str(int(v)))

    def add_real_array(self, key, vals):
        self._array(key, vals, format_number)

    def add_real_array_2d(self, key, rows):
        self._key(key)
        rows = [list(r) for r in rows]
        if not rows or not rows[0]:
            self._put("[]")
            return
        self._put("[")
        for i, row in enumerate(rows):
            if i:
                self._put(",")
            self._put("\n  [")
            for j, v in enumerate(row):
                if j:
                    self._put(",")
                self._put("\n    " + format_number(v))
            self._put("\n  ]")
        self._put("\n]")

    # ---------- output ----------
    def get_text(self) -> str:
        return "".join(self._parts)
