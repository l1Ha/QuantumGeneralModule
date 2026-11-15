"""
JSON experiment-record writer — Python mirror of `mod_json`
(src/mod_json.f90).

Produces byte-identical documents to the Fortran writer under the canonical
format contract (fixed by tests/fixtures/canonical_record.json and validated
by the stdlib `json` module in the tests):

  * 2-space indentation, one object/array element per line;
  * empty objects ``{}`` / arrays ``[]`` inline;
  * comma + newline between elements, no trailing comma;
  * numbers as ``d.dddddddddddddddddE±ddd`` (17 fractional digits, exactly
    three exponent digits; matches Fortran ``ES26.17E3`` trimmed);
  * strings escaped per RFC 8259 (``"``, ``\\``, ``\\b``, ``\\t``, ``\\n``,
    ``\\f``, ``\\r``, ``\\u00XX`` for other control characters).
"""
import math

__all__ = ["JsonWriter", "format_number"]


def format_number(v: float) -> str:
    """Canonical number formatting: ``f"{v:.17E}"`` with a 3-digit exponent."""
    s = f"{float(v):.17E}"
    mant, exp = s.split("E")
    sign = exp[0]
    digits = exp[1:].rjust(3, "0")
    return f"{mant}E{sign}{digits}"


def _escape(s: str) -> str:
    out = []
    for ch in str(s):
        code = ord(ch)
        if ch == '"':
            out.append('\\"')
        elif ch == "\\":
            out.append("\\\\")
        elif ch == "\b":
            out.append("\\b")
        elif ch == "\t":
            out.append("\\t")
        elif ch == "\n":
            out.append("\\n")
        elif ch == "\f":
            out.append("\\f")
        elif ch == "\r":
            out.append("\\r")
        elif code < 32:
            out.append(f"\\u00{code:02X}")
        else:
            out.append(ch)
    return "".join(out)


class JsonWriter:
    """Streaming JSON writer with the canonical formatting contract."""

    def __init__(self):
        self._parts = []
        self._depth = 0
        self._after_key = False
        self._has_elems = [False] * 128

    # ---------- low level ----------
    def _put(self, txt):
        self._parts.append(txt)

    def _indent(self, level):
        self._put("\n" + "  " * level)

    def _pre_value(self):
        if self._after_key:
            self._after_key = False
            return
        if not self._parts:
            return
        if self._has_elems[self._depth]:
            self._put(",")
        self._indent(self._depth)

    def _post_value(self):
        self._has_elems[self._depth] = True

    def _key(self, key):
        self._put('"' + _escape(key) + '": ')

    # ---------- containers ----------
    def begin_object(self, key=None):
        self._pre_value()
        if key is not None:
            self._key(key)
        self._put("{")
        self._depth += 1
        self._has_elems[self._depth] = False

    def end_object(self):
        had = self._has_elems[self._depth]
        self._depth -= 1
        if had:
            self._indent(self._depth)
        self._put("}")
        self._post_value()

    def begin_array(self, key=None):
        self._pre_value()
        if key is not None:
            self._key(key)
        self._put("[")
        self._depth += 1
        self._has_elems[self._depth] = False

    def end_array(self):
        had = self._has_elems[self._depth]
        self._depth -= 1
        if had:
            self._indent(self._depth)
        self._put("]")
        self._post_value()

    # ---------- values ----------
    def add_string(self, key, val):
        self._pre_value()
        if key is not None:
            self._key(key)
        self._put('"' + _escape(val) + '"')
        self._post_value()

    def add_int(self, key, val):
        self._pre_value()
        if key is not None:
            self._key(key)
        self._put(str(int(val)))
        self._post_value()

    def add_real(self, key, val):
        self._pre_value()
        if key is not None:
            self._key(key)
        self._put(format_number(val))
        self._post_value()

    def add_bool(self, key, val):
        self._pre_value()
        if key is not None:
            self._key(key)
        self._put("true" if val else "false")
        self._post_value()

    def add_null(self, key=None):
        self._pre_value()
        if key is not None:
            self._key(key)
        self._put("null")
        self._post_value()

    # ---------- arrays ----------
    def add_int_array(self, key, vals):
        self._pre_value()
        if key is not None:
            self._key(key)
        vals = list(vals)
        if not vals:
            self._put("[]")
        else:
            self._put("[")
            for i, v in enumerate(vals):
                if i:
                    self._put(",")
                self._indent(self._depth + 1)
                self._put(str(int(v)))
            self._indent(self._depth)
            self._put("]")
        self._post_value()

    def add_real_array(self, key, vals):
        self._pre_value()
        if key is not None:
            self._key(key)
        vals = list(vals)
        if not vals:
            self._put("[]")
        else:
            self._put("[")
            for i, v in enumerate(vals):
                if i:
                    self._put(",")
                self._indent(self._depth + 1)
                self._put(format_number(v))
            self._indent(self._depth)
            self._put("]")
        self._post_value()

    def add_real_array_2d(self, key, rows):
        """2D array: outer array with one inner array per row."""
        self._pre_value()
        if key is not None:
            self._key(key)
        rows = [list(r) for r in rows]
        if not rows or not rows[0]:
            self._put("[]")
        else:
            self._put("[")
            for i, row in enumerate(rows):
                if i:
                    self._put(",")
                self._indent(self._depth + 1)
                self._put("[")
                for j, v in enumerate(row):
                    if j:
                        self._put(",")
                    self._indent(self._depth + 2)
                    self._put(format_number(v))
                self._indent(self._depth + 1)
                self._put("]")
            self._indent(self._depth)
            self._put("]")
        self._post_value()

    # ---------- output ----------
    def get_text(self) -> str:
        return "".join(self._parts)
