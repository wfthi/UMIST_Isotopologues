"""
umist_isotopologues.py -- Python port of the Mathematica package UMIST`
(source: isotopologue_umist.txt).

Reads a UMIST-style reaction network, generates all 13C/18O isotopologue
reactions, prunes unreasonable permutations with chemistry heuristics and
rescales rates for new product branches.

Row layout (0-based):  [n, t, a, b, c, d, e, f, g, alpha, beta, gamma, *extra]
Rows must have at least 13 columns.
"""

from __future__ import annotations

import math
import re
from itertools import permutations

__all__ = ["read_umist_csv", "insert_isotopes", "rescale_rates",
           "fortran_out", "frac_out", "eformat", "printf2", "pad_string"]

# --------------------------------------------------------------------------- #
# constants
# --------------------------------------------------------------------------- #

PHOTONS = ("CRPHOT", "CRP", "PHOTON")

# single-character placeholders for multi-char element symbols
# (Mathematica used double-struck letters; any reserved chars will do)
_HE, _SI, _CL = "\u210d", "\U0001d54a", "\u2102"
_NA, _MG, _FE = "\u2115", "\U0001d544", "\U0001d53d"
_C10, _C11 = "\U0001d54f", "\U0001d54e"

NORM2DOUBLE = (("C11", _C11), ("C10", _C10), ("Cl", _CL), ("Si", _SI),
               ("He", _HE), ("Na", _NA), ("Mg", _MG), ("Fe", _FE))
DOUBLE2NORM = tuple((b, a) for a, b in NORM2DOUBLE)

# element tokeniser used by insertOne (order matters: He before H, ...)
_INSERT_RE = re.compile("(" + "|".join(
    [_HE, "H", "C", _C10, _C11, "O", "F", "S", _SI, "N", "P",
     _NA, _MG, _FE, _CL, "X", "Y", "e"]) + r")(\d)?")

def _mask(s: str, pairs) -> str:
    for old, new in pairs:
        s = s.replace(old, new)
    return s

def _norm_to_double(s): return _mask(s, NORM2DOUBLE)
def _double_to_norm(s): return _mask(s, DOUBLE2NORM)

# masking sets used before regex group tests / counting (verbatim from source)
_CR = (("CRPHOT", "crphot"), ("CRP", "crp"), ("PHOTON", "photon"))
_MASK_CN = (("Cl", "cl"),) + _CR + (("Na", "na"),)                 # no Si (sic)
_MASK_CO = (("Si", "si"), ("Cl", "cl")) + _CR + (("Na", "na"),)
_MASK_CH = (("He", "he"), ("Cl", "cl")) + _CR + (("Na", "na"),)
_MASK_ISO = (("Si", "si"), ("Cl", "cl")) + _CR + (("Na", "na"),)

# --------------------------------------------------------------------------- #
# functional-group regexes (transcribed verbatim from the Mathematica source)
# --------------------------------------------------------------------------- #

RX_CN_BOUND = (r"\b(C|13C)N|(?<=(H\d))(C|13C)N[+\-]?\b|H[23]?(C|13C)N"
               r"|H[23]?N(C|13C)|NCCN|N(C13C|13CC)N|S(i)?N(C|13C)|O(C|13C)N")
RX_CN_ISOTOP = (r"\b13CN|(?<=(H\d))13CN[+\-]?\b|H[23]?13CN|H[23]?N13C"
                r"|N(C13C|13CC)N|S(i)?N13C|O13CN")
RX_CP_BOUND = (r"^(C|13C)(C|13C)P[+]?$|^(C|13C)P[+]?$|^H(C|13C)P[+]?$"
               r"|^(C|13C)H2PH$|^P(C|13C)H[234]{1}[+]?$")
RX_CP_ISOTOP = (r"^(C|13C)13CP[+]?$|^13CP[+]?$|^H13CP[+]?$"
                r"|^13CH2PH$|^P13CH[234]{1}[+]?$")
RX_CS_BOUND = r"CS|S[2-9]*(C|13C)\D*|CH[2-9]*S"
RX_CS_ISOTOP = (r"13C[2-9]*S|S[2-9]*13C|13C(H|H[2-9])*(D|D[2-9])*S"
                r"|13C(D|D[2-9])*(H|H[2-9])*S")
RX_CS_NOTISOTOP = r"(?<!13)CS|SC"
RX_CO_BOUND = (r"^(C|13C)H[2-3]{1}(C|13C)[H]{0,1}(O|18O)|^(C|13C)H[23]{1}(O|18O)H"
               r"|^(C|13C)(O|18O)|^H[23]?(C|13C)(O|18O)"
               r"|^H[23]{1}(C|13C)(C|13C)(O|18O)|^HN(C|13C)(O|18O)[+]?$"
               r"|^H(O|18O)(C|13C)|^(O|18O)(C|13C)"
               r"|^(C|13C)H3(O|18O)(C|13C)H[34]{1}")
RX_C18O_BOUND = (r"^(C|13C)H[2-3]{1}(C|13C)[H]{0,1}18O|^(C|13C)H[23]{1}18OH"
                 r"|^(C|13C)18O|^H[23]?(C|13C)18O|^H[23]{1}(C|13C)(C|13C)18O"
                 r"|^HN(C|13C)18O[+]?$|^H18O(C|13C)|^18O(C|13C)"
                 r"|^H[23]?(C|13C)(O18O|18OO|18O18O)(C|13C)"
                 r"|^(C|13C)(O18O|18OO|18O18O)(C|13C)"
                 r"|^H(C|13C)O18OH|^H(C|13C)O18O"
                 r"|^(C|13C)H318O(C|13C)H[34]{1}|^(C|13C)O18O[+]?$")
# note: the original isNotC18Obound lost a "|" before the last alternative
RX_C18O_NOTBOUND = (r"^(C|13C)H[2-3]{1}(C|13C)[H]{0,1}18O|^(C|13C)H[23]{1}18OH"
                    r"|^(C|13C)18O|^H[23]?(C|13C)18O|^H[23]{1}(C|13C)(C|13C)18O"
                    r"|^HN(C|13C)18O[+]?$|^H18O(C|13C)|^18O(C|13C)"
                    r"|^H[23]?(C|13C)(O18O|18OO|18O18O)(C|13C)"
                    r"|^(C|13C)(O18O|18OO|18O18O)(C|13C)"
                    r"|^H(C|13C)O18OH|^H(C|13C)O18O"
                    r"|^(C|13C)H318O(C|13C)H[34]{1}^(C|13C)O18O[+]?$")
RX_CO_ISOTOP = (r"^(C|13C)H[2-3]{1}13C[H]{0,1}(O|18O)|^13CH[23]{1}(O|18O)H"
                r"|^13C(O|18O)|^H[23]?13C(O|18O)"
                r"|^H[23]{1}(C|13C)13C(O|18O)|^HN13C(O|18O)[+]?$"
                r"|^H(O|18O)13C|^(O|18O)13C"
                r"|^H[23]?C(O|18O)(O|18O)13C|^C(O|18O)(O|18O)13C"
                r"|^13CH3(O|18O)(C|13C)H[34]{1}|^(C|13C)H3(O|18O)13CH[34]{1}")
RX_CO_NOTISOTOP = (r"^(C|13C)H[2-3]{1}13C[H]{0,1}(O|18O)|^13CH[23]{1}(O|18O)H"
                   r"|^13C(O|18O)|^H[23]?13C(O|18O)"
                   r"|^H[23]{1}(C|13C)13C(O|18O)|^HN13C(O|18O)[+]?$"
                   r"|^H(O|18O)13C|^(O|18O)13C"
                   r"|^H[23]?C(O|18O)(O|18O)13C|^C(O|18O)(O|18O)13C")
RX_13C = r"13C"
RX_18O = r"18O"
RX_LEFT18O = r"^H(C|13C)18OO|^(C|13C)18OO"
RX_RIGHT18O = r"^H(C|13C)O18O|^(C|13C)O18O"
RX_LEFT13CN = r"^N13CCN"
RX_CH_BOUND = r"C[2-9]*H[2-5]*|\b(?!=N)(?<!C[2-9]{1})H[2-5]*(C|13C)(?![2-9]{1})"
RX_CH_ISOTOP = r"13CH[2-5]*|^H[2-5]*13C(?![2-9]{1})"
RX_CH_NOTISOTOP = r"(C|13C)H[2-5]*|^H[2-5]*(C|13C)(?![2-9]{1})"
RX_CHN_BOUND = (r"^(C|13C)H[2-5]{1,1}|[abd-z|ABD-Z](13C|C)H[2-4]{1,1}[+]?$"
                r"|^(13C|C)H3(13C|C)H3|^H[2-3]{1}C(?!\d{1,1})"
                r"|^H2(C13C|13CC)|^H[2-3]{1}13C[abd-zABD-Z]"
                r"|^H[2-3]{1}(C18O|13C18O)")
RX_CHN_ISOTOP = (r"^13CH[2-5]{1,1}|[abd-z|ABD-Z]13CH[2-4]{1,1}[+]?$"
                 r"|^13CH3CH3|^CH313CH3|^H[2-3]{1}13C[abd-zABD-Z]"
                 r"|^H[2-3]{1}13C18O|^H[2-3]{1}13CCl"
                 r"|^H[2-3]{1}13C(C|13C)(O|18O)[+]?"
                 r"|^H[2-3]{1}13C(C|13C)(C|13C)(C|13C)?[+]?")
RX_COG_MAIN = r"(C|13C)H[2-3]*(O|18O)|(C|13C)(O|18O)|(O|18O)(C|13C)"
RX_COG_S1 = r"H3(18O|O)(13C|C)H[34]"
RX_COG_S2 = r"^[H]?(13C|C)O18O[+-]?$|^[H]?(13C|C)O2[+-]?$"
RX_COG_S3 = (r"^H(13C|C)O18OH[2]?[+-]?$|^H(13C|C)18OOH[2]?[+-]?$"
             r"|^H(13C|C)OOH[2]?[+-]?$")
RX_COG_S4 = r"^H(13C|C)(O18O|18OO|OO)(C|13C)H3$"
RX_COG_NEG = (r"^C13C|^C[2-9]?13C|^H[35]{1}C[2-9]{1}|^H[35]{1}C13C"
              r"|^HC2(O|18O)|^HC13C(O|18O)|^HC3(O|18O)|^HC213C(O|18O)")
RX_CHNG_S1 = r"^(C|13C)H3(C|13C)H3$"
RX_STICKY13C = r"^H2(C|13C)(C|13C)(C|13C)$|^H2(C|13C)(C|13C)(C|13C)(C|13C)$"
RX_STICKY18O = r"^H?(C|13C)(O|18O)(O|18O)(C|13C)H[34]{1}[+]?$"
RX_STICKY18O_CORE = r"(C|13C)(O|18O)(O|18O)(C|13C)"

# --------------------------------------------------------------------------- #
# small logic helpers (Mathematica n-ary logic / chained equalities)
# --------------------------------------------------------------------------- #

def implies(p, q): return (not p) or q
def _when(cond, val): return val if cond else True
def any_of(xs): return any(bool(x) for x in xs)
def all_of(xs): return all(bool(x) for x in xs)
def xor_of(xs): return sum(1 for x in xs if x) % 2 == 1
def nand_of(xs): return not all_of(xs)
def nor_of(xs): return not any_of(xs)
def mapthread_and(xs, ys): return [bool(x) and bool(y) for x, y in zip(xs, ys)]

# --------------------------------------------------------------------------- #
# species string helpers
# --------------------------------------------------------------------------- #

def is_ion(s: str) -> bool:
    return bool(s) and s[-1] in "+-"

def charge_of(s: str) -> str:
    return s[-1] if is_ion(s) else ""

def neutral(s: str) -> str:
    return s[:-1] if is_ion(s) else s

def strip_charge(s: str) -> str:
    return s[:-1] if len(s) > 1 and s[-1] in "+-" else s

def on_grain(s: str) -> bool:
    return bool(s) and s[0] == "J"

def from_j(s: str) -> str:
    return s[1:] if s.startswith("J") else s

def to_j(s: str) -> str:
    return s if s.startswith("J") else "J" + s

# --------------------------------------------------------------------------- #
# atom counting (algorithm identical to the original, including its quirks)
# --------------------------------------------------------------------------- #

def _count_symbol(s: str, sym: str, c10_special: bool = False) -> int:
    count = 0
    if s.endswith(sym):                     # trailing element without number
        count += 1
        s = s[:-1]
    if c10_special:                         # only countCatoms has this line
        count += 9 * s.count("10")
    i = s.find(sym)
    while i >= 0:
        nxt = s[i + 1] if i + 1 < len(s) else ""
        count += int(nxt) if nxt.isdigit() else 1   # ONE digit only (sic)
        i = s.find(sym, i + 1)
    return count

def count_hatoms(s):
    return _count_symbol(_mask(s, (("He", _HE), ("Cl", "cl")) + _CR), "H")

def count_heatoms(s):
    return _count_symbol(_mask(s, (("He", _HE), ("Si", "si"), ("Cl", "cl")) + _CR), _HE)

def count_catoms(s):
    return _count_symbol(_mask(s, (("Cl", "cl"), ("CRPHOT", "crphot"), ("CRP", "crp"))),
                         "C", c10_special=True)

def count_natoms(s):   # no Na mask in the original (quirk kept)
    return _count_symbol(_mask(s, (("Cl", "cl"),) + _CR), "N")

def count_oatoms(s):
    return _count_symbol(_mask(s, (("Si", "si"), ("Cl", "cl")) + _CR), "O")

def count_fatoms(s):
    return _count_symbol(_mask(s, (("Si", "si"), ("Cl", "cl")) + _CR), "F")

def count_naatoms(s):
    return _count_symbol(_mask(s, (("Si", _SI), ("Cl", _CL), ("He", _HE),
                                    ("Mg", _MG), ("Na", _NA), ("Fe", _FE)) + _CR), _NA)

def count_mgatoms(s):
    return _count_symbol(_mask(s, (("Si", _SI), ("Cl", _CL), ("He", _HE),
                                    ("Mg", _MG), ("Na", _NA), ("Fe", _FE)) + _CR), _MG)

def count_siatoms(s):
    return _count_symbol(_mask(s, (("Si", _SI), ("Cl", "cl")) + _CR), _SI)

def count_patoms(s):
    return _count_symbol(_mask(s, (("Si", "si"), ("Cl", "cl")) + _CR), "P")

def count_satoms(s):   # no Si mask in the original (quirk kept)
    return _count_symbol(_mask(s, (("Cl", "cl"),) + _CR), "S")

def count_clatoms(s):
    return _count_symbol(_mask(s, (("Si", "\u03a3"), ("Cl", _CL)) + _CR), _CL)

def count_featoms(s):
    return _count_symbol(_mask(s, (("Si", _SI), ("Cl", _CL), ("He", _HE),
                                    ("Mg", _MG), ("Na", _NA), ("Fe", _FE)) + _CR), _FE)

def count_elements(s):
    """Inventory in the order {H,He,C,N,O,F,Na,Mg,Si,P,S,Cl,Fe}."""
    return [count_hatoms(s), count_heatoms(s), count_catoms(s), count_natoms(s),
            count_oatoms(s), count_fatoms(s), count_naatoms(s), count_mgatoms(s),
            count_siatoms(s), count_patoms(s), count_satoms(s),
            count_clatoms(s), count_featoms(s)]

def hcnos(s):
    """{H, He, C, N, O, Si, P, S, Cl} -- the 9-vector used by all deltas."""
    return [count_hatoms(s), count_heatoms(s), count_catoms(s), count_natoms(s),
            count_oatoms(s), count_siatoms(s), count_patoms(s),
            count_satoms(s), count_clatoms(s)]

def is_atom(s):
    return sum(hcnos(s)) == 1

# --------------------------------------------------------------------------- #
# transfer vectors / partner search
# --------------------------------------------------------------------------- #

def _vsub(u, v): return [a - b for a, b in zip(u, v)]
def _vadd(u, v): return [a + b for a, b in zip(u, v)]
def mixed_signs(v): return any(x > 0 for x in v) and any(x < 0 for x in v)
def total_abs(v): return sum(abs(x) for x in v)
def _norm(v): return math.sqrt(sum(x * x for x in v))
def _sgn(x): return 1 if x > 0 else 0       # Sign of a non-negative count

def delta_ad(a, b, d, e): return _vsub(hcnos(a), hcnos(d))
def delta_ae(a, b, d, e): return _vsub(hcnos(a), hcnos(e))
def delta_bd(a, b, d, e): return _vsub(hcnos(b), hcnos(d))
def delta_be(a, b, d, e): return _vsub(hcnos(b), hcnos(e))

def delta_ax(a, b, d, e):
    ae, ad = _vsub(hcnos(a), hcnos(e)), _vsub(hcnos(a), hcnos(d))
    if not mixed_signs(ae) and not mixed_signs(ad) and total_abs(ae) < total_abs(ad):
        return ae
    if not mixed_signs(ae) and not mixed_signs(ad) and total_abs(ae) > total_abs(ad):
        return ad
    if not mixed_signs(ae) and mixed_signs(ad): return ae
    if not mixed_signs(ad) and mixed_signs(ae): return ad
    return []

def delta_bx(a, b, d, e):
    be, bd = _vsub(hcnos(b), hcnos(e)), _vsub(hcnos(b), hcnos(d))
    if not mixed_signs(be) and not mixed_signs(bd) and total_abs(be) < total_abs(bd):
        return be
    if not mixed_signs(be) and not mixed_signs(bd) and total_abs(be) > total_abs(bd):
        return bd
    if not mixed_signs(be) and mixed_signs(bd): return be
    if not mixed_signs(bd) and mixed_signs(be): return bd
    return []

def find_partner_a(a, b, d, e):
    ae, ad = _vsub(hcnos(a), hcnos(e)), _vsub(hcnos(a), hcnos(d))
    if not mixed_signs(ae) and not mixed_signs(ad) and total_abs(ae) < total_abs(ad): return 2
    if not mixed_signs(ae) and not mixed_signs(ad) and total_abs(ae) > total_abs(ad): return 1
    if not mixed_signs(ae) and mixed_signs(ad): return 2
    if not mixed_signs(ad) and mixed_signs(ae): return 1
    return 0

def find_partner_b(a, b, d, e):
    be, bd = _vsub(hcnos(b), hcnos(e)), _vsub(hcnos(b), hcnos(d))
    if not mixed_signs(be) and not mixed_signs(bd) and total_abs(be) < total_abs(bd): return 2
    if not mixed_signs(be) and not mixed_signs(bd) and total_abs(be) > total_abs(bd): return 1
    if not mixed_signs(be) and mixed_signs(bd): return 2
    if not mixed_signs(bd) and mixed_signs(be): return 1
    return 0

def _match_delta(vec, spec, length=9):
    """MatchQ[vec, {_,_,0,...}] style: exact length + fixed positions.
    NOTE: the source contains 8-slot patterns for 9-vectors; they never
    match, which is reproduced here via the length check."""
    if len(vec) != length:
        return False
    return all(vec[i] == val for i, val in spec.items())

# --------------------------------------------------------------------------- #
# functional group tests
# --------------------------------------------------------------------------- #

def _hit(rx, s, masks):
    return re.search(rx, _mask(s, masks)) is not None

def has_13c(s):        return _hit(RX_13C, s, _MASK_ISO)
def has_18o(s):        return _hit(RX_18O, s, _MASK_ISO)
def cn_bound(s):       return _hit(RX_CN_BOUND, s, _MASK_CN)
def cn_isotop(s):      return _hit(RX_CN_ISOTOP, s, _MASK_CN)
def no_cn_isotop(s):   return not _hit(RX_CN_ISOTOP, s, _MASK_CN)     # StringFreeQ
def cp_bound(s):       return _hit(RX_CP_BOUND, s, _MASK_CN)
def cp_isotop(s):      return _hit(RX_CP_ISOTOP, s, _MASK_CN)
def cs_bound(s):       return _hit(RX_CS_BOUND, s, _MASK_CO)
def cs_isotop(s):      return _hit(RX_CS_ISOTOP, s, _MASK_CO)
def no_cs_isotop(s):   return _hit(RX_CS_NOTISOTOP, s, _MASK_CO)    # Not@StringFreeQ (sic)
def co_bound(s):       return _hit(RX_CO_BOUND, s, _MASK_CO)
def co_isotop(s):      return _hit(RX_CO_ISOTOP, s, _MASK_CO)
def no_co_isotop(s):   return not _hit(RX_CO_ISOTOP, s, _MASK_CO)
def c18o_bound(s):     return _hit(RX_C18O_BOUND, s, _MASK_CO)
def ch_bound(s):       return _hit(RX_CH_BOUND, s, _MASK_CH)
def ch_isotop(s):      return _hit(RX_CH_ISOTOP, s, _MASK_CH)
def no_ch_isotop(s):   return not _hit(RX_CH_NOTISOTOP, s, _MASK_CH)
def chn_bound(s):      return _hit(RX_CHN_BOUND, s, _MASK_CH)
def chn_isotop(s):     return _hit(RX_CHN_ISOTOP, s, _MASK_CH)
def no_chn_isotop(s):  return not _hit(RX_CHN_ISOTOP, s, _MASK_CH)
def left18o_bound(s):  return _hit(RX_LEFT18O, s, _MASK_CO)
def right18o_bound(s): return _hit(RX_RIGHT18O, s, _MASK_CO)
def left13cn_bound(s): return _hit(RX_LEFT13CN, s, _MASK_CO)

def n_cogroups(s):
    ins = _mask(s, (("He", "he"), ("Cl", "cl")) + _CR + (("Na", "na"),))
    count = len(re.findall(RX_COG_MAIN, ins))          # Overlaps -> False
    if re.search(RX_COG_S1, ins): count += 1
    if re.search(RX_COG_S2, ins): count += 1
    if re.search(RX_COG_S3, ins): count += 1
    if re.search(RX_COG_S4, ins): count += 1
    if re.search(RX_COG_NEG, ins): count = 0
    return count

def n_chngroups(s):
    ins = _mask(s, (("He", "he"), ("Cl", "cl")) + _CR + (("Na", "na"),))
    count = len(re.findall(RX_CHN_BOUND, ins))
    if re.search(RX_CHNG_S1, ins):
        count += 1
    return count

def coch_bound(s):
    if not chn_bound(s):
        return False
    return n_cogroups(s) + n_chngroups(s) > count_catoms(s)

def sticky_13c(species):
    rx = re.compile(RX_STICKY13C)
    pos = [i for i, s in enumerate(species) if rx.match(s)]
    if len(pos) == 2:
        return species[pos[0]] == species[pos[1]]
    return True

def sticky_18o(species):
    rx, core = re.compile(RX_STICKY18O), re.compile(RX_STICKY18O_CORE)
    pos = [i for i, s in enumerate(species) if rx.match(s)]
    if len(pos) == 2:
        m1 = [m.group(0) for m in core.finditer(species[pos[0]])]
        m2 = [m.group(0) for m in core.finditer(species[pos[1]])]
        return m1 == m2
    return True

def chemical_palindrome(word):
    rest = neutral(word)
    n = len(rest)
    if n >= 6:
        return rest[:3] == rest[-3:] or rest[:3] == rest[-3:][::-1]
    if n >= 4:
        return rest[:2] == rest[-2:] or rest[:2] == rest[-2:][::-1]
    if n == 3:
        return rest[0] == rest[-1]
    return False

# --------------------------------------------------------------------------- #
# isotopologue generation
# --------------------------------------------------------------------------- #

# slot indices within a row
N, T, A, B, C, D, E, F, G, ALPHA, BETA, GAMMA = range(12)

def insert_one(x: str) -> str:
    """Enumerate every element symbol with an explicit count (x -> x1)."""
    return "".join(sym + (num if num else "1")
                   for sym, num in _INSERT_RE.findall(x))

def _string_replace_list(s, rules):
    """All distinct results of applying one rule at one position
    (StringReplaceList semantics for our use case); [s] if nothing matches."""
    results = []
    for pat, rep in rules:
        start = 0
        while True:
            idx = s.find(pat, start)
            if idx < 0:
                break
            results.append(s[:idx] + rep + s[idx + len(pat):])
            start = idx + 1
    if not results:
        return [s]
    seen, out = set(), []
    for r in results:
        if r not in seen:
            seen.add(r)
            out.append(r)
    return out

RULE_C = [("C2", "C1X1"), ("C3", "C2X1"), ("C4", "C3X1"), ("C5", "C4X1"),
          ("C6", "C5X1"), ("C7", "C6X1"), ("C8", "C7X1"), ("C9", "C8X1"),
          (_C10 + "1", "C9X1"), (_C11 + "1", _C10 + "X1"), ("C1", "X1")]  # X = 13C
RULE_O = [("O2", "O1Y1"), ("O1", "Y1")]                                     # Y = 18O
RULE_H = [("H2", "H1D1"), ("H3", "H2D1"), ("H4", "H3D1"), ("H5", "H4D1"),
          ("H6", "H5D1"), ("H7", "H6D1"), ("H1", "D1")]

def carbonize_single(x):
    return [v.replace("1", "").replace("X", "13C")
            for v in _string_replace_list(insert_one(x), RULE_C)]

def oxygenize_single(x):
    return [v.replace("1", "").replace("Y", "18O")
            for v in _string_replace_list(insert_one(x), RULE_O)]

def carbonize_chain(x):     # keeps the X marker (cleaned later)
    return [v.replace("1", "") for v in _string_replace_list(insert_one(x), RULE_C)]

def deuterate(x):
    return [v.replace("1", "") for v in _string_replace_list(insert_one(x), RULE_H)]

def _with_charge_and_grain(variants, original):
    out = [v + charge_of(original) for v in variants]
    if on_grain(original):
        out = [to_j(v) for v in out]
    return out

def carbonize_plus_single(y):  return _with_charge_and_grain(carbonize_single(y), y)
def oxygenize_plus_single(y):  return _with_charge_and_grain(oxygenize_single(y), y)
def carbonize_plus_chain(y):   return _with_charge_and_grain(carbonize_chain(y), y)

# second-isotope machinery (insertOne2 / ruleC2 / carbonize2 in the original)
_INSERT_RE2 = re.compile("(He|D|H|C|O|S|Si|N|P|Na|Mg|Fe|Cl|X|Y|e)(\\d)?")
RULE_C2 = [("C1X1", "X2"), ("C2X1", "C1X2"), ("C3X1", "C2X2"), ("C4X1", "C3X2"),
           ("C5X1", "C4X2"), ("C6X1", "C5X2"), ("C7X", "C6X2"), ("C8X1", "C7X2"),
           ("C2", "C1X1"), ("C3", "C2X1"), ("C4", "C3X1"), ("C5", "C4X1"),
           ("C6", "C5X1"), ("C7", "C6X1"), ("C8", "C7X1"), ("C9", "C8X1"),
           ("C1", "X1")]

def insert_one2(x):
    return "".join(sym + (num if num else "1")
                   for sym, num in _INSERT_RE2.findall(_mask(x, (("13C", "X"), ("18O", "Y")))))

def carbonize2(x):
    variants = [v.replace("1", "") for v in _string_replace_list(insert_one2(x), RULE_C2)]
    return [to_j(v) if on_grain(x) else v for v in variants]

# ---- the 13 slot-pair rules ----------------------------------------------- #

SLOT_PAIRS = {1: (A, D), 2: (A, E), 3: (A, F), 4: (A, G),
              5: (B, D), 6: (B, E), 7: (B, F), 8: (B, G),
              9: (C, D), 10: (C, E), 11: (C, F), 12: (C, G)}

def _rule_applies(i, row, fkt):
    ia, ib = SLOT_PAIRS[i]
    va, vb = row[ia], row[ib]
    if va in PHOTONS or vb in PHOTONS:
        return False
    if not fkt(va) or not fkt(vb):
        return False
    if i == 2:  return row[D] != row[E]                       # d =!= e
    if i == 3:  return row[E] != row[F]                       # e =!= f
    if i == 5:  return row[A] != row[B]                       # a =!= b
    if i == 6:  return row[D] != row[E] and row[A] != row[B]
    if i == 7:  return row[E] != row[F] and row[A] != row[B]
    if i == 8:  return row[A] != row[B]
    if i == 10: return row[C] != row[B]                       # c =!= b (sic)
    if i == 11: return row[E] != row[F] and row[C] != row[B]
    if i == 12: return row[C] != row[B]
    return True                                               # rules 1,4,9

def _apply_rule(i, row, fkt):
    if not _rule_applies(i, row, fkt):
        return []
    ia, ib = SLOT_PAIRS[i]
    out = []
    for x in fkt(row[ia]):
        for y in fkt(row[ib]):
            new = list(row)
            new[ia], new[ib] = x, y
            out.append(new)
    return out

def _dedupe(rows):
    seen, out = set(), []
    for r in rows:
        key = tuple(r)
        if key not in seen:
            seen.add(key)
            out.append(list(key))
    return out

def _row_key(r):
    return tuple(str(v) for v in r)

def _expand_family(row, fkt):
    rows = [list(row)]                                        # rule 0: identity
    for i in range(1, 13):
        rows.extend(_apply_rule(i, row, fkt))
    return _dedupe(rows)

def _clean_rule_val(v):
    return v.replace("X", "13C").replace("Y", "18O") if isinstance(v, str) else v

def insert_13c_o18(rows):
    """in13C18O: 13C chain pass, then 18O pass, then X/Y -> 13C/18O."""
    stage1 = []
    for row in rows:
        stage1.extend(_expand_family(row, carbonize_plus_chain))
    stage1 = _dedupe(stage1)
    stage2 = []
    for row in stage1:
        stage2.extend(_expand_family(row, oxygenize_plus_single))
    return [[_clean_rule_val(v) for v in row] for row in _dedupe(stage2)]

# --------------------------------------------------------------------------- #
# small predicate helpers shared by the clean-up functions
# --------------------------------------------------------------------------- #

def _perm_joins(species):
    return {"".join(p) for p in permutations(species)}

def _neutral_perm_joins(species):
    return _perm_joins([neutral(s) for s in species])

def _charged_perm_joins(species):
    return {j + "+" for j in _neutral_perm_joins(species)}

def _ext_pref(s, prefix):
    """startswith + strictly longer + next char not a digit
    ('CO' is not a prefix of 'CO2')."""
    return (s.startswith(prefix) and len(s) > len(prefix)
            and not s[len(prefix)].isdigit())

def _ext_pref_loose(s, prefix):
    return s.startswith(prefix) and len(s) > len(prefix)

def _ch3_shrink_match(src, dst):
    """reference-line check of the CH3->CH2 cases (Min[4,len(src),len(dst)])."""
    mod = src[:2] + "2" + src[3:]
    k = min(4, len(src), len(dst))
    return mod[:k] == dst[:k]

def _any3(pred, *slots):
    return any_of([pred(s) for s in slots])

def _co_variant(base):
    """(first-4/last-n overlap) && (any C18O == any) && (any COisotop == any)"""
    def pred(a, b, c, d, e, f, g):
        return (base(a, b, c, d, e, f, g)
            and _any3(c18o_bound, a, b, c) == _any3(c18o_bound, d, e, f)
            and _any3(co_isotop, a, b, c) == _any3(co_isotop, d, e, f))
    return pred

def _c18o_together(a, b, c, d, e, f):
    """13C and 18O stay together as a 13C18O group."""
    return (_any3(co_isotop, a, b, c) == _any3(co_isotop, d, e, f)
        and any_of(mapthread_and([co_isotop(s) for s in (a, b, c)],
                                 [c18o_bound(s) for s in (a, b, c)]))
        == any_of(mapthread_and([co_isotop(s) for s in (d, e, f)],
                                [c18o_bound(s) for s in (d, e, f)])))

def _proton_transfer_pred(a, b, c, d, e, f, g, ae_is_unit):
    """shared by cleanUpIN and cleanUpNN ('H/Proton-Transfer' case)"""
    ok = sticky_13c([a, b, d, e])
    if ae_is_unit:                      # a<->e , b<->d keep their atoms
        ok &= _when(cn_bound(b) == cn_bound(d),
                    cn_isotop(b) == cn_isotop(d) and left13cn_bound(b) == left13cn_bound(d))
        ok &= _when(cn_bound(a) == cn_bound(e),
                    cn_isotop(a) == cn_isotop(e) and left13cn_bound(a) == left13cn_bound(e))
        ok &= _when(co_bound(b) == co_bound(d), left18o_bound(b) == left18o_bound(d))
        ok &= _when(co_bound(a) == co_bound(e), left18o_bound(a) == left18o_bound(e))
        cob_b, cob_d = coch_bound(b), coch_bound(d)
        ok &= _when(cob_b == cob_d and cob_d is False, co_isotop(b) == co_isotop(d))
        # original: (cob_b != cob_d) == False  ->  cob_b == cob_d (left-assoc parse)
        ok &= _when(cob_b == cob_d, implies(co_isotop(d), co_isotop(b)))
        cob_a, cob_e = coch_bound(a), coch_bound(e)
        ok &= _when(cob_a == cob_e and cob_e is False, co_isotop(a) == co_isotop(e))
        ok &= _when(cob_a == cob_e, implies(co_isotop(a), co_isotop(e)))
        ok &= has_13c(b) == has_13c(d) and has_13c(a) == has_13c(e)
        ok &= has_18o(b) == has_18o(d) and has_18o(a) == has_18o(e)
        ok &= _when(chn_bound(a) == chn_bound(e) and coch_bound(a) == coch_bound(e),
                    chn_isotop(a) == chn_isotop(e))
        ok &= _when(chn_bound(b) == chn_bound(d) and coch_bound(b) == coch_bound(d),
                    chn_isotop(b) == chn_isotop(d))
    else:                               # a<->d , b<->e
        ok &= _when(cn_bound(b) == cn_bound(e),
                    cn_isotop(b) == cn_isotop(e) and left13cn_bound(b) == left13cn_bound(e))
        ok &= _when(cn_bound(a) == cn_bound(d),
                    cn_isotop(a) == cn_isotop(d) and left13cn_bound(a) == left13cn_bound(d))
        ok &= _when(co_bound(b) == co_bound(e), left18o_bound(b) == left18o_bound(e))
        ok &= _when(co_bound(a) == co_bound(d), left18o_bound(a) == left18o_bound(d))
        ok &= co_isotop(b) == co_isotop(e) and co_isotop(a) == co_isotop(d)
        ok &= has_13c(b) == has_13c(e) and has_13c(a) == has_13c(d)
        ok &= has_18o(b) == has_18o(e) and has_18o(a) == has_18o(d)
        ok &= _when(chn_bound(a) == chn_bound(d), chn_isotop(a) == chn_isotop(d))
        ok &= _when(chn_bound(b) == chn_bound(e), chn_isotop(b) == chn_isotop(e))
    return bool(ok)

# ---- predicates shared verbatim between IN and NN ------------------------- #

def _pred_0b(a, b, c, d, e, f, g):
    return ((_any3(co_isotop, a, b, c) and xor_of([co_isotop(s) for s in (d, e, f)]))
         or (_any3(cn_isotop, a, b, c) and xor_of([cn_isotop(s) for s in (d, e, f)]))
         or (all_of([no_co_isotop(s) for s in (a, b, c)])
             and all_of([no_cn_isotop(s) for s in (a, b, c)])
             and all_of([no_co_isotop(s) for s in (d, e, f)])
             and all_of([no_cn_isotop(s) for s in (d, e, f)])))

def _pred_0b1(a, b, c, d, e, f, g):
    cond18 = ((any_of([c18o_bound(s) for s in (a, b, c)])
               and any_of([c18o_bound(s) for s in (d, e, f)]))
           or (nor_of([c18o_bound(s) for s in (a, b, c)])
               and nor_of([c18o_bound(s) for s in (d, e, f)])))
    return (((_any3(co_isotop, a, b, c) and xor_of([co_isotop(s) for s in (d, e, f)])) and cond18)
         or ((_any3(cn_isotop, a, b, c) and xor_of([cn_isotop(s) for s in (d, e, f)])) and cond18)
         or ((all_of([no_co_isotop(s) for s in (a, b, c)])
              and all_of([no_cn_isotop(s) for s in (a, b, c)])
              and all_of([no_co_isotop(s) for s in (d, e, f)])
              and all_of([no_cn_isotop(s) for s in (d, e, f)])) and cond18))

def _pred_0b2(a, b, c, d, e, f, g):
    return ((_any3(co_isotop, a, b, c) and xor_of([co_isotop(s) for s in (d, e, f)]))
         or (_any3(cn_isotop, a, b, c) and all_of([no_co_isotop(s) for s in (d, e, f)]))
         or (all_of([no_co_isotop(s) for s in (a, b, c)])
             and all_of([no_cn_isotop(s) for s in (a, b, c)])
             and all_of([no_co_isotop(s) for s in (d, e, f)])))

def _pred_0c(a, b, c, d, e, f, g):
    return ((_any3(co_isotop, a, b, c) and xor_of([co_isotop(s) for s in (d, e, f)]))
         or (all_of([no_co_isotop(s) for s in (a, b, c)])
             and all_of([no_co_isotop(s) for s in (d, e, f)])))

def _pred_1a(a, b, c, d, e, f, g):
    return ((_any3(cn_isotop, a, b, c) and xor_of([cn_isotop(s) for s in (d, e, f)]))
         or (all_of([no_cn_isotop(s) for s in (a, b, c)])
             and all_of([no_cn_isotop(s) for s in (d, e, f)])))

def _pred_1b(a, b, c, d, e, f, g):
    return ((cn_isotop(b) and xor_of([cn_isotop(s) for s in (d, e, f)]))
         or (no_cn_isotop(b) and any_of([no_cn_isotop(s) for s in (d, e, f)])))

def _pred_15_4chn(a, b, c, d, e, f, g):
    if [n_chngroups(a), n_chngroups(b)] == [n_chngroups(e), n_chngroups(d)]:
        return (chn_isotop(a) == chn_isotop(e) and chn_isotop(b) == chn_isotop(d)
            and co_isotop(a) == co_isotop(e) and co_isotop(b) == co_isotop(d)
            and c18o_bound(a) == c18o_bound(e) and c18o_bound(b) == c18o_bound(d))
    return (implies(chn_isotop(a), chn_isotop(e))
        and co_isotop(a) == co_isotop(e) and co_isotop(b) == co_isotop(d)
        and c18o_bound(a) == c18o_bound(e) and c18o_bound(b) == c18o_bound(d))

def _pred_15_2co(a, b, c, d, e, f, g):
    return (_any3(chn_isotop, a, b, c) == _any3(chn_isotop, d, e, f)
        and _any3(c18o_bound, a, b, c) == _any3(c18o_bound, d, e, f)
        and _any3(co_isotop, a, b, c) == _any3(co_isotop, d, e, f))

def _pred_2(a, b, c, d, e, f, g):
    return neutral(d)[:4] == neutral(b)[:4]

def _pred_2b(a, b, c, d, e, f, g):
    return neutral(e)[:4] == neutral(b)[:4]

def _pred_3(a, b, c, d, e, f, g):
    return neutral(d)[:3] == neutral(a)[:3]

def _pred_3c(a, b, c, d, e, f, g):
    return neutral(e)[:2] == neutral(a)[:2]

def _pred_3d(a, b, c, d, e, f, g):
    return neutral(e)[:2] == neutral(b)[:2]

def _pred_3e(a, b, c, d, e, f, g):
    if is_ion(e) == is_ion(a):
        return neutral(e)[:2] == neutral(a)[:2]
    return neutral(d)[:2] == neutral(a)[:2]

def _pred_4(a, b, c, d, e, f, g):
    return neutral(d)[-3:] == neutral(a)[-3:]

def _pred_4b(a, b, c, d, e, f, g):
    return neutral(d)[-3:] == neutral(b)[-3:]

def _pred_5(a, b, c, d, e, f, g):
    na = neutral(a)
    return _ext_pref_loose(na, neutral(e)) or _ext_pref_loose(na, neutral(d))

def _pred_5b(a, b, c, d, e, f, g):
    nb = neutral(b)
    return _ext_pref_loose(nb, neutral(e)) or _ext_pref_loose(nb, neutral(d))

def _pred_6(a, b, c, d, e, f, g):
    na = neutral(a)
    repl = na.replace("CH3", "CH2")
    def hit(x):
        nx = neutral(x)
        k = min(5, len(na), len(nx))
        return repl[:k] == nx[:k]
    return hit(d) or hit(e) or (neutral(f) != "" and hit(f))

def _pred_6b(a, b, c, d, e, f, g):
    nb = neutral(b)
    repl = nb.replace("CH3", "CH2")
    def hit(x):
        nx = neutral(x)
        k = min(5, len(nb), len(nx))
        return repl[:k] == nx[:k]
    return hit(d) or hit(e) or (neutral(f) != "" and hit(f))

# --------------------------------------------------------------------------- #
# clean-up functions (one per reaction type)
# --------------------------------------------------------------------------- #

def _filtered(ref, rest, pred):
    return [ref] + [r for r in rest if pred(*r[2:9])]

def clean_up_photo(family):
    ref, rest = family[0], family[1:]
    t, ra, rb, rc, rd, re_, rf, rg = ref[1:9]
    ph = t in ("PH", "CP", "CR") and rb in PHOTONS
    if ph and re_ == "e-" and neutral(rd) == neutral(ra):            # pure ionization
        pred = lambda a, b, c, d, e, f, g: d == a + "+" or d == a + "-"
    elif ph and re_ == "e-" and rf == "" and neutral(rd) != neutral(ra):
        pred = lambda a, b, c, d, e, f, g: (                          # restructuring
            sticky_13c([a, b, c, d, e, f, g]) and sticky_18o([a, b, c, d, e, f, g])
            and _when(co_bound(a) == co_bound(d) and chn_bound(a) == chn_bound(d)
                      and n_cogroups(a) == n_cogroups(d) == 1,
                      chn_isotop(a) == chn_isotop(d) and co_isotop(a) == co_isotop(d))
            and _when(cn_bound(a) == cn_bound(d) and chn_bound(a) == chn_bound(d),
                      chn_isotop(a) == chn_isotop(d) and cn_isotop(a) == cn_isotop(d))
            and _when(cn_bound(a) == cn_bound(d) and chn_bound(a) != chn_bound(d),
                      cn_isotop(a) == cn_isotop(d)))
    elif ph and ra in _perm_joins([rd, re_, rf]):                    # 3-part dissociation
        pred = lambda a, b, c, d, e, f, g: any(
            a == "".join(p) for p in permutations([d, e, f]))
    elif ph and (ra == rd + re_ or ra == re_ + rd):                  # 2-part dissociation
        pred = lambda a, b, c, d, e, f, g: d + e == a or e + d == a
    elif (ph and chn_bound(ra) and co_bound(ra)
          and count_catoms(ra) >= n_cogroups(ra) + n_chngroups(ra)):
        pred = lambda a, b, c, d, e, f, g: (                          # CO + CHn
            sticky_18o([a, b, c, d, e, f, g])
            and _when(co_bound(d) and chn_bound(e),
                      co_isotop(d) == co_isotop(a) and chn_isotop(e) == chn_isotop(a))
            and _when(co_bound(e) and chn_bound(d),
                      co_isotop(e) == co_isotop(a) and chn_isotop(d) == chn_isotop(a)))
    elif ph and n_chngroups(ra) == n_chngroups(rd) == 1:
        pred = lambda a, b, c, d, e, f, g: (sticky_18o([a, b, c, d, e, f, g])
                                            and chn_isotop(a) == chn_isotop(d))
    elif ph and n_chngroups(ra) == n_chngroups(re_) == 1:
        pred = lambda a, b, c, d, e, f, g: (sticky_18o([a, b, c, d, e, f, g])
                                            and chn_isotop(a) == chn_isotop(e))
    elif ph and cn_bound(ra) and cn_bound(re_) and not cn_bound(rd):
        pred = lambda a, b, c, d, e, f, g: (sticky_18o([a, b, c, d, e, f, g])
                                            and cn_isotop(a) == cn_isotop(e))
    else:
        return list(family)                                          # no match
    return _filtered(ref, rest, pred)

def clean_up_ra(family):
    ref, rest = family[0], family[1:]
    t, ra, rb, rc, rd, re_, rf, rg = ref[1:9]
    if t == "RA" and is_ion(rd) and rd in _charged_perm_joins([ra, rb, rc]):
        pred = lambda a, b, c, d, e, f, g: any(
            d == j for j in _charged_perm_joins([a, b, c]))
    elif t == "RA" and co_bound(rb) and co_bound(rd):                # keep CO binding
        pred = lambda a, b, c, d, e, f, g: co_isotop(d) == co_isotop(b)
    else:
        return list(family)
    return _filtered(ref, rest, pred)

def clean_up_ad(family):
    ref, rest = family[0], family[1:]
    t, ra, rb, rc, rd, re_, rf, rg = ref[1:9]
    if t == "AD" and re_ == "e-" and rd in _neutral_perm_joins([ra, rb, rc]):
        pred = lambda a, b, c, d, e, f, g: any(
            d == j for j in _neutral_perm_joins([a, b, c]))
    else:
        return list(family)
    return _filtered(ref, rest, pred)

def clean_up_dr(family):
    ref, rest = family[0], family[1:]
    t, ra, rb, rc, rd, re_, rf, rg = ref[1:9]
    d_ad_ref = delta_ad(ra, rb, rd, re_)
    d_ax_ref = delta_ax(ra, rb, rd, re_)
    chn_r = _any3(chn_bound, ra, rb, rc)
    chn_p = _any3(chn_bound, rd, re_, rf)
    nc_r = sum(count_catoms(s) for s in (ra, rb, rc))
    nch_r = sum(n_chngroups(s) for s in (ra, rb, rc))
    nch_p = sum(n_chngroups(s) for s in (rd, re_, rf))
    nco_r = sum(n_cogroups(s) for s in (ra, rb, rc))
    nco_p = sum(n_cogroups(s) for s in (rd, re_, rf))

    if (t == "DR" and chn_r and chn_p and nc_r > 1 and nch_r == nch_p
            and nco_r == nco_p == 2
            and not (co_bound(rd) and co_bound(re_))):               # case 1.5, 2CO
        def pred(a, b, c, d, e, f, g):
            if n_cogroups(d) == 2 and n_chngroups(d) == 0 and n_chngroups(e) == 1:
                return (implies(co_isotop(d), co_isotop(a) and not chn_isotop(a))
                    and implies(chn_isotop(e), chn_isotop(a)))
            return _pred_15_2co(a, b, c, d, e, f, g)
    elif (t == "DR" and chn_r and chn_p and nc_r > 1
          and nch_r == nch_p < 3 and nco_r == nco_p <= 1):           # case 1.5
        def pred(a, b, c, d, e, f, g):
            ok = True
            if _match_delta(d_ad_ref, {2: 0, 4: -1}):
                ok &= [has_13c(a), has_13c(b)] == [has_13c(d), has_13c(e)]
            if _match_delta(d_ad_ref, {2: 0, 4: 1}):
                ok &= [has_13c(a), has_13c(b)] == [has_13c(e), has_13c(d)]
            if d_ax_ref == [0, 0, 1, 0, 0, 0, 0, 0, 0]:
                ok &= implies(has_13c(a), has_13c(d))
            if not _match_delta([abs(x) for x in d_ad_ref], {4: 1}):
                ok &= _pred_15_2co(a, b, c, d, e, f, g)
            return bool(ok)
    elif (t == "DR" and rb == "e-" and is_ion(ra)
          and neutral(ra) in _perm_joins([rd, re_, rf])):            # case 1
        pred = lambda a, b, c, d, e, f, g: any(
            neutral(a) == "".join(p) for p in permutations([d, e, f]))
    elif (t == "DR" and rb == "e-" and is_ion(ra)
          and (neutral(ra) == rd + re_ or neutral(ra) == re_ + rd)):  # case 2
        pred = lambda a, b, c, d, e, f, g: any(
            neutral(a) == "".join(p) for p in permutations([d, e, f]))
    elif (t == "DR" and rb == "e-" and is_ion(ra) and len(rd) >= 4 and len(ra) >= 4
          and rd[:4] == ra[:4]):                                     # case 3
        pred = lambda a, b, c, d, e, f, g: d[:4] == a[:4]
    elif (t == "DR" and rb == "e-" and is_ion(ra) and len(rd) >= 3 and len(ra) >= 3
          and rd[-3:] == neutral(ra)[-3:]):                          # case 4
        pred = lambda a, b, c, d, e, f, g: d[-3:] == neutral(a)[-3:]
    elif (t == "DR" and rb == "e-" and is_ion(ra)
          and (_ext_pref(neutral(ra), rd) or _ext_pref(neutral(ra), re_))):  # case 5
        pred = lambda a, b, c, d, e, f, g: (
            _ext_pref_loose(a[:-1], e) or _ext_pref_loose(a[:-1], d))
    elif (t == "DR" and rb == "e-" and is_ion(ra) and len(ra) > 3
          and ra[:3] == "CH3"
          and ((len(rd) > 2 and _ch3_shrink_match(ra, rd))
               or (len(re_) > 2 and _ch3_shrink_match(ra, re_))
               or (len(rf) > 2 and _ch3_shrink_match(ra, rf)))):     # case 6
        def pred(a, b, c, d, e, f, g):
            repl = a.replace("CH3", "CH2")
            def m(x):
                k = min(5, len(x))
                return repl[:k] == x[:k]
            return m(d) or m(e) or (f != "" and m(f))
    elif (t == "DR" and rb == "e-" and is_ion(ra) and len(ra) > 3
          and ra[:3] == "CH3"
          and ((len(rd) > 2 and _ch3_shrink_match(ra, rd))
               or (len(re_) > 2 and _ch3_shrink_match(ra, re_))
               or (len(rf) > 2 and _ch3_shrink_match(ra, rf)))):
        # NOTE: the original's case 7 (CH3->CH4) shares the ref condition of
        # case 6 in this transcription order and is therefore unreachable
        # after case 6 -- kept separate here for auditability:
        def pred(a, b, c, d, e, f, g):                               # CH3 -> CH4
            repl = a.replace("CH3", "CH4")
            def m(x):
                k = min(4, len(x))
                return repl[:k] == x[:k]
            return m(d) or m(e) or (f != "" and m(f))
    elif (t == "DR" and rb == "e-" and is_ion(ra) and len(ra) > 3
          and ra[:3] == "CH2"
          and any(len(x) > 0 for x in (rd, re_, rf))):               # case 8
        def pred(a, b, c, d, e, f, g):
            repl = a.replace("CH2", "CH")
            def m(x):
                k = min(5, len(x))
                return repl[:k] == x[:k]
            return m(d) or m(e) or (f != "" and m(f))
    else:
        return list(family)
    return _filtered(ref, rest, pred)

def _ce_mn_case12_pred(a, b, c, d, e, f, g, comp_a_d):
    """shared body of cleanUpCE/cleanUpMN 'restructuring' case"""
    ok = sticky_18o([a, b, c, d, e, f, g]) and sticky_13c([a, b, c, d, e, f, g])
    if comp_a_d:                                # a<->d and b<->e
        ok &= _when(co_bound(b) == co_bound(e), co_isotop(b) == co_isotop(e))
        ok &= _when(co_bound(a) == co_bound(d), co_isotop(a) == co_isotop(d))
        ok &= _when(chn_bound(b) == chn_bound(e), chn_isotop(b) == chn_isotop(e))
        ok &= _when(chn_bound(a) == chn_bound(d), chn_isotop(a) == chn_isotop(d))
        ok &= has_13c(a) == has_13c(d) and has_13c(b) == has_13c(e)
        ok &= has_18o(a) == has_18o(d) and has_18o(b) == has_18o(e)
    else:                                       # a<->e and b<->d
        ok &= _when(co_bound(b) == co_bound(d), co_isotop(b) == co_isotop(d))
        ok &= _when(co_bound(a) == co_bound(e), co_isotop(a) == co_isotop(e))
        ok &= _when(chn_bound(b) == chn_bound(d), chn_isotop(b) == chn_isotop(d))
        ok &= _when(chn_bound(a) == chn_bound(e), chn_isotop(a) == chn_isotop(e))
        ok &= has_13c(a) == has_13c(e) and has_13c(b) == has_13c(d)
        ok &= has_18o(a) == has_18o(e) and has_18o(b) == has_18o(d)
    return bool(ok)

def clean_up_ce(family):
    ref, rest = family[0], family[1:]
    t, ra, rb, rc, rd, re_, rf, rg = ref[1:9]
    if (t == "CE" and rf == ""
            and ((neutral(ra) == neutral(re_) and neutral(rb) == neutral(rd))
                 or (neutral(ra) == neutral(rd) and neutral(rb) == neutral(re_)))):
        pred = lambda a, b, c, d, e, f, g: (                     # places are kept
            (neutral(a) == neutral(d) or neutral(a) == neutral(e))
            and (neutral(b) == neutral(d) or neutral(b) == neutral(e)))
    elif t == "CE" and rf == "" and (hcnos(ra) == hcnos(rd) or hcnos(ra) == hcnos(re_)):
        comp = hcnos(ra) == hcnos(rd)
        pred = lambda a, b, c, d, e, f, g: _ce_mn_case12_pred(
            a, b, c, d, e, f, g, comp)
    elif t == "CE" and rf == "" and (hcnos(ra) != hcnos(rd) or hcnos(ra) != hcnos(re_)):
        pred = lambda a, b, c, d, e, f, g: (sticky_18o([a, b, c, d, e, f, g])
                                            and sticky_13c([a, b, c, d, e, f, g]))
    else:
        return list(family)
    return _filtered(ref, rest, pred)

def clean_up_mn(family):
    ref, rest = family[0], family[1:]
    t, ra, rb, rc, rd, re_, rf, rg = ref[1:9]
    if (t == "MN" and rf == ""
            and ((neutral(ra) == neutral(re_) and neutral(rb) == neutral(rd))
                 or (neutral(ra) == neutral(rd) and neutral(rb) == neutral(re_)))):
        pred = lambda a, b, c, d, e, f, g: (
            (neutral(a) == neutral(d) or neutral(a) == neutral(e))
            and (neutral(b) == neutral(d) or neutral(b) == neutral(e)))
    elif t == "MN" and rf == "" and (hcnos(ra) == hcnos(rd) or hcnos(ra) == hcnos(re_)):
        comp = hcnos(ra) == hcnos(rd)
        pred = lambda a, b, c, d, e, f, g: _ce_mn_case12_pred(a, b, c, d, e, f, g, comp)
    elif t == "MN" and rf == "" and (hcnos(ra) != hcnos(rd) or hcnos(ra) != hcnos(re_)):
        pred = lambda a, b, c, d, e, f, g: (sticky_18o([a, b, c, d, e, f, g])
                                            and sticky_13c([a, b, c, d, e, f, g]))
    elif t == "MN" and rf != "":               # dissociative-recombination-like
        def pred(a, b, c, d, e, f, g):
            ha, hb = hcnos(a), hcnos(b)
            hd, he, hf = hcnos(d), hcnos(e), hcnos(f)
            return ((neutral(a) == neutral(d) and hb == _vadd(he, hf))
                 or (neutral(b) == neutral(d) and ha == _vadd(he, hf))
                 or (neutral(a) == neutral(e) and hb == _vadd(hd, hf))
                 or (neutral(a) == neutral(f) and hb == _vadd(hd, he))
                 or (neutral(b) == neutral(e) and ha == _vadd(hd, hf))
                 or (neutral(b) == neutral(f) and ha == _vadd(hd, he)))
    else:
        return list(family)
    return _filtered(ref, rest, pred)

def clean_up_cl(family):
    ref, rest = family[0], family[1:]
    if ref[1] != "CL":
        return list(family)
    pred = lambda a, b, c, d, e, f, g: (
        sticky_13c([a, b, c, d, e, f, g]) and sticky_18o([a, b, c, d, e, f, g])
        and _c18o_together(a, b, c, d, e, f))
    return _filtered(ref, rest, pred)

# --------------------------------------------------------------------------- #
# cleanUpIN and cleanUpNN -- the two large heuristic blocks
# --------------------------------------------------------------------------- #

def _common_cases(kind, ref, rest, ra, rb, rc, rd, re_, rf, rg):
    """Cases shared (condition text and predicate) between IN and NN.
    Returns (handled, result)."""
    t = ref[1]
    d_ad_ref = delta_ad(ra, rb, rd, re_)
    d_ae_ref = delta_ae(ra, rb, rd, re_)
    nc_r = sum(count_catoms(s) for s in (ra, rb, rc))
    no_r = sum(count_oatoms(s) for s in (ra, rb, rc))

    if nc_r <= 1 and no_r <= 1:                                   # trivial
        return True, [ref] + list(rest)
    if (rf == "" and not mixed_signs(d_ad_ref) and not mixed_signs(d_ae_ref)
            and d_ad_ref == [-x for x in d_ae_ref]):              # arbitrary
        def pred(a, b, c, d, e, f, g):
            if kind == "NN":
                if [abs(x) for x in delta_ad(a, b, d, e)] == [1, 0, 0, 0, 0, 0, 0, 0, 0]:
                    return _c18o_together(a, b, c, d, e, f)
                return True
            # IN: the two 8-slot Sign@Abs patterns never match (length check)
            if (all(co_bound(s) for s in (ra, rb, rd, re_))
                    and not _match_delta([_sgn(abs(x)) for x in d_ad_ref], {2: 1}, length=8)
                    and not _match_delta([_sgn(abs(x)) for x in d_ad_ref], {4: 1})):
                return (sorted(mapthread_and([has_13c(a), has_13c(b)],
                                             [c18o_bound(a), c18o_bound(b)]))
                    == sorted(mapthread_and([has_13c(d), has_13c(e)],
                                            [c18o_bound(d), c18o_bound(e)])))
            return True
        return True, _filtered(ref, rest, pred)
    ae_abs = [abs(x) for x in d_ae_ref]
    ad_abs = [abs(x) for x in d_ad_ref]
    unit = [1, 0, 0, 0, 0, 0, 0, 0, 0]
    if (rf == "" and nc_r > 1 and ae_abs != ad_abs
            and (ae_abs == unit or ad_abs == unit)):              # H/proton transfer
        pred = lambda a, b, c, d, e, f, g: _proton_transfer_pred(
            a, b, c, d, e, f, g, ae_abs == unit)
        return True, _filtered(ref, rest, pred)
    if (cn_bound(ra) and co_bound(ra)) or (cn_bound(rb) and co_bound(rb)):
        if kind == "IN":                                          # case 0a
            pred = lambda a, b, c, d, e, f, g: _pred_1a(a, b, c, d, e, f, g)
        else:                                                     # NN 0a (extra Ifs)
            def pred(a, b, c, d, e, f, g):
                base = _pred_1a(a, b, c, d, e, f, g)
                dax = delta_ax(a, b, d, e)
                fpA = find_partner_a(a, b, d, e)
                allC = sum(_sgn(count_catoms(s)) for s in (a, b, d, e)) == 4
                allO = sum(_sgn(count_oatoms(s)) for s in (a, b, d, e)) == 4
                ok = base
                ok &= _when(allC, has_13c(b) == has_13c(e) and has_13c(a) == has_13c(d))
                ok &= _when(allO and fpA == 0,
                            implies(has_18o(a), has_18o(d))
                            and implies(has_18o(b), has_18o(e)))
                ok &= _when(allO and fpA == 1 and [abs(x) for x in dax][4] == 0,
                            implies(has_18o(a), has_18o(d))
                            and implies(has_18o(b), has_18o(e)))
                ok &= _when(allO and fpA == 2 and [abs(x) for x in dax][4] == 0,
                            implies(has_18o(a), has_18o(e))
                            and implies(has_18o(b), has_18o(d)))
                return bool(ok)
        return True, _filtered(ref, rest, pred)
    if (nc_r >= 3 and no_r < 2
            and ((cn_bound(ra) and co_bound(rb)) or (cn_bound(rb) and co_bound(ra)))
            and xor_of([co_bound(s) for s in (rd, re_, rf)])
            and xor_of([cn_bound(s) for s in (rd, re_, rf)])):    # case 0b
        return True, _filtered(ref, rest, _pred_0b)
    if (nc_r >= 3 and no_r > 1
            and ((cn_bound(ra) and co_bound(rb)) or (cn_bound(rb) and co_bound(ra)))
            and xor_of([co_bound(s) for s in (rd, re_, rf)])
            and xor_of([cn_bound(s) for s in (rd, re_, rf)])):    # case 0b1
        return True, _filtered(ref, rest, _pred_0b1)
    if (nc_r >= 3
            and ((cn_bound(ra) and co_bound(rb)) or (cn_bound(rb) and co_bound(ra)))
            and xor_of([co_bound(s) for s in (rd, re_, rf)])):    # case 0b2
        return True, _filtered(ref, rest, _pred_0b2)
    if (nc_r < 3
            and ((cn_bound(ra) and co_bound(rb)) or (cn_bound(rb) and co_bound(ra)))
            and xor_of([co_bound(s) for s in (rd, re_, rf)])):    # case 0c
        return True, _filtered(ref, rest, _pred_0c)
    if (_any3(cn_bound, ra, rb, rc) and _any3(cn_bound, rd, re_, rf)
            and nc_r > 1 and sum(count_natoms(s) for s in (ra, rb, rc)) == 1):
        return True, _filtered(ref, rest, _pred_1a)               # case 1a
    if (cn_bound(rb) and _any3(cn_bound, rd, re_, rf)
            and nc_r > 1 and count_natoms(rb) == 1):
        return True, _filtered(ref, rest, _pred_1b)               # case 1b
    return False, None

def _string_overlap_cases(ref, rest, ra, rb, rc, rd, re_, rf, rg):
    """cases 2..6b shared by IN and NN (all based on neutralised strings)"""
    nrd, nre, nra, nrb, nrf = (neutral(x) for x in (rd, re_, ra, rb, rf))
    co_r = _any3(co_bound, ra, rb, rc)
    co_p = _any3(co_bound, rd, re_, rf)
    if (len(nrd) >= 5 and len(nrb) >= 5 and len(nre) >= 5 and len(nra) >= 5
            and nrd[:5] == nrb[:5] and nre[:5] == nra[:5]):       # case 2-0
        pred = lambda a, b, c, d, e, f, g: (
            neutral(d)[:5] == neutral(b)[:5] and neutral(e)[:5] == neutral(a)[:5])
    elif (len(nrd) >= 4 and len(nrb) >= 4
          and (nrd[:4] == nrb[:4] or nrd[-3:] == nrb[-3:])
          and co_r and co_p):                                     # case 2-CO
        pred = _co_variant(lambda a, b, c, d, e, f, g:
                           neutral(d)[:4] == neutral(b)[:4]
                           or neutral(d)[-3:] == neutral(b)[-3:])
    elif len(nrd) >= 4 and len(nrb) >= 4 and nrd[:4] == nrb[:4]:  # case 2
        pred = _pred_2
    elif (len(nre) >= 4 and len(nrb) >= 4
          and (nre[:4] == nrb[:4] or nre[-3:] == nrb[-3:])
          and co_r and co_p):                                     # case 2b-CO
        pred = _co_variant(lambda a, b, c, d, e, f, g:
                           neutral(e)[:4] == neutral(b)[:4]
                           or neutral(e)[-3:] == neutral(b)[-3:])
    elif len(nre) >= 4 and len(nrb) >= 4 and nre[:4] == nrb[:4]:  # case 2b
        pred = _pred_2b
    elif (len(nrd) >= 4 and len(nra) >= 4 and nrd[:4] == nra[:4]
          and len(nre) >= 2 and nrd[:2] != nre[:2]):              # case 3
        pred = _pred_3
    elif len(nrd) >= 4 and len(nrb) >= 4 and nrd[:4] == nrb[:4]:  # case 3b
        pred = _pred_2
    elif (len(nre) >= 3 and len(nra) >= 3 and len(nrd) >= 3
          and nre[:3] == nra[:3] and nre[:2] != nrd[:2]):         # case 3c
        pred = _pred_3c
    elif (len(nre) >= 3 and len(nrb) >= 3 and len(nra) >= 2
          and nre[:3] == nrb[:3] and nra[:2] != nrb[:2]):         # case 3d
        pred = _pred_3d
    elif (len(nre) >= 3 and len(nra) >= 3 and len(nrd) >= 2
          and nre[:3] == nra[:3] and nre[:2] == nrd[:2]
          and sum(count_catoms(s) for s in (ra, rb, rc)) > 1):    # case 3e
        pred = _pred_3e
    elif (len(nrd) >= 3 and len(nra) >= 3 and nrd[-3:] == nra[-3:]
          and co_r and co_p):                                     # case 4-CO
        pred = _co_variant(_pred_4)
    elif len(nrd) >= 3 and len(nra) >= 3 and nrd[-3:] == nra[-3:]:
        pred = _pred_4                                            # case 4
    elif (len(nrd) >= 3 and len(nrb) >= 3 and nrd[-3:] == nrb[-3:]
          and co_r and co_p):                                     # case 4b-CO
        pred = _co_variant(_pred_4b)
    elif len(nrd) >= 3 and len(nrb) >= 3 and nrd[-3:] == nrb[-3:]:
        pred = _pred_4b                                           # case 4b
    elif (_ext_pref(nra, nre) or _ext_pref(nra, nrd)):            # case 5
        pred = _pred_5
    elif (_ext_pref(nrb, nre) or _ext_pref(nrb, nrd)):            # case 5b
        pred = _pred_5b
    elif (len(nra) >= 3 and nra[:3] == "CH3"
          and ((len(nrd) > 2 and _ch3_shrink_match(nra, nrd))
               or (len(nre) > 2 and _ch3_shrink_match(nra, nre))
               or (len(nrf) > 2 and _ch3_shrink_match(nra, nrf)))):  # case 6
        pred = _pred_6
    elif (len(nrb) >= 3 and nrb[:3] == "CH3"
          and ((len(nrd) > 2 and _ch3_shrink_match(nrb, nrd))
               or (len(nre) > 2 and _ch3_shrink_match(nrb, nre))
               or (len(nrf) > 2 and _ch3_shrink_match(nrb, nrf)))):  # case 6b
        pred = _pred_6b
    else:
        return False, None
    return True, _filtered(ref, rest, pred)

def clean_up_in(family):
    ref, rest = family[0], family[1:]
    t, ra, rb, rc, rd, re_, rf, rg = ref[1:9]
    d_ad_ref = delta_ad(ra, rb, rd, re_)
    d_ae_ref = delta_ae(ra, rb, rd, re_)
    d_bd_ref = delta_bd(ra, rb, rd, re_)
    d_be_ref = delta_be(ra, rb, rd, re_)
    d_bx_ref = delta_bx(ra, rb, rd, re_)
    d_ax_ref = delta_ax(ra, rb, rd, re_)
    fpA_ref = find_partner_a(ra, rb, rd, re_)
    fpB_ref = find_partner_b(ra, rb, rd, re_)
    nc_r = sum(count_catoms(s) for s in (ra, rb, rc))

    if ref[1:9] == ["IN", "H3+", "CH3COCH3", "", "C2H3+", "CH4", "H2O", ""]:
        return _filtered(ref, rest,
                         lambda a, b, c, d, e, f, g: implies(co_isotop(b), has_13c(d)))
    if ref[1:9] == ["IN", "H3+", "CH3COCH3", "", "CH3+", "C2H5OH", "", ""]:
        return _filtered(ref, rest,
                         lambda a, b, c, d, e, f, g: implies(co_isotop(b), has_13c(e)))
    if ref[1:9] == ["IN", "CH3OH", "CH3OH2+", "", "CH3OCH4+", "H2O", "", ""]:
        return _filtered(ref, rest,
                         lambda a, b, c, d, e, f, g:
                         has_18o(a) == has_18o(d) and has_18o(b) == has_18o(e))

    handled, result = _common_cases("IN", ref, rest, ra, rb, rc, rd, re_, rf, rg)
    if handled:
        return result

    chn_r = _any3(chn_bound, ra, rb, rc)
    chn_p = _any3(chn_bound, rd, re_, rf)
    nch_r = sum(n_chngroups(s) for s in (ra, rb, rc))
    nch_p = sum(n_chngroups(s) for s in (rd, re_, rf))
    nco_r = sum(n_cogroups(s) for s in (ra, rb, rc))
    nco_p = sum(n_cogroups(s) for s in (rd, re_, rf))

    if t == "IN" and all_of([chn_bound(s) for s in (ra, rb, rd, re_)]):   # 1.5 4CHn
        return _filtered(ref, rest, _pred_15_4chn)
    if (t == "IN" and chn_r and chn_p and nc_r > 1 and nch_r == nch_p
            and nco_r == nco_p == 2
            and not (co_bound(rd) and co_bound(re_))):            # 1.5 2CO
        return _filtered(ref, rest, _pred_15_2co)
    if (t == "IN" and chn_r and chn_p and nc_r > 1
            and nch_r == nch_p < 3 and nco_r == nco_p <= 1):      # 1.5
        def pred(a, b, c, d, e, f, g):
            ok = True
            if _match_delta(d_ad_ref, {2: 0, 4: -1}):
                ok &= [has_13c(a), has_13c(b)] == [has_13c(d), has_13c(e)]
            if _match_delta(d_ad_ref, {2: 0, 4: 1}):
                ok &= [has_13c(a), has_13c(b)] == [has_13c(e), has_13c(d)]
            if d_ax_ref == [0, 0, 1, 0, 0, 0, 0, 0, 0]:
                ok &= implies(has_13c(a), has_13c(d))
            if f != "":
                ok &= _when(co_bound(b) and co_bound(d), co_isotop(d) == co_isotop(b))
            if not _match_delta([abs(x) for x in d_ad_ref], {4: 1}):
                ok &= _pred_15_2co(a, b, c, d, e, f, g)
            return bool(ok)
        return _filtered(ref, rest, pred)
    if (t == "IN" and chn_r and chn_p and nc_r > 1 and nch_r == nch_p
            and (nco_r == nco_p == 3 or nco_r != nco_p)
            and fpA_ref != 0 and fpB_ref != 0):                   # 1.5 3/!=CO
        def pred(a, b, c, d, e, f, g):
            if fpA_ref == 1:
                if d_ax_ref[0] != 0:
                    return (chn_isotop(a) == chn_isotop(d)
                        and co_isotop(a) == co_isotop(d)
                        and has_18o(a) == has_18o(d)
                        and left18o_bound(a) == left18o_bound(d)
                        and left18o_bound(b) == left18o_bound(e)
                        and cn_isotop(a) == cn_isotop(d))
                if d_ax_ref[4] != 0:
                    return (chn_isotop(a) == chn_isotop(d)
                        and cn_isotop(a) == cn_isotop(d)
                        and left18o_bound(a) == left18o_bound(d)
                        and left18o_bound(b) == left18o_bound(e))
                return True
            else:
                if d_ax_ref[0] != 0:
                    return (chn_isotop(a) == chn_isotop(e)
                        and co_isotop(a) == co_isotop(e)
                        and has_18o(a) == has_18o(e)
                        and cn_isotop(a) == cn_isotop(e)
                        and left18o_bound(a) == left18o_bound(e)
                        and left18o_bound(b) == left18o_bound(d))
                if d_ax_ref[4] != 0:
                    return (chn_isotop(a) == chn_isotop(e)
                        and cn_isotop(a) == cn_isotop(e)
                        and left18o_bound(a) == left18o_bound(e)
                        and left18o_bound(b) == left18o_bound(d))
                return True
        return _filtered(ref, rest, pred)
    if (t == "IN" and chn_r and _any3(chn_bound, rd, re_, rf, rg) and nc_r > 1
            and nch_r == sum(n_chngroups(s) for s in (rd, re_, rf, rg)) < 3
            and nco_r > sum(n_cogroups(s) for s in (rd, re_, rf, rg)) <= 1
            and rg != ""):                                        # 1.5 4 products
        def pred(a, b, c, d, e, f, g):
            return (any_of([chn_isotop(s) for s in (a, b, c)])
                == any_of([chn_isotop(s) for s in (d, e, f, g)])
                and any_of([c18o_bound(s) for s in (a, b, c)])
                == any_of([c18o_bound(s) for s in (d, e, f, g)]))
        return _filtered(ref, rest, pred)
    if (t == "IN" and chn_r and chn_p and nc_r > 1 and rf == ""
            and _any3(co_bound, ra, rb, rc) == _any3(co_bound, rd, re_, rf)
            and nch_r != nch_p and fpA_ref != 0):                 # case 1.6
        def pred(a, b, c, d, e, f, g):
            dax = delta_ax(a, b, d, e)
            axabs = [abs(x) for x in dax]
            fpA = find_partner_a(a, b, d, e)
            ok = True
            if dax[2] == 0 and axabs != [3, 0, 0, 0, 0, 0, 0, 0, 0] \
                    and axabs != [2, 0, 0, 0, 0, 0, 0, 0, 0]:
                ok &= ([has_13c(a), has_13c(b)] == [has_13c(d), has_13c(e)] if fpA == 1
                       else [has_13c(a), has_13c(b)] == [has_13c(e), has_13c(d)])
            if dax in ([-3, 0, 0, 0, 0, 0, 0, 0, 0], [-2, 0, 0, 0, 0, 0, 0, 0, 0]) \
                    and n_chngroups(rb) == 2:
                ok &= _when(chn_bound(b) == chn_bound(e),
                            implies(chn_isotop(e), chn_isotop(b)))
            if dax[4] == 0:
                ok &= ([has_18o(a), has_18o(b)] == [has_18o(d), has_18o(e)] if fpA == 1
                       else [has_18o(a), has_18o(b)] == [has_18o(e), has_18o(d)])
            if fpA == 1 and axabs != [3, 0, 0, 0, 0, 0, 0, 0, 0] \
                    and axabs != [2, 0, 0, 0, 0, 0, 0, 0, 0]:
                dbx = delta_bx(a, b, d, e)
                h1, h3 = dax[0], dax[2]
                a1, a3 = axabs[0], axabs[2]
                ok &= _when([h1, h3] == [-2, -1] and chn_bound(a) == chn_bound(d),
                            implies(chn_isotop(d), chn_isotop(a)))
                ok &= _when([h1, h3] == [-2, -1] and chn_bound(b) == chn_bound(e),
                            implies(chn_isotop(e), chn_isotop(b)))
                ok &= _when([a1, a3] == [3, 1] and chn_bound(a) == chn_bound(d),
                            implies(chn_isotop(d), chn_isotop(a)))
                ok &= _when([a1, a3] == [3, 1] and chn_bound(b) == chn_bound(e),
                            implies(chn_isotop(e), chn_isotop(b)))
                ok &= _when([h1, h3] == [-3, -1] and co_bound(a) == co_bound(d),
                            co_isotop(d) == co_isotop(a))
                ok &= _when([dbx[0], dbx[2]] == [3, 1] and co_bound(b) == co_bound(e),
                            co_isotop(b) == co_isotop(e))
                ok &= _when([h1, h3] == [3, 1] and co_bound(a) == co_bound(d),
                            co_isotop(d) == co_isotop(a))
                ok &= _when([dbx[0], dbx[2]] == [-3, -1] and co_bound(b) == co_bound(e),
                            co_isotop(b) == co_isotop(e))
                ok &= _when([a1, a3] != [3, 1] and [a1, a3] != [2, 1]
                            and chn_bound(a) == chn_bound(d),
                            chn_isotop(a) == chn_isotop(d))
                ok &= _when([a1, a3] != [3, 1] and [a1, a3] != [2, 1]
                            and chn_bound(b) == chn_bound(e),
                            chn_isotop(b) == chn_isotop(e))
            elif fpA == 2 and axabs != [3, 0, 0, 0, 0, 0, 0, 0, 0] \
                    and axabs != [2, 0, 0, 0, 0, 0, 0, 0, 0]:
                dbx = delta_bx(a, b, d, e)
                h1, h3 = dax[0], dax[2]
                a1, a3 = axabs[0], axabs[2]
                ok &= _when([h1, h3] == [-2, -1] and chn_bound(a) == chn_bound(e),
                            implies(chn_isotop(e), chn_isotop(a)))
                ok &= _when([h1, h3] == [-2, -1] and chn_bound(b) == chn_bound(d),
                            implies(chn_isotop(d), chn_isotop(b)))
                ok &= _when([a1, a3] != [3, 1] and [a1, a3] != [2, 1]
                            and chn_bound(a) == chn_bound(e),
                            chn_isotop(a) == chn_isotop(e))
                ok &= _when([a1, a3] != [3, 1] and [a1, a3] != [2, 1]
                            and chn_bound(b) == chn_bound(d),
                            chn_isotop(b) == chn_isotop(d))
                ok &= _when([a1, a3] == [3, 1] and chn_bound(a) == chn_bound(e),
                            implies(chn_isotop(e), chn_isotop(a)))
                ok &= _when([a1, a3] == [3, 1] and chn_bound(b) == chn_bound(d),
                            implies(chn_isotop(d), chn_isotop(b)))
                ok &= _when([h1, h3] == [-3, -1] and co_bound(a) == co_bound(e),
                            co_isotop(e) == co_isotop(a))
                ok &= _when([dbx[0], dbx[2]] == [3, 1] and co_bound(b) == co_bound(d),
                            co_isotop(b) == co_isotop(d))
                ok &= _when([h1, h3] == [3, 1] and co_bound(a) == co_bound(e),
                            co_isotop(e) == co_isotop(a))
                ok &= _when([dbx[0], dbx[2]] == [-3, -1] and co_bound(b) == co_bound(d),
                            co_isotop(b) == co_isotop(d))
            return bool(ok)
        return _filtered(ref, rest, pred)

    handled, result = _string_overlap_cases(ref, rest, ra, rb, rc, rd, re_, rf, rg)
    if handled:
        return result

    # ---- remaining IN-only cases ---------------------------------------- #
    co_r = _any3(co_bound, ra, rb, rc)
    co_p = _any3(co_bound, rd, re_, rf)
    if co_r and co_p and nco_r == nco_p:                          # X iso-CO
        def pred(a, b, c, d, e, f, g):
            dax = delta_ax(a, b, d, e); axabs = [abs(x) for x in dax]
            dae = delta_ae(a, b, d, e); daeabs = [abs(x) for x in dae]
            dad = delta_ad(a, b, d, e); dadabs = [abs(x) for x in dad]
            dbx = delta_bx(a, b, d, e); dbxabs = [abs(x) for x in dbx]
            dbe = delta_be(a, b, d, e); dbeabs = [abs(x) for x in dbe]
            fpA = find_partner_a(a, b, d, e)
            unit = [1, 0, 0, 0, 0, 0, 0, 0, 0]
            coga = sorted([n_cogroups(a), n_cogroups(b)])
            cogd = sorted([n_cogroups(d), n_cogroups(e)])
            ca, cb, cd, ce = (count_catoms(s) for s in (a, b, d, e))
            oa, ob, od, oe = (count_oatoms(s) for s in (a, b, d, e))
            ok = True
            if f == "" and coga == cogd and _match_delta(dbxabs, {4: 0}):
                ok &= _c18o_together(a, b, c, d, e, f)
            if f == "" and coga != cogd and fpA == 1 and _match_delta(dbeabs, {4: 1}):
                ok &= (implies(has_18o(e), has_18o(b))
                    and implies(has_13c(e), has_13c(b))
                    and implies(has_13c(d), has_13c(a)))
            if (f == "" and coga == cogd and fpA == 0
                    and daeabs != dadabs and daeabs == unit):
                ok &= (has_18o(a) == has_18o(e) and has_18o(b) == has_18o(d)
                    and implies(has_13c(e), has_13c(a))
                    and implies(has_13c(d), has_13c(b)))
            if (f == "" and coga == cogd and fpA == 0
                    and daeabs == dadabs and daeabs == unit):
                if is_ion(a) == is_ion(e):
                    ok &= (_when(ca > 0 and ce > 0, has_13c(a) == has_13c(e))
                        and _when(oa > 0 and oe > 0, has_18o(a) == has_18o(e))
                        and _when(cd == 1 and cb > 0, has_13c(b) == has_13c(d))
                        and _when(od > 0 and ob > 0, has_18o(b) == has_18o(d)))
                else:
                    ok &= (_when(cb > 0 and ce > 0, has_13c(b) == has_13c(e))
                        and _when(ob > 0 and oe > 0, has_18o(b) == has_18o(e))
                        and _when(cd == 1 and ca > 0, has_13c(a) == has_13c(d))
                        and _when(od > 0 and oa > 0, has_18o(a) == has_18o(d)))
            if (f == "" and sum([n_cogroups(a), n_cogroups(b)])
                    == sum([n_cogroups(d), n_cogroups(e)])
                    and fpA == 0 and daeabs == [0, 0, 0, 0, 1, 0, 0, 0, 0]):
                ok &= (has_13c(a) == has_13c(e) and has_13c(b) == has_13c(d)
                       if is_ion(a) == is_ion(e)
                       else has_13c(b) == has_13c(e) and has_13c(a) == has_13c(d))
            if (f == "" and coga == cogd and fpA == 0 and daeabs != dadabs
                    and not _match_delta(dadabs, {2: 1})):
                if is_ion(a) == is_ion(e):
                    ok &= (_when(ca > 0 and ce > 0, has_13c(a) == has_13c(e))
                        and _when(cd == 1 and cb > 0, has_13c(b) == has_13c(d)))
                else:
                    ok &= (_when(cb > 0 and ce > 0, has_13c(b) == has_13c(e))
                        and _when(cd == 1 and ca > 0, has_13c(a) == has_13c(d)))
            if f != "" and f != "e-":
                ok &= _c18o_together(a, b, c, d, e, f)
            return bool(ok)
        return _filtered(ref, rest, pred)
    if co_r and co_p and nco_r - nco_p == 1:                      # X non-iso-CO
        def pred(a, b, c, d, e, f, g):
            ok = True
            if f == "":
                if count_catoms(a) == count_catoms(d):
                    ok &= has_13c(a) == has_13c(d) and has_13c(b) == has_13c(e)
                else:
                    ok &= has_13c(a) == has_13c(e) and has_13c(b) == has_13c(d)
                ok &= _when(count_catoms(b) == count_catoms(e)
                            and co_bound(b) and co_bound(e)
                            and [d_be_ref[2], d_be_ref[4]] == [0, 1],
                            implies(c18o_bound(e), c18o_bound(b)))
            if (a == "He+" and f == "He") or (a == _HE + "+" and f == _HE):
                ok &= _when(count_oatoms(b) == count_oatoms(d),
                            has_18o(b) == has_18o(d))
                ok &= _when(chn_bound(b) == chn_bound(d), chn_isotop(b) == chn_isotop(d))
                ok &= _when(chn_bound(b) == chn_bound(e), chn_isotop(b) == chn_isotop(e))
            if f != "" and sum(count_catoms(s) for s in (a, b, c)) == 1:
                ok &= (implies(has_18o(e), has_18o(b))
                       if count_catoms(b) == count_catoms(e)
                       else implies(has_18o(d), has_18o(b)))
            return bool(ok)
        return _filtered(ref, rest, pred)
    if (sum(_sgn(count_catoms(s)) for s in (ra, rb, rc)) == 1
            or sum(_sgn(count_catoms(s)) for s in (rd, re_, rf, rg)) == 1):
        def pred(a, b, c, d, e, f, g):                            # bundled C atoms
            ok = True
            if _any3(cp_bound, a, b, c) and _any3(cp_bound, d, e, f, g):
                ok &= (any_of([cp_isotop(s) for s in (a, b, c)])
                    == any_of([cp_isotop(s) for s in (d, e, f, g)]))
            if (co_bound(a) and co_bound(d) and not co_bound(b) and not co_bound(e)
                    and [d_ad_ref[2], d_ad_ref[4]] == [0, -1]):
                ok &= implies(c18o_bound(a), c18o_bound(d))
            if (co_bound(a) and co_bound(e) and not co_bound(b) and not co_bound(d)
                    and [d_ae_ref[2], d_ae_ref[4]] == [0, -1]):
                ok &= implies(c18o_bound(a), c18o_bound(e))
            return bool(ok)
        return _filtered(ref, rest, pred)
    if sum(_sgn(count_catoms(s)) for s in (ra, rb, rc, rd, re_, rf, rg)) >= 4:
        def pred(a, b, c, d, e, f, g):                            # C atoms everywhere
            ok = sticky_13c([a, b, d, e])
            if fpB_ref == 1:
                if d_bx_ref[2] > 0:  ok &= implies(has_13c(d), has_13c(b))
                if d_bx_ref[2] < 0:  ok &= implies(has_13c(b), has_13c(d))
                if d_bx_ref[2] == 0:
                    ok &= has_13c(b) == has_13c(d) and has_13c(a) == has_13c(e)
            elif fpB_ref == 2:
                if d_bx_ref[2] > 0:  ok &= implies(has_13c(e), has_13c(b))
                if d_bx_ref[2] < 0:  ok &= implies(has_13c(b), has_13c(e))
                if d_bx_ref[2] == 0:
                    ok &= has_13c(b) == has_13c(e) and has_13c(a) == has_13c(d)
            else:                                                  # fpB == 0
                dbd = delta_bd(a, b, d, e)
                dbe = delta_be(a, b, d, e)
                ok &= _when(cn_bound(e) and count_natoms(b) == 1
                            and [dbd[2], dbd[3]] != [-1, 1],
                            implies(cn_isotop(e), has_13c(b)))
                ok &= _when(count_natoms(d) == count_natoms(b) == 1
                            and count_catoms(b) >= count_catoms(d),
                            implies(has_13c(d), has_13c(b)))
                ok &= _when(cp_bound(d) and cp_bound(b)
                            and count_catoms(b) >= count_catoms(d),
                            implies(cp_isotop(d), cp_isotop(b)))
                ok &= _when(not cp_bound(b) and cp_bound(d)
                            and count_patoms(b) == count_patoms(d)
                            and count_catoms(b) >= count_catoms(d),
                            implies(cp_isotop(d), has_13c(b)))
                m1 = _match_delta(dbd, {2: 1}) or _match_delta(dbd, {2: 2})
                m2 = _match_delta(dbd, {2: -1}) or _match_delta(dbd, {2: -2})
                ok &= _when(m1 and _norm(dbd) < _norm(dbe),
                            implies(has_13c(d), has_13c(b))
                            and implies(has_13c(a), has_13c(e)))
                ok &= _when(m2 and _norm(dbd) < _norm(dbe),
                            implies(has_13c(b), has_13c(d))
                            and implies(has_13c(e), has_13c(a)))
                ok &= _when(m1 and _norm(dbd) > _norm(dbe),
                            implies(has_13c(e), has_13c(b))
                            and implies(has_13c(a), has_13c(d)))
                ok &= _when(m2 and _norm(dbd) > _norm(dbe),
                            implies(has_13c(e), has_13c(b))
                            and implies(has_13c(a), has_13c(d)))
            return bool(ok)
        return _filtered(ref, rest, pred)
    if rf == "" and fpA_ref == 0:                                 # last filter
        def pred(a, b, c, d, e, f, g):
            ok = sticky_13c([a, b, d, e])
            if d_ae_ref == [-3, 0, -1, 0, 0, 0, 0, 0, 0]:
                if count_catoms(ra) > 0 and count_catoms(re_) > 1:
                    ok &= implies(has_13c(a), has_13c(e))
            return bool(ok)
        return _filtered(ref, rest, pred)
    return list(family)                                           # no match

def clean_up_nn(family):
    ref, rest = family[0], family[1:]
    t, ra, rb, rc, rd, re_, rf, rg = ref[1:9]
    d_ad_ref = delta_ad(ra, rb, rd, re_)
    d_ae_ref = delta_ae(ra, rb, rd, re_)
    d_bd_ref = delta_bd(ra, rb, rd, re_)
    d_be_ref = delta_be(ra, rb, rd, re_)
    d_bx_ref = delta_bx(ra, rb, rd, re_)
    d_ax_ref = delta_ax(ra, rb, rd, re_)
    fpA_ref = find_partner_a(ra, rb, rd, re_)
    fpB_ref = find_partner_b(ra, rb, rd, re_)
    nc_r = sum(count_catoms(s) for s in (ra, rb, rc))

    if ref[1:9] == ["NN", "C", "NCCN", "", "CN", "C2N", "", ""]:  # custom #118
        return _filtered(ref, rest,
                         lambda a, b, c, d, e, f, g: implies(has_13c(a), has_13c(e)))

    handled, result = _common_cases("NN", ref, rest, ra, rb, rc, rd, re_, rf, rg)
    if handled:
        return result

    chn_r = _any3(chn_bound, ra, rb, rc)
    chn_p = _any3(chn_bound, rd, re_, rf)
    nch_r = sum(n_chngroups(s) for s in (ra, rb, rc))
    nch_p = sum(n_chngroups(s) for s in (rd, re_, rf))
    nco_r = sum(n_cogroups(s) for s in (ra, rb, rc))
    nco_p = sum(n_cogroups(s) for s in (rd, re_, rf))

    if (all_of([chn_bound(s) for s in (ra, rb, rd, re_)]) and fpA_ref != 0):
        return _filtered(ref, rest, _pred_15_4chn)                # 1.5 4CHn
    if (chn_r and chn_p and nc_r > 1 and nch_r == nch_p and nco_r == nco_p == 2
            and not (co_bound(rd) and co_bound(re_))):
        return _filtered(ref, rest, _pred_15_2co)                 # 1.5 2CO
    if chn_r and chn_p and nc_r > 1 and nch_r == nch_p < 3 and nco_r == nco_p <= 1:
        return _filtered(ref, rest, _pred_15_2co)                 # 1.5 (simple)
    if (chn_r and chn_p and nc_r > 1 and nch_r == nch_p
            and (nco_r == nco_p == 3 or nco_r != nco_p)
            and fpA_ref != 0 and fpB_ref != 0):                   # 1.5 3/!=CO
        def pred(a, b, c, d, e, f, g):
            if fpA_ref == 1:
                if d_ax_ref[0] != 0:
                    return (chn_isotop(a) == chn_isotop(d)
                        and co_isotop(a) == co_isotop(d)
                        and has_18o(a) == has_18o(d)
                        and left18o_bound(a) == left18o_bound(d)
                        and left18o_bound(b) == left18o_bound(e)
                        and cn_isotop(a) == cn_isotop(d))
                if d_ax_ref[4] != 0:
                    return (chn_isotop(a) == chn_isotop(d)
                        and cn_isotop(a) == cn_isotop(d)
                        and left18o_bound(a) == left18o_bound(d)
                        and left18o_bound(b) == left18o_bound(e))
                return True
            else:
                if d_ax_ref[0] != 0:
                    return (chn_isotop(a) == chn_isotop(e)
                        and co_isotop(a) == co_isotop(e)
                        and has_18o(a) == has_18o(e)
                        and cn_isotop(a) == cn_isotop(e)
                        and left18o_bound(a) == left18o_bound(e)
                        and left18o_bound(b) == left18o_bound(d))
                if d_ax_ref[4] != 0:
                    return (chn_isotop(a) == chn_isotop(e)
                        and cn_isotop(a) == cn_isotop(e)
                        and left18o_bound(a) == left18o_bound(e)
                        and left18o_bound(b) == left18o_bound(d))
                return True
        return _filtered(ref, rest, pred)
    if (chn_r and chn_p and nc_r > 1 and rf == ""
            and _any3(co_bound, ra, rb, rc) == _any3(co_bound, rd, re_, rf)
            and nch_r != nch_p and fpA_ref != 0):                 # case 1.6
        def pred(a, b, c, d, e, f, g):
            dax = delta_ax(a, b, d, e)
            axabs = [abs(x) for x in dax]
            fpA = find_partner_a(a, b, d, e)
            ok = True
            if dax[2] == 0:
                ok &= ([has_13c(a), has_13c(b)] == [has_13c(d), has_13c(e)] if fpA == 1
                       else [has_13c(a), has_13c(b)] == [has_13c(e), has_13c(d)])
            if dax[4] == 0:
                ok &= ([has_18o(a), has_18o(b)] == [has_18o(d), has_18o(e)] if fpA == 1
                       else [has_18o(a), has_18o(b)] == [has_18o(e), has_18o(d)])
            if fpA == 1:
                dbx = delta_bx(a, b, d, e)
                h1, h3 = dax[0], dax[2]
                a1, a3 = axabs[0], axabs[2]
                ok &= _when([a1, a3] != [3, 1] and chn_bound(a) == chn_bound(d),
                            chn_isotop(a) == chn_isotop(d))
                ok &= _when([a1, a3] != [3, 1] and chn_bound(b) == chn_bound(e),
                            chn_isotop(b) == chn_isotop(e))
                ok &= _when([a1, a3] == [3, 1] and chn_bound(a) == chn_bound(d),
                            implies(chn_isotop(d), chn_isotop(a)))
                ok &= _when([a1, a3] == [3, 1] and chn_bound(b) == chn_bound(e),
                            implies(chn_isotop(e), chn_isotop(b)))
                ok &= _when([h1, h3] == [-3, -1] and co_bound(a) == co_bound(d),
                            co_isotop(d) == co_isotop(a))
                ok &= _when([dbx[0], dbx[2]] == [3, 1] and co_bound(b) == co_bound(e),
                            co_isotop(b) == co_isotop(e))
                ok &= _when([h1, h3] == [3, 1] and co_bound(a) == co_bound(d),
                            co_isotop(d) == co_isotop(a))
                ok &= _when([dbx[0], dbx[2]] == [-3, -1] and co_bound(b) == co_bound(e),
                            co_isotop(b) == co_isotop(e))
            elif fpA == 2:
                dbx = delta_bx(a, b, d, e)
                h1, h3 = dax[0], dax[2]
                a1, a3 = axabs[0], axabs[2]
                ok &= _when([a1, a3] != [3, 1] and chn_bound(a) == chn_bound(e),
                            chn_isotop(a) == chn_isotop(e))
                ok &= _when([a1, a3] != [3, 1] and chn_bound(b) == chn_bound(d),
                            chn_isotop(b) == chn_isotop(d))
                ok &= _when([a1, a3] == [3, 1] and chn_bound(a) == chn_bound(e),
                            implies(chn_isotop(e), chn_isotop(a)))
                ok &= _when([a1, a3] == [3, 1] and chn_bound(b) == chn_bound(d),
                            implies(chn_isotop(d), chn_isotop(b)))
                ok &= _when([h1, h3] == [-3, -1] and co_bound(a) == co_bound(e),
                            co_isotop(e) == co_isotop(a))
                ok &= _when([dbx[0], dbx[2]] == [3, 1] and co_bound(b) == co_bound(d),
                            co_isotop(b) == co_isotop(d))
                ok &= _when([h1, h3] == [3, 1] and co_bound(a) == co_bound(e),
                            co_isotop(e) == co_isotop(a))
                ok &= _when([dbx[0], dbx[2]] == [-3, -1] and co_bound(b) == co_bound(d),
                            co_isotop(b) == co_isotop(d))
            return bool(ok)
        return _filtered(ref, rest, pred)

    handled, result = _string_overlap_cases(ref, rest, ra, rb, rc, rd, re_, rf, rg)
    if handled:
        return result

    # ---- remaining NN-only cases ---------------------------------------- #
    nrd, nre, nra, nrb = neutral(rd), neutral(re_), neutral(ra), neutral(rb)
    if "C2N" in (ra, rb, rc) or "C2O" in (ra, rb, rc):            # 5b C2X
        pred = lambda a, b, c, d, e, f, g: (
            implies(has_13c(e), has_13c(b)) and implies(has_18o(a), has_18o(d))
            and implies(has_13c(a), has_13c(d)) and implies(has_18o(e), has_18o(b)))
        return _filtered(ref, rest, pred)
    if ((len(rb) >= 3 and len(re_) >= 3 and _ext_pref(nrb, nre))
            or (len(rb) >= 3 and len(rd) >= 3 and _ext_pref(nrb, nrd))):  # 5b
        return _filtered(ref, rest, _pred_5b)
    co_r = _any3(co_bound, ra, rb, rc)
    co_p = _any3(co_bound, rd, re_, rf)
    nco_r = sum(n_cogroups(s) for s in (ra, rb, rc))
    nco_p = sum(n_cogroups(s) for s in (rd, re_, rf))
    if co_r and co_p and nco_r == nco_p:                          # X iso-CO
        def pred(a, b, c, d, e, f, g):
            dbx = delta_bx(a, b, d, e); dbxabs = [abs(x) for x in dbx]
            dae = delta_ae(a, b, d, e); daeabs = [abs(x) for x in dae]
            fpA = find_partner_a(a, b, d, e)
            coga = sorted([n_cogroups(a), n_cogroups(b)])
            cogd = sorted([n_cogroups(d), n_cogroups(e)])
            unit = [1, 0, 0, 0, 0, 0, 0, 0, 0]
            ok = True
            if f == "" and coga == cogd and _match_delta(dbxabs, {4: 0}):
                ok &= _c18o_together(a, b, c, d, e, f)
            if f == "" and coga != cogd and fpA == 1 and _match_delta(
                    [abs(x) for x in d_be_ref], {4: 1}):
                ok &= (implies(has_18o(e), has_18o(b))
                    and implies(has_13c(e), has_13c(b))
                    and implies(has_13c(d), has_13c(a)))
            if f == "" and coga == cogd and fpA == 0 and daeabs == unit:
                ok &= (has_18o(a) == has_18o(e) and has_18o(b) == has_18o(d)
                    and implies(has_13c(e), has_13c(a))
                    and implies(has_13c(d), has_13c(b)))
            if f != "":
                ok &= _c18o_together(a, b, c, d, e, f)
            return bool(ok)
        return _filtered(ref, rest, pred)
    if co_r and co_p and nco_r - nco_p == 1:                      # X non-iso-CO
        def pred(a, b, c, d, e, f, g):
            # NOTE: the original's second If has no else-branch, so rows with
            # f != "" drop out of the condition chain (And[..., Null, ...])
            if f != "":
                return False
            ok = True
            if count_catoms(a) == count_catoms(d):
                ok &= has_13c(a) == has_13c(d) and has_13c(b) == has_13c(e)
            else:
                ok &= has_13c(a) == has_13c(e) and has_13c(b) == has_13c(d)
            ok &= _when(count_oatoms(a) + 1 == count_oatoms(d),
                        implies(has_18o(a), has_18o(d)))
            if (a == "He+" and f == "He") or (a == _HE + "+" and f == _HE):
                ok &= _when(count_oatoms(b) == count_oatoms(d),
                            has_18o(b) == has_18o(d))
                ok &= _when(chn_bound(b) == chn_bound(d), chn_isotop(b) == chn_isotop(d))
                ok &= _when(chn_bound(b) == chn_bound(e), chn_isotop(b) == chn_isotop(e))
            return bool(ok)
        return _filtered(ref, rest, pred)
    if (sum(_sgn(count_catoms(s)) for s in (ra, rb, rc)) == 1
            or sum(_sgn(count_catoms(s)) for s in (rd, re_, rf, rg)) == 1):
        def pred(a, b, c, d, e, f, g):                            # bundled C atoms
            ok = sticky_18o([a, b, c, d, e, f, g])
            if _any3(cp_bound, a, b, c) and _any3(cp_bound, d, e, f, g):
                ok &= (any_of([cp_isotop(s) for s in (a, b, c)])
                    == any_of([cp_isotop(s) for s in (d, e, f, g)]))
            if (co_bound(a) and co_bound(d) and not co_bound(b) and not co_bound(e)
                    and [d_ad_ref[2], d_ad_ref[4]] == [0, -1]):
                ok &= implies(c18o_bound(a), c18o_bound(d))
            return bool(ok)
        return _filtered(ref, rest, pred)
    if (sum(_sgn(count_catoms(s)) for s in (ra, rb, rc, rd, re_, rf, rg)) >= 4
            and (fpA_ref != 0 or fpB_ref != 0)):
        def pred(a, b, c, d, e, f, g):                            # C atoms everywhere
            ok = sticky_18o([a, b, c, d, e, f, g]) and sticky_13c([a, b, d, e])
            if fpB_ref == 1:
                if d_bx_ref[2] > 0:  ok &= implies(has_13c(d), has_13c(b))
                if d_bx_ref[2] < 0:  ok &= implies(has_13c(b), has_13c(d))
                if d_bx_ref[2] == 0:
                    ok &= has_13c(b) == has_13c(d) and has_13c(a) == has_13c(e)
            else:                                                  # fpB == 2
                if d_bx_ref[2] > 0:  ok &= implies(has_13c(e), has_13c(b))
                if d_bx_ref[2] < 0:  ok &= implies(has_13c(b), has_13c(e))
                if d_bx_ref[2] == 0:
                    ok &= has_13c(b) == has_13c(e) and has_13c(a) == has_13c(d)
            return bool(ok)
        return _filtered(ref, rest, pred)
    if rf == "" and fpA_ref == 0:                                 # last filter
        def pred(a, b, c, d, e, f, g):
            ok = sticky_18o([a, b, c, d, e, f, g]) and sticky_13c([a, b, d, e])
            if d_ae_ref == [-3, 0, -1, 0, 0, 0, 0, 0, 0]:
                if count_catoms(ra) > 0 and count_catoms(re_) > 1:
                    ok &= implies(has_13c(a), has_13c(e))
            return bool(ok)
        return _filtered(ref, rest, pred)
    return list(family)                                           # no match

# --------------------------------------------------------------------------- #
# pipeline
# --------------------------------------------------------------------------- #

_CLEANERS = {"NN": clean_up_nn, "MN": clean_up_mn, "IN": clean_up_in,
             "RA": clean_up_ra, "AD": clean_up_ad, "DR": clean_up_dr,
             "CE": clean_up_ce, "CL": clean_up_cl,
             "PH": clean_up_photo, "CR": clean_up_photo, "CP": clean_up_photo}

def insert_isotopes(rows, verbose=False):
    """InsertIsotopes: generate 13C/18O families and prune them per type."""
    updated = [[r[0], r[1]] + [_norm_to_double(x) for x in r[2:9]] + list(r[9:])
               for r in rows]
    families = {}
    for row in insert_13c_o18(updated):
        families.setdefault(row[0], []).append(row)
    out = []
    numbers = sorted(families)
    for idx, num in enumerate(numbers):
        if verbose and idx % 500 == 0:
            print(f"Processing entry {idx + 1} of {len(numbers)}")
        fam = sorted(families[num], key=_row_key, reverse=True)   # original first
        cleaner = _CLEANERS.get(fam[0][1])
        if cleaner is None:
            print("Unknown reaction", num)
            cleaned = fam
        else:
            cleaned = cleaner(fam)
        out.extend(cleaned)
    return [[r[0], r[1]] + [_double_to_norm(x) for x in r[2:9]] + list(r[9:])
            for r in out]

def rescale_rates(rows):
    """Uniform branching: alpha is divided by the number of consecutive rows
    sharing identical {n, t, reactants} (i.e. product branches of the same
    reactant-side isotopologue)."""
    out = [list(r) for r in rows]
    i = 0
    while i < len(out):
        j = i
        while j + 1 < len(out) and out[j + 1][:5] == out[i][:5]:
            j += 1
        branches = j - i + 1
        for k in range(i, j + 1):
            out[k][ALPHA] = out[k][ALPHA] / branches
        i = j + 1
    return out

# --------------------------------------------------------------------------- #
# FORTRAN77 output helpers
# --------------------------------------------------------------------------- #

def pad_string(s: str, width: int) -> str:
    if len(s) <= width:
        return s.ljust(width)
    print("Warning! Input String is too large.")
    return s   # deviation: original returned Null here

def _round_half_away(x: float, nd: int) -> float:
    c = 10 ** nd
    v = abs(x) * c
    f = math.floor(v + 0.5)
    return math.copysign(f / c, x) if x else 0.0

def printf2(x, m: int, n: int) -> str:
    """FORMAT(Fm.n)-style fixed point (NumberForm semantics: m-2 total
    digits, n after the point -> integer part padded with spaces)."""
    r = _round_half_away(float(x), n)
    s = f"{r:.{n}f}"
    ip, _, fp = s.partition(".")
    return ip.rjust(max(0, m - 2 - n)) + "." + fp

def _mantissa_exponent10(x: float):
    if x == 0:
        return 0.0, 0
    exp = math.floor(math.log10(abs(x)))
    mant = x / 10 ** exp
    if abs(mant) >= 10:
        mant, exp = mant / 10, exp + 1
    elif abs(mant) < 1:
        mant, exp = mant * 10, exp - 1
    return mant, exp

def eformat(x, tot: int = 9, ndigits: int = 2) -> str:
    """Fortran E representation (faithful incl. the mantissa*10 quirk and
    the zero test against 0 instead of $MachineEpsilon)."""
    if tot - ndigits < 7:
        raise ValueError("Eformat: too many digits to the right of the "
                         "decimal point requested")
    x = float(x)
    mant, exp = _mantissa_exponent10(x)
    if mant != 1.0:                    # original transformation (sic)
        mant, exp = mant * 10.0, exp - 1
    lead = "" if x < 0 else " "
    s = f"{mant:.{ndigits}f}"
    man_str = s.rjust(tot - ndigits) + "." + s.partition(".")[2] if "." in s else s
    sign = "+" if exp >= 0 else "-"
    return f"{lead}{man_str}E{sign}{abs(exp)}"

def fortran_out(row, duplicated=()) -> str:
    """fortranOut: fixed format for cmread.f (widened reactant fields)."""
    n, t, a, b, c, d, e, f, g, h, i, j, *rest = row
    lead = "*" if n in duplicated else " "
    return (lead + str(n).rjust(5) + " " + pad_string(a, 15) + " "
            + pad_string(b, 15) + " " + pad_string(d, 15) + " "
            + pad_string(e, 15) + " " + pad_string(f, 10) + " "
            + pad_string(g, 9) + eformat(h)[1:] + " "          # Rest@Eformat
            + printf2(i, 5, 2) + " " + printf2(j, 8, 1))

def frac_out(row, duplicated=()) -> str:
    """fracOut: same layout for fractionation-reaction files."""
    n, a, b, d, e, h, i, j, *rest = row
    lead = "*" if n in duplicated else " "
    return (lead + str(n).rjust(5) + " " + pad_string(a, 15) + " "
            + pad_string(b, 15) + " " + pad_string(d, 15) + " "
            + pad_string(e, 15) + " " + pad_string("", 10) + " "
            + pad_string("", 9) + eformat(h)[1:] + " "
            + printf2(i, 5, 2) + " " + printf2(j, 8, 1))

# --------------------------------------------------------------------------- #
# CSV I/O
# --------------------------------------------------------------------------- #

def read_umist_csv(path):
    """Minimal reader for UMIST-style CSV files (>= 13 columns)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if not parts[0].lstrip("-").isdigit():
                continue                       # header / malformed line
            row = [int(parts[0]), parts[1]] + parts[2:9]
            row += [float(parts[9]), float(parts[10]), float(parts[11])]
            row += parts[12:]
            rows.append(row)
    return rows

# --------------------------------------------------------------------------- #
# demo
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    demo = [
        [1, "NN", "C",   "OH", "", "CO",  "H",  "",   "", 1.0e-10, 0.0,  0.0],
        [2, "IN", "H3+", "CO", "", "HCO+", "H2", "",   "", 1.7e-9,  0.3,  0.0],
        [3, "DR", "HCO+","e-", "", "CO",  "H",  "",   "", 1.1e-7, -0.7, 0.0],
        [4, "NN", "CH3", "CH3","", "CH4", "CH2","",   "", 6.1e-11, 0.0,  0.0],
    ]
    result = rescale_rates(insert_isotopes(demo, verbose=True))
    for r in result:
        print(fortran_out(r))

