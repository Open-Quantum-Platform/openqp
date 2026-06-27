#!/usr/bin/env python3
"""Deterministic code<->manual consistency check for input keywords.

Compares the ``default`` value of every keyword declared in
``pyoqp/oqp/molecule/oqpdata.py`` (OQP_CONFIG_SCHEMA) against the value
documented in the separate user manual (openqp-docs, ``docs/keywords/<section>.md``),
and flags keywords that exist in code but are undocumented.

This is the reliable, non-AI backstop for the most mechanical kind of drift
(a changed default, or a new keyword that nobody documented). It complements the
agent PR review, which handles the semantic/cross-cutting cases.

Usage:
    python tools/check_doc_consistency.py --docs /path/to/openqp-docs
    python tools/check_doc_consistency.py            # auto-locate sibling openqp-docs

Exit status: 0 if consistent (or only advisory notes), 1 if a concrete default
mismatch or an undocumented keyword is found.
"""
import argparse
import ast
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SCHEMA_FILE = os.path.join(REPO, 'pyoqp', 'oqp', 'molecule', 'oqpdata.py')

# Sections whose keyword pages we do not mechanically diff (free-form / no table).
SKIP_SECTIONS = {'tests'}


def schema_defaults():
    """{section: {option: default_string}} parsed from OQP_CONFIG_SCHEMA via AST."""
    tree = ast.parse(open(SCHEMA_FILE).read())
    schema_node = None
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'OQP_CONFIG_SCHEMA'
                        for t in node.targets)):
            schema_node = node.value
            break
    if schema_node is None or not isinstance(schema_node, ast.Dict):
        raise SystemExit('could not find OQP_CONFIG_SCHEMA dict in ' + SCHEMA_FILE)

    out = {}
    for sec_key, sec_val in zip(schema_node.keys, schema_node.values):
        section = ast.literal_eval(sec_key)
        if not isinstance(sec_val, ast.Dict):
            continue
        opts = {}
        for opt_key, opt_val in zip(sec_val.keys, sec_val.values):
            option = ast.literal_eval(opt_key)
            if not isinstance(opt_val, ast.Dict):
                continue
            for k, v in zip(opt_val.keys, opt_val.values):
                if ast.literal_eval(k) == 'default':
                    try:
                        opts[option] = str(ast.literal_eval(v))
                    except Exception:
                        opts[option] = None
        out[section] = opts
    return out


def manual_defaults(docs_dir):
    """{section: {keyword: documented_default_or_None}} from docs/keywords/*.md."""
    kdir = os.path.join(docs_dir, 'docs', 'keywords')
    out = {}
    if not os.path.isdir(kdir):
        raise SystemExit('keyword docs not found: ' + kdir)
    # A heading may document several keywords at once, e.g.
    # "### `cam_alpha`, `cam_beta`, `cam_mu`" or "### `step_size`, `step_tol`".
    is_heading = re.compile(r'^#{2,4}\s+`')
    backticked = re.compile(r'`([^`]+)`')
    default_row = re.compile(r'^\|\s*Defaults?\s*\|\s*(.+?)\s*\|', re.I)
    for fname in os.listdir(kdir):
        if not fname.endswith('.md') or fname == 'index.md':
            continue
        section = fname[:-3]
        kw = {}
        current = []
        for line in open(os.path.join(kdir, fname)):
            if is_heading.match(line):
                current = [k.strip() for k in backticked.findall(line)]
                for k in current:
                    kw.setdefault(k, None)
                continue
            if current:
                m = default_row.match(line)
                if m:
                    # combined heading -> the same Default row applies to each;
                    # a single token gets the value, multiple get "documented".
                    val = m.group(1).strip()
                    for k in current:
                        if kw.get(k) is None:
                            kw[k] = val if len(current) == 1 else '<documented>'
        out[section] = kw
    return out


def normalize(val):
    """Normalize a default for comparison; return None if not concretely comparable."""
    if val is None:
        return None
    if val == '<documented>':
        return None  # present via a combined heading; specific default not parseable
    v = val.strip().strip('`').strip()
    low = v.lower()
    # documented "empty"/"none" prose <-> code ''
    if low.startswith('empty') or low in ('', 'none', 'empty string'):
        return ''
    if low in ('true', 'false'):
        return low
    # numbers: compare by value so 1e-6 == 1.0e-6, 0 == 0.0, etc.
    try:
        return ('num', float(v))
    except ValueError:
        pass
    # descriptive prose (spaces, parentheses) is not a concrete token -> skip
    if ' ' in v or '(' in v:
        return None
    return low


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--docs', help='path to the openqp-docs checkout')
    ap.add_argument('--strict', action='store_true',
                    help='also fail when a code keyword is undocumented '
                         '(default: undocumented keywords are advisory only, '
                         'since the manual has known pre-existing gaps)')
    args = ap.parse_args()

    docs = args.docs
    if not docs:
        for cand in (os.path.join(REPO, os.pardir, 'openqp-docs'),
                     os.path.join(REPO, os.pardir, os.pardir, 'openqp-docs')):
            if os.path.isdir(cand):
                docs = cand
                break
    if not docs or not os.path.isdir(docs):
        raise SystemExit('openqp-docs not found; pass --docs /path/to/openqp-docs')

    schema = schema_defaults()
    manual = manual_defaults(docs)

    mismatches, undocumented = [], []
    for section, opts in schema.items():
        if section in SKIP_SECTIONS:
            continue
        doc_kw = manual.get(section)
        if doc_kw is None:
            continue  # no dedicated keyword page (e.g. merged/aliased section)
        for option, code_default in opts.items():
            if option not in doc_kw:
                undocumented.append(f'[{section}] {option}  (default={code_default!r})')
                continue
            cd, dd = normalize(code_default), normalize(doc_kw[option])
            if cd is not None and dd is not None and cd != dd:
                mismatches.append(
                    f'[{section}] {option}: code default {code_default!r} != '
                    f'manual {doc_kw[option]!r}')

    if mismatches:
        print('Default-value mismatches (code vs manual):')
        for m in mismatches:
            print('  FAIL ' + m)
    if undocumented:
        label = 'FAIL' if args.strict else 'note'
        print(f'Keywords in code schema but not documented ({label}):')
        for u in undocumented:
            print(f'  {label} ' + u)
    if not mismatches and not undocumented:
        print('OK: every schema keyword is documented and defaults agree.')
        return 0
    if not mismatches and not args.strict:
        print('OK: defaults agree (undocumented keywords above are advisory).')
        return 0
    return 1


if __name__ == '__main__':
    sys.exit(main())
