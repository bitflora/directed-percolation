"""Compare a wall series (crt.py output of dpwall) with Table 1 of Essam,
Guttmann, Jensen and TanlaKishani, J. Phys. A 29 (1996) 1619.

Usage: python check_wall_table.py results/series/wall110.json [key]
"""
import json
import sys

COLS = ["T", "S", "mt1", "mt2", "mx1", "mx2"]


def main():
    d = json.load(open(sys.argv[1]))
    s = d[sys.argv[2]] if len(sys.argv) > 2 else next(iter(d.values()))
    table = __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0] + "/egjt96_wall_table.txt"
    bad, rows = [], 0
    for line in open(table):
        if line.startswith("#") or not line.strip():
            continue
        v = [int(x) for x in line.split()]
        n = v[0]
        if n > s["order"]:
            break
        rows += 1
        for c, x in zip(COLS, v[1:]):
            if s[c][n] != x:
                bad.append((c, n, x, s[c][n]))
    print(f"{rows} table rows compared; {'all entries match' if not bad else f'{len(bad)} differ: {bad[:5]}'}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
