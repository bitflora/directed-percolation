"""Validation shared by exact-series analysis entry points."""
import json


def load_record(path, key, allow_unverified=False):
    with open(path) as f:
        ser = json.load(f)[key]
    if ser.get("verified") is not True and not allow_unverified:
        raise ValueError("series is not CRT-verified; use --allow-unverified for exploratory analysis")
    return ser


def validate_order(ser, order, names):
    if order < 1 or order > ser["order"]:
        raise ValueError(f"requested order {order}; available order is {ser['order']}")
    for name in names:
        if len(ser[name]) < order + 1:
            raise ValueError(f"{name}: need {order + 1} coefficients, have {len(ser[name])}")
