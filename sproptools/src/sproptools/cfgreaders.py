#! /usr/bin/env python3
import argparse
import configparser
import re
from typing import Dict, Tuple, Union

Number = Union[int, float]

"""
Tools for parsing configuration files used in solarprop analysis:

  - parameter scans
  - ...

"""

def parse_range(value: str) -> Tuple[Number, Number]:
    """
    Parse a 'range' string into a 2-tuple of numbers (int or float).

    Accepted formats (whitespace optional):
      - "a, b"
      - "(a, b)"
      - "[a, b]"
      - "a b"

    Examples:
      "1e-4, 1e-2" -> (0.0001, 0.01)
      "(16, 256)"  -> (16, 256)
      "[0.0 0.99]" -> (0.0, 0.99)
    """
    s = value.strip()
    # Strip optional brackets/parentheses
    if (s.startswith('(') and s.endswith(')')) or (s.startswith('[') and s.endswith(']')):
        s = s[1:-1].strip()

    # Split by comma or whitespace
    parts = [p for p in re.split(r'[,\s]+', s) if p]
    if len(parts) != 2:
        raise ValueError(f'range must contain exactly two numbers, got: {value!r}')

    def to_number(x: str) -> Number:
        try:
            return int(x)
        except ValueError:
            return float(x)

    a, b = map(to_number, parts)

    if a > b:
        raise ValueError(f'range lower bound must be <= upper bound, got: {a} > {b}')
    return (a, b)


def read_parameter_scan_configuration(config_path: str) -> Tuple[Dict[str, Tuple[Number, Number]], Dict[str, float]]:
    """
    Parse the configuration file located at `config_path` for parameter scan.

    Read all sections as parameter names and parse their 'range' into tuples.
    Also read start values (for emcee), or calculate them as the mid-point of the range.

    Example cfg file:

    ```
    [kappaScaling]
    range = [0.6, 2.2]
    start_value = 1.6

    [driftFactor]
    range = [0., 1.]

    ```

    Returns two dicts: {parameter_name: (low, high)} and {parameter_name: start_value}.
    """

    cfg = configparser.ConfigParser()
    with open(config_path, 'r', encoding='utf-8') as f:
        cfg.read_file(f)
        cfg.optionxform = str

    params_ranges: Dict[str, Tuple[Number, Number]] = {}
    params_startvals: Dict[str, float] = {}
    for section in cfg.sections():
        if not cfg.has_option(section, 'range'):
            raise KeyError(f'Section [{section}] is missing required "range" option')
        p_range = parse_range(cfg.get(section, 'range'))
        params_ranges[section] = p_range
        if 'start_value' in cfg[section]:
            params_startvals[section] = cfg.getfloat(section, 'start_value')
        else:
            params_startvals[section] = 0.5 * float(p_range[0] + p_range[1])
            print(f'Using midpoint start value of {params_startvals[section]:g} for parameter "{section}".')

    return params_ranges, params_startvals


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Read parameter ranges and start values from an INI config.')
    parser.add_argument('config', help='Path to INI file (e.g., params.ini)')
    args = parser.parse_args()

    params_ranges, start_values = read_parameter_scan_configuration(args.config)
    for name, p_range in params_ranges.items():
        print(f'{name}: {p_range} {start_values[name]}')
