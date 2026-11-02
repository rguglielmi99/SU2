#!/usr/bin/env python3
"""Export Forward Euler temporal-order diagnostics from the Catch2 log."""
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'Docs')
    args = parser.parse_args()
    lines = args.log.read_text().splitlines()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for tag, suffix in [('TEMPORAL_MMS', 'temporal'), ('COUPLED_TEMPORAL_MMS', 'coupled_temporal')]:
        records = [line[len(tag)+1:] for line in lines if line.startswith(tag+',')]
        if len(records) < 2:
            raise ValueError('Missing complete data for '+tag)
        (args.output_dir / ('primitive_euler_phase5_'+suffix+'.csv')).write_text('\n'.join(records)+'\n')


if __name__ == '__main__':
    main()
