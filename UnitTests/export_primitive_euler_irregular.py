#!/usr/bin/env python3
"""Export the phase-4 irregular-mesh audit from a focused/full Catch2 log."""
import argparse
import csv
import io
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'Docs')
    args = parser.parse_args()
    lines = args.log.read_text().splitlines()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for tag, suffix in [('IRREGULAR_MMS', 'mms'), ('IRREGULAR_GRAD', 'gradients'),
                        ('CONTACT_SPLIT', 'contact'), ('CONTACT_WEAK', 'weak')]:
        prefix = tag + ','
        records = [line[len(prefix):] for line in lines if line.startswith(prefix)]
        rows = list(csv.DictReader(io.StringIO('\n'.join(records))))
        if not rows:
            raise ValueError('Missing audit records: ' + tag)
        target = args.output_dir / ('primitive_euler_phase4_irregular_' + suffix + '.csv')
        target.write_text('\n'.join(records) + '\n')
        print(str(target) + ': ' + str(len(rows)) + ' data rows')


if __name__ == '__main__':
    main()
