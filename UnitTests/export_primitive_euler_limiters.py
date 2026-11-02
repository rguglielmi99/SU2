#!/usr/bin/env python3
"""Export the production-limiter audit from a Catch2 log; no plotting dependency."""
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
    for tag, suffix in [('LIMITER_MMS', 'mms'), ('LIMITER_ACTIVITY', 'activity'),
                        ('LIMITER_EXTREMA', 'extrema'), ('LIMITER_STRESS', 'stress')]:
        records = [line[len(tag)+1:] for line in lines if line.startswith(tag+',')]
        rows = list(csv.DictReader(io.StringIO('\n'.join(records))))
        if not rows:
            raise ValueError('Missing audit records: '+tag)
        target = args.output_dir / ('primitive_euler_phase4_limiter_'+suffix+'.csv')
        target.write_text('\n'.join(records)+'\n')
        print(str(target)+': '+str(len(rows))+' data rows')


if __name__ == '__main__':
    main()
