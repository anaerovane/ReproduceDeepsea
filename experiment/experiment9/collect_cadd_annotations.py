#!/usr/bin/env python3
"""Poll and collect finished CADD GRCh37-v1.4 inclAnno batches."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import re
import time
from urllib.parse import urljoin

import requests

HERE = Path(__file__).resolve().parent
OUTPUTS = HERE / "outputs"


class VisibleText(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts = []

    def handle_data(self, data):
        data = ' '.join(data.split())
        if data:
            self.parts.append(data)


def save_status(path, payload):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(payload, indent=2) + '\n')
    temp.replace(path)


def check_job(job):
    index = job['index']
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    dest = OUTPUTS / f"cadd_inclAnno_batch_{index:02d}.tsv.gz"
    if dest.exists() and dest.stat().st_size:
        try:
            with gzip.open(dest, 'rb') as f:
                f.read(1)
            return str(index), {'state': 'downloaded', 'bytes': dest.stat().st_size}
        except (OSError, EOFError):
            dest.unlink(missing_ok=True)

    try:
        response = requests.get(job['path'], timeout=(10, 40))
        body = response.content
        if response.status_code != 200:
            return str(index), {'state': 'http_error', 'http': response.status_code,
                                'detail': body[:300].decode('utf-8', 'replace')}

        if body.startswith(b'\x1f\x8b'):

            with gzip.GzipFile(fileobj=io.BytesIO(body)) as f:
                f.read(1)
            temp = dest.with_suffix(dest.suffix + '.tmp')
            temp.write_bytes(body)
            temp.replace(dest)
            return str(index), {'state': 'downloaded', 'bytes': dest.stat().st_size}

        parser = VisibleText()
        parser.feed(body.decode('utf-8', 'replace'))
        detail = ' '.join(parser.parts)
        finished = re.search(r'href=["\']([^"\']*static/finished/[^"\']+\.tsv\.gz)["\']',
                             body.decode('utf-8', 'replace'), re.IGNORECASE)
        if finished:
            result_url = urljoin(response.url, finished.group(1))
            result = requests.get(result_url, timeout=(10, 120))
            result.raise_for_status()
            payload = result.content
            if payload.startswith(b'\x1f\x8b'):
                with gzip.GzipFile(fileobj=io.BytesIO(payload)) as f:
                    f.read(1)
            elif b'\t' in payload[:1000] and not payload.lstrip().startswith(b'<'):
                payload = gzip.compress(payload)
            else:
                return str(index), {'state': 'download_error', 'http': result.status_code,
                                    'detail': f'Finished page linked to non-TSV result: {result_url}'}
            temp = dest.with_suffix(dest.suffix + '.tmp')
            temp.write_bytes(payload)
            temp.replace(dest)
            return str(index), {'state': 'downloaded', 'bytes': dest.stat().st_size,
                                'source_url': result_url}
        if 'not yet finished' in detail.lower() or 'recheck in a few minutes' in detail.lower():
            marker = detail.lower().find('successfully uploaded variants')
            if marker < 0:
                marker = detail.lower().find('not yet finished')
            detail = detail[marker:marker + 500] if marker >= 0 else detail[-500:]
            return str(index), {'state': 'processing', 'http': 200,
                                'detail': detail[:700]}
        return str(index), {'state': 'unexpected_response', 'http': 200,
                            'content_type': response.headers.get('content-type'),
                            'detail': detail[:700]}
    except Exception as exc:
        return str(index), {'state': 'error', 'error': f'{type(exc).__name__}: {exc}'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--once', action='store_true', help='check all jobs once and exit')
    parser.add_argument('--interval', type=int, default=300,
                        help='seconds between polling rounds (default: 300)')
    parser.add_argument('--workers', type=int, default=4,
                        help='parallel availability checks (default: 4)')
    args = parser.parse_args()
    if args.interval < 10 or args.workers < 1:
        parser.error('--interval must be >=10 and --workers must be >=1')

    manifest = json.loads((HERE / 'cadd_annotation_jobs.json').read_text())
    status_path = HERE / 'cadd_annotation_status.json'
    jobs = manifest['jobs']

    while True:
        statuses = {}
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(check_job, job) for job in jobs]
            for future in as_completed(futures):
                key, value = future.result()
                statuses[key] = value
                save_status(status_path, {
                    'source': manifest['source'],
                    'last_checked_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
                    'completed': sum(v['state'] == 'downloaded' for v in statuses.values()),
                    'checked_this_round': len(statuses), 'total': len(jobs),
                    'jobs': {**statuses},
                })

        completed = sum(v['state'] == 'downloaded' for v in statuses.values())
        processing = sum(v['state'] == 'processing' for v in statuses.values())
        errors = sum(v['state'] in {'error', 'http_error', 'unexpected_response', 'download_error'}
                     for v in statuses.values())
        print(f"{datetime.now(timezone.utc).isoformat(timespec='seconds')} "
              f"CADD batches: {completed}/{len(jobs)} downloaded; "
              f'{processing} processing; {errors} errors', flush=True)
        if processing:
            example = next(v['detail'] for v in statuses.values()
                           if v['state'] == 'processing')
            print(f'Server status: {example[:350]}', flush=True)
        if completed == len(jobs) or args.once:
            break
        time.sleep(args.interval)


if __name__ == '__main__':
    main()
