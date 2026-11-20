#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch runner for bolum_block_v2.py.

Place this script in /bolum/ alongside bolum_block_v2.py.
It discovers event folders in /bolum/Events/, runs each event independently,
and writes batch diagnostics and a CSV summary.

Usage:
    python3 run_all_events.py
    python3 run_all_events.py --no-plots
    python3 run_all_events.py --quiet-model
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import subprocess
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
MODEL = SCRIPT_DIR / "bolum_block_v2.py"
EVENTS_DIR = SCRIPT_DIR / "Events"
LOG_DIR = SCRIPT_DIR / "batch_logs"


def timestamp():
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


def main():
    parser = argparse.ArgumentParser(description="Run BOLUM for every event folder in Events/.")
    parser.add_argument("--events-dir", type=Path, default=EVENTS_DIR,
                        help=f"Events directory (default: {EVENTS_DIR})")
    parser.add_argument("--model", type=Path, default=MODEL,
                        help=f"Model script (default: {MODEL})")
    parser.add_argument("--no-plots", action="store_true",
                        help="Pass --no-plots to the model (faster; no PDF output).")
    parser.add_argument("--quiet-model", action="store_true",
                        help="Pass --quiet to the model (batch runner still reports diagnostics).")
    parser.add_argument("--sat-W", type=float, default=None,
                        help="Optional sensor saturation level passed to the model.")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip event folders that already contain summary.dat.")
    args = parser.parse_args()

    model = args.model.resolve()
    events_dir = args.events_dir.resolve()
    if not model.is_file():
        print(f"FATAL: Model script not found: {model}", file=sys.stderr)
        return 2
    if not events_dir.is_dir():
        print(f"FATAL: Events directory not found: {events_dir}", file=sys.stderr)
        return 2

    event_dirs = sorted(p for p in events_dir.iterdir() if p.is_dir())
    candidates = [p for p in event_dirs if (p / "bolum_in.dat").is_file()]
    missing_decks = [p for p in event_dirs if not (p / "bolum_in.dat").is_file()]

    LOG_DIR.mkdir(exist_ok=True)
    batch_stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_csv = SCRIPT_DIR / f"batch_summary_{batch_stamp}.csv"
    batch_log = LOG_DIR / f"batch_{batch_stamp}.log"

    print("=" * 72)
    print("BOLUM BATCH RUN")
    print(f"Start:       {timestamp()}")
    print(f"Model:       {model}")
    print(f"Events dir:  {events_dir}")
    print(f"Found:       {len(event_dirs)} event folders")
    print(f"Runnable:    {len(candidates)} (input decks present)")
    print(f"Batch log:   {batch_log}")
    print("=" * 72)

    records = []
    batch_start = time.perf_counter()

    with batch_log.open("w", encoding="utf-8") as log:
        log.write(f"BOLUM batch started: {timestamp()}\n")
        log.write(f"Model: {model}\nEvents directory: {events_dir}\n\n")

        for skipped_dir in missing_decks:
            msg = f"[SKIP] {skipped_dir.name}: bolum_in.dat not found"
            print(msg)
            log.write(msg + "\n")

        for idx, event_dir in enumerate(candidates, 1):
            event = event_dir.name
            deck = event_dir / "bolum_in.dat"
            event_log = LOG_DIR / f"{event}_{batch_stamp}.log"

            if args.skip_existing and (event_dir / "summary.dat").exists():
                msg = f"[{idx}/{len(candidates)}] SKIP {event}: summary.dat already exists"
                print(msg)
                log.write(msg + "\n")
                records.append({
                    "event": event, "status": "skipped", "return_code": "",
                    "elapsed_seconds": 0.0, "start": "", "finish": "",
                    "event_log": str(event_log), "message": "summary.dat already exists",
                })
                continue

            cmd = [sys.executable, str(model), event, "--base", str(events_dir)]
            if args.no_plots:
                cmd.append("--no-plots")
            if args.quiet_model:
                cmd.append("--quiet")
            if args.sat_W is not None:
                cmd.extend(["--sat-W", str(args.sat_W)])

            start = timestamp()
            t0 = time.perf_counter()
            msg = f"\n[{idx}/{len(candidates)}] START {event} | {start}"
            print(msg)
            log.write(msg + "\n")
            log.write("Command: " + " ".join(cmd) + "\n")
            log.flush()

            status = "success"
            return_code = None
            error_msg = ""
            try:
                # Capture model output per-event to avoid interleaving, then display it.
                proc = subprocess.run(
                    cmd, cwd=str(SCRIPT_DIR), text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    check=False,
                )
                return_code = proc.returncode
                output = proc.stdout or ""
                event_log.write_text(output, encoding="utf-8")
                if output:
                    print(output.rstrip())
                    log.write(output)
                if return_code != 0:
                    status = "failed"
                    error_msg = f"model exited with status {return_code}"
            except Exception as exc:
                status = "failed"
                error_msg = f"{type(exc).__name__}: {exc}"
                output = error_msg + "\n"
                event_log.write_text(output, encoding="utf-8")
                print(f"[ERROR] {error_msg}")
                log.write(output)

            elapsed = time.perf_counter() - t0
            finish = timestamp()
            msg = f"[{status.upper()}] {event} | elapsed {elapsed:.2f} s | exit={return_code} | {finish}"
            print(msg)
            log.write(msg + "\n")
            log.flush()

            records.append({
                "event": event, "status": status, "return_code": return_code,
                "elapsed_seconds": round(elapsed, 3), "start": start,
                "finish": finish, "event_log": str(event_log), "message": error_msg,
            })

    total_elapsed = time.perf_counter() - batch_start
    successes = sum(r["status"] == "success" for r in records)
    failures = sum(r["status"] == "failed" for r in records)
    skips = sum(r["status"] == "skipped" for r in records)

    fields = ["event", "status", "return_code", "elapsed_seconds", "start", "finish", "event_log", "message"]
    with summary_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)

    summary = [
        "", "=" * 72, "BATCH SUMMARY",
        f"Finished: {timestamp()}",
        f"Events with input decks: {len(candidates)}",
        f"Successful: {successes}",
        f"Failed:    {failures}",
        f"Skipped:   {skips}",
        f"Total batch elapsed: {total_elapsed:.2f} s ({total_elapsed / 60:.2f} min)",
        f"CSV summary: {summary_csv}",
        f"Batch log:   {batch_log}",
        "=" * 72,
    ]
    summary.extend(f"{r['event']:<30} {r['status']:<10} {r['elapsed_seconds']:>9.2f} s"
                   for r in records)
    summary_text = "\n".join(summary)
    print(summary_text)
    with batch_log.open("a", encoding="utf-8") as log:
        log.write(summary_text + "\n")

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
