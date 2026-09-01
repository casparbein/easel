#!/usr/bin/env python3
"""How much wall-clock does checkpoint re-planning cost a run?

Snakemake's scheduler loop does two expensive things per tick, in this order
(snakemake/scheduling/job_scheduler.py):

    with self._lock:
        self._finish_jobs()     # 1. finish jobs. For a CHECKPOINT job this
                                #    re-plans the ENTIRE DAG (dag.py: finish ->
                                #    update_checkpoint_dependencies ->
                                #    postprocess), logging "Updating checkpoint
                                #    dependencies.". Inherent to gating on a
                                #    checkpoint.
    ...
    run = self.job_selector(needrun)  # 2. pick what to submit. The ONLY step
                                      #    --scheduler controls. Fix:
                                      #    'scheduler: greedy' or
                                      #    'scheduler-subsample: N'.

Both stall submission, so both look identical from the outside -- jobs draining
away and no new ones going out.

This measures (1) and deliberately does not try to measure (2). Snakemake
timestamps its job events but not the phase boundaries, so the interval between
two timestamps contains job runtime, re-planning and selection all mixed
together. What makes (1) measurable anyway is a control: compare the gap after a
job completion that triggered a re-plan against the gap after one that did not.
Selection happens in both cases, so its cost cancels; the difference is the
re-plan. That same cancellation is why (2) cannot be read off one log -- to test
it, set 'scheduler: greedy' in the profile and compare two runs.

    diagnose_scheduler.py .snakemake/log/<run>.snakemake.log
"""
import re
import sys
import time
from collections import Counter

TS = re.compile(r"^\[(\w{3} \w{3} [ \d]\d \d\d:\d\d:\d\d \d{4})\]\s*$")
CHECKPOINT = "Updating checkpoint dependencies."
SELECT = "Select jobs to execute..."
EXECUTE = re.compile(r"^Execute (\d+) jobs?\.\.\.")
FINISHED = re.compile(r"^Finished jobid: \d+ \(Rule: (\S+)\)")


def parse(path):
    """-> ordered list of (epoch_or_None, kind, detail).

    kind is one of: ts, checkpoint, selector, execute, finished.
    """
    events = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            m = TS.match(line)
            if m:
                events.append((time.mktime(time.strptime(m.group(1))), "ts", ""))
                continue
            stripped = line.strip()
            if stripped == CHECKPOINT:
                events.append((None, "checkpoint", ""))
            elif stripped == SELECT:
                events.append((None, "selector", ""))
                continue
            m = EXECUTE.match(line)
            if m:
                events.append((None, "execute", m.group(1)))
                continue
            m = FINISHED.match(line)
            if m:
                events.append((None, "finished", m.group(1)))
    return events


def gaps_after_completions(events):
    """-> (gaps_with_replan, gaps_without).

    For each "Finished jobid" line, the time until the next timestamp. Split by
    whether a checkpoint re-plan was logged directly after that completion.
    """
    ## The timestamp of a job event is printed immediately BEFORE it, so the
    ## timestamp in force at a completion is the last one seen.
    with_replan, without = [], []
    last_ts = None
    for i, (ts, kind, _) in enumerate(events):
        if kind == "ts":
            last_ts = ts
            continue
        if kind != "finished" or last_ts is None:
            continue
        replanned = any(events[j][1] == "checkpoint"
                        for j in range(i + 1, min(i + 3, len(events))))
        nxt = next((events[j][0] for j in range(i + 1, len(events))
                    if events[j][1] == "ts"), None)
        if nxt is None or nxt < last_ts:
            continue
        (with_replan if replanned else without).append(nxt - last_ts)
    return with_replan, without


def describe(label, gaps):
    if not gaps:
        print(f"  {label:<36} (none)")
        return
    gaps = sorted(gaps)
    print(f"  {label:<36} n={len(gaps):>7,}  median {gaps[len(gaps)//2]:>6.1f}s"
          f"  p90 {gaps[int(len(gaps)*0.9)]:>6.1f}s"
          f"  total {sum(gaps)/3600:>6.2f} h")


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    events = parse(sys.argv[1])
    if not events:
        sys.exit("No recognisable snakemake events in that file.")

    stamps = [e[0] for e in events if e[1] == "ts"]
    span = (stamps[-1] - stamps[0]) if len(stamps) > 1 else 0.0
    finished = Counter(e[2] for e in events if e[1] == "finished")
    sizes = [int(e[2]) for e in events if e[1] == "execute"]
    n_replan = sum(1 for e in events if e[1] == "checkpoint")
    n_select = sum(1 for e in events if e[1] == "selector")

    print(f"log span                      : {span/3600:.2f} h")
    print(f"full-DAG re-plans             : {n_replan:,}"
          "   ('Updating checkpoint dependencies.')")
    print(f"job-selection rounds          : {n_select:,}"
          "   ('Select jobs to execute...')")
    print(f"jobs finished                 : {sum(finished.values()):,}")
    if sizes:
        print(f"submission rounds             : {len(sizes):,}"
              f"   jobs per round: mean {sum(sizes)/len(sizes):.1f}, "
              f"min {min(sizes)}, max {max(sizes)}")

    print("\nfinished jobs by rule (top 8):")
    for rule, n in finished.most_common(8):
        print(f"  {rule:<34} {n:>8,}")

    with_replan, without = gaps_after_completions(events)
    print("\ntime from a job finishing to the next timestamped event:")
    describe("completions that re-planned the DAG", with_replan)
    describe("completions that did not (control)", without)

    if not with_replan:
        print("\nNo checkpoint re-plans in this log: either there are no "
              "checkpoints\nleft in the DAG, or none finished during it.")
    elif not without:
        print("\nNo control population, so the re-plan cost cannot be "
              "isolated.")
    else:
        base = sorted(without)[len(without) // 2]
        excess = sum(max(0.0, g - base) for g in with_replan)
        share = (100 * excess / span) if span else 0.0
        print(f"\n  excess over the control median ({base:.1f}s): "
              f"{excess/3600:.2f} h"
              + (f", {share:.0f}% of the log span" if span else ""))
        ## Only say something when there is something to say: on a run whose
        ## re-planning is already cheap, pointing at it would send the reader
        ## after the wrong cost.
        if excess >= 300 and share >= 5:
            print("  That is the scheduler holding its lock in a full-DAG "
                  "postprocess,\n  submitting nothing. It is inherent to "
                  "gating downstream work on a\n  checkpoint verdict -- one "
                  "re-plan per transcript. See 'Large runs'\n  in the README "
                  "for the costs that ARE tunable.")
        else:
            print("  Small enough not to be worth chasing. If throughput still "
                  "sags, the\n  cost is elsewhere -- job selection, or the "
                  "cluster itself.")
        print("\n  This is a LOWER bound, and a rough one: re-planning "
              "completions are\n  mostly fast rules while the control is "
              "whatever else finished, so their\n  own runtimes differ. It is "
              "sound for 'is this worth fixing', not for\n  costing a single "
              "re-plan.")

    print("\n  The job-selection cost is NOT measured here -- it happens in "
          "every tick,\n  so it cancels out of the comparison above. To test "
          "it, put\n  'scheduler: greedy' in the profile and compare two runs.")


if __name__ == "__main__":
    main()
