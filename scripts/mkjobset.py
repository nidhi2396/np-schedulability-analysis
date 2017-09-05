#!/usr/bin/env python3

import argparse

import csv
import re

import fractions

from collections import namedtuple

def lcm(a,b):
    return abs(a * b) / fractions.gcd(a,b) if a and b else 0

def parse_range(s):
    m = re.fullmatch('uniform params: \[([0-9.]+) to ([0-9.]+)\]', s)
    return (float(m.group(1)), float(m.group(2)))

def parse_jitter(s):
    if s.strip() == 'none':
        return (0.0, 0.0)
    else:
        return parse_range(s)

def parse_task_id(s):
    m = re.fullmatch('T\(([0-9]+) [0-9]+\)', s)
    return int(m.group(1))

Task = namedtuple('Task', ['id', 'period', 'phase', 'relative_deadline', 'cost_range', 'jitter_range'])
Job  = namedtuple('Job', ['id', 'tid', 'arrival_range', 'absolute_deadline', 'cost_range'])

def parse_task_set(fname):
    f = open(fname, 'r')
    data = csv.DictReader(f, skipinitialspace=True)
    for row in data:
        id = parse_task_id(row['Job Name'])
        period = float(row['Period'])
        phase  = float(row['Release Phase'])
        deadline = float(row['Deadline'])
        exec_cost = float(row['Initial Execution Time'])
        wcet_range = parse_range(row['WCET Jitter'])
        jitter_range = parse_jitter(row['Release Jitter'])

        yield Task(id, period, phase, deadline, wcet_range, jitter_range)

def hyperperiod(tasks):
    h = 1
    for t in tasks:
        h = lcm(h, t.period)
    return h

def generate_jobs(t, hyperperiod):
    assert t.phase == 0 # can't yet handle non-zero phase

    rel = 0.0
    id = 1
    while rel < hyperperiod:
        yield Job(id, t.id, (rel + t.jitter_range[0], rel + t.jitter_range[1]),
                  rel + t.relative_deadline, t.cost_range)
        id  += 1
        rel += t.period

def assign_edf_priority(job):
    return job.absolute_deadline

def mk_rate_monotonic(tasks):
    by_period = sorted(tasks, key=lambda t: t.period)
    prios = {}
    for p, t in enumerate(by_period):
        prios[t.id] = p + 1
    return lambda j: prios[j.tid]

def format_job(j, prio_fn):
    #Job<dtime_t> j1{9,  I(0, 0), I(3, 13), 60, 60};
    return "Job({0}, range({1}, {2}), range({3}, {4}), {5}, {6}, {7})".format(
        j.id, j.arrival_range[0], j.arrival_range[1],
        j.cost_range[0], j.cost_range[1], j.absolute_deadline, prio_fn(j), j.tid)

def format_csv(j, prio_fn):
    return "{7:10}, {0:10}, {1:20}, {2:20}, {3:20}, {4:20}, {5:20}, {6:20}".format(
        j.id, j.arrival_range[0], j.arrival_range[1],
        j.cost_range[0], j.cost_range[1], j.absolute_deadline, prio_fn(j), j.tid)

def csv_header():
    return "{7:>10s}, {0:>10s}, {1:>20s}, {2:>20s}, {3:>20s}, {4:>20s}, {5:>20s}, {6:>20s}".format(
        "Job ID", "Arrival min", "Arrival max", "Cost min", "Cost max",
        "Deadline", "Priority", "Task ID")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert task specification to job sets")

    parser.add_argument('input_files', nargs='*',
        help='the task sets that should contribute to this job set')

    parser.add_argument('--rate-monotonic', dest='priority_policy', default='EDF',
                        action='store_const', const='RM', help='how to assing job prios')

    return parser.parse_args()

def main():
    opts = parse_args()

    tasks = []
    for f in opts.input_files:
        tasks += list(parse_task_set(f))

    if opts.priority_policy == 'EDF':
        prio = assign_edf_priority
    elif opts.priority_policy == 'RM':
        prio = mk_rate_monotonic(tasks)
    else:
        assert False # unsupported prio policy

    h = hyperperiod(tasks)

    print(csv_header())
    for t in tasks:
        for j in generate_jobs(t, h):
            print(format_csv(j, prio))


if __name__ == '__main__':
    main()
