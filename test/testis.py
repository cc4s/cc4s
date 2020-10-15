#!/usr/bin/env python3
import argparse
import hashlib
import os
import os.path as op
import sys
import logging
from collections import namedtuple
import subprocess as sp
import shlex
import urllib.request
import json
import functools
import re


def get_store_folder():
    f = __file__
    while os.path.islink(f):
        f = os.readlink(f)
    return op.join(op.dirname(op.abspath(f)), "testis-store")


__version__ = "0.0.1"
__author__ = "Alejandro Gallo"
__email__ = "aamsgallo@gmail.com"
__license__ = "GPLv3"
STORE_FOLDER = get_store_folder()
INFO_FILE_NAME = "test.json"
DEFAULT_RUN_SCRIPT = "./run.py"
DEFAULT_CHECK_SCRIPT = "./check.py"
CLEAR = "\x1b[0m"
GREEN = "\x1b[32m"
RED = "\x1b[31m"
MAGENTA = "\x1b[35m"

TestCase = namedtuple("TestCase", "tags name path outpath "
                                  "description run check resources")
Resource = namedtuple("Resource", "url out hash")


def get_tests_in_dir(folder):
    return [p for p, _, _ in os.walk(folder)
              if op.exists(op.join(p, INFO_FILE_NAME))]


def folder_to_test(folder, outname):
    info = op.join(folder, INFO_FILE_NAME)

    with open(info) as f:
        _data = json.load(f)

    path = op.abspath(folder)
    outpath = op.join(path, outname)

    return TestCase( name=_data.get("name") or os.basename(folder)
                   , tags=_data.get("tags", "").split(" ")
                   , path=path
                   , outpath=outpath
                   , description=_data.get("description", "")
                   , run=_data.get("run", DEFAULT_RUN_SCRIPT)
                   , check=_data.get("check", DEFAULT_CHECK_SCRIPT)
                   , resources=parse_resources(_data.get("resources", []))
                   )


def parse_resources(raw_resources):
    assert isinstance(raw_resources, list)
    return [ Resource( url=r["url"]
                     , out=r["out"]
                     , hash=hash_string(r["url"])
                     ) for r in raw_resources ]


def hash_string(r):
    assert isinstance(r, str)
    m = hashlib.md5()
    m.update(r.encode())
    return m.hexdigest()


def get_resource(r):
    assert isinstance(r, Resource)
    out = op.join(STORE_FOLDER, r.hash, r.out)
    if not op.exists(out):
        os.makedirs(op.dirname(out), exist_ok=True)
        logging.info("%sdownlo.%s %s => %s",
                     MAGENTA, CLEAR, r.out, op.relpath(out))
        data = urllib.request.urlopen(r.url).read()
        with open(out, 'wb+') as f:
            f.write(data)


def link_resource(r, basedir):
    assert isinstance(r, Resource)
    assert isinstance(basedir, str)
    i = op.abspath(op.join(STORE_FOLDER, r.hash, r.out))
    o = op.join(basedir, r.out)
    if op.exists(o):
        return
    logging.info("%ssymlink%s %s <= %s", MAGENTA, CLEAR, r.out, op.relpath(i))
    os.symlink(i, o)


def read_yaml(filepath):
    assert isinstance(filepath, str)
    try:
        import yaml
    except ImportError:
        raise Exception("Your python can not find the python package"
                        "pyyaml! Ask your system administrator to install it")

    with open(filepath) as f:
        return yaml.load(f, Loader=yaml.Loader) or {}


def get_and_link_resources(test):
    assert isinstance(test, TestCase)
    for r in test.resources:
        get_resource(r)
        link_resource(r, test.outpath)


def link_inputfiles_to_outputpath(test):
    assert isinstance(test, TestCase)
    os.makedirs(test.outpath, exist_ok=True)
    for f in os.listdir(test.path):
        # do not link the info file name so that it is not recognised
        # as a test
        if op.basename(f) == INFO_FILE_NAME:
            continue
        out = op.join(test.outpath, op.basename(f))
        inf = op.join(test.path, op.basename(f))
        if not op.exists(out):
            os.symlink(inf, out)


def run_in_test(test, script_file):
    assert isinstance(test, TestCase)
    assert isinstance(script_file, str)
    cmd = shlex.split(script_file)
    logging.debug("running: %s", cmd)
    cwd = os.getcwd()

    logging.debug("chdir: %s", test.outpath)
    os.chdir(test.outpath)

    with sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE) as p:
        p.wait()
        os.chdir(cwd)
        return dict( stdout=p.stdout.read().decode()
                   , stderr=p.stderr.read().decode()
                   , returncode=p.returncode
                   )


def call(cmd):
    assert isinstance(cmd, str)
    try:
        cmd = shlex.split(cmd.format(**os.environ))
    except KeyError as e:
        raise Exception("Environment variable {}{}{} not known,"
                        " please set it and run the test again"
                        .format(RED, e, CLEAR))

    with sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE) as p:
        p.wait()
        ret = p.returncode
        out = p.stdout.read().decode()
        err = p.stdout.read().decode()
        if p.returncode != 0:
            raise Exception(out + "\n" + err)
        else:
            sys.stdout.write(out)
            sys.stderr.write(err)


def tail(xs, n):
    return xs[-min(len(xs), n): -1]


def filter_tags(tests, tags):
    tags_code = tags
    tags_code = re.sub(r".and.", r"&", tags_code)
    tags_code = re.sub(r".or.", r"|", tags_code)
    tags_code = re.sub(r"(\w+)", r"('\1' in __x_)", tags_code)
    tags_code = re.sub(r"\&", r"and", tags_code)
    tags_code = re.sub(r"\|", r"or", tags_code)
    tags_code = re.sub(r"^", r"lambda __x_: ", tags_code)
    tags_lambda = eval(tags_code)
    logging.info("tags λ: %s", tags_code)

    return [t for t in tests if tags_lambda(t.tags)]


def main():
    parser = argparse.ArgumentParser("testis", "A minimal test run")
    parser.add_argument("folders",
                        help="Test folder where to find tests",
                        nargs="+")
    parser.add_argument("-d", help="Debug mode", action="store_true")
    parser.add_argument("--list-tags", help="List tags", action="store_true")
    parser.add_argument("-t", "--tags", default=None, type=str,
                        help="Run only tags matching, for instance: "
                             "--tags 'essential .and. (mem .or. ccsd)'")
    parser.add_argument("--rx", type=str,
                        help="Run tests matching regular expression",
                        default=None)
    parser.add_argument("-r", "--run",
                        help="Just run the 'run' phase",
                        action="store_true",
                        default=None)
    parser.add_argument("-c", "--check",
                        help="Just run the 'check' phase",
                        action="store_true")
    parser.add_argument("-n", "--name",
                        help="Name for the results folder of the test",
                        default="test-results",
                        type=str)
    parser.add_argument(
        "--tail",
        help="How many lines to see from output when error occurs",
        type=int,
        default=10)
    args = parser.parse_args()

    os.makedirs(STORE_FOLDER, exist_ok=True)

    if args.d:
        logging.basicConfig(
            level=logging.DEBUG,
                format="(%(relativeCreated)d) %(levelname)s:%(message)s")
    else:
        logging.basicConfig(level=logging.INFO,
                            format="(%(relativeCreated)d) %(message)s")

    # get only folders
    args.folders = [f for f in args.folders if os.path.isdir(f)]

    logging.info("Looking for tests in %s folder(s)", len(args.folders))

    tests = []
    logging.info("Finding tests")
    for f in args.folders:
        logging.debug("in %s", f)
        tests.extend(map(functools.partial(folder_to_test, outname=args.name),
                         get_tests_in_dir(f)))

    logging.info("Found %s test folders", len(tests))
    all_tags = set(sum([t.tags for t in tests], []))

    if args.list_tags:
        for t in all_tags:
            print(t)
        return

    if args.rx:
        rx = re.compile(args.rx)
        tests = [t for t in tests if rx.match(t.name)]
        logging.info("restrict to %s tests due to regex", len(tests))

    if args.tags:
        tests = filter_tags(tests, args.tags)
        logging.info("restrict to %s tests due to tags", len(tests))

    if args.run: logging.info("will run   %s test", len(tests))
    if args.run: logging.info("will check %s test", len(tests))



    cwd = os.getcwd()
    logging.info("Running tests (in {}{}{})".format(MAGENTA, args.name, CLEAR))
    for test in tests:
        logging.info("{}∷{} {}".format(GREEN, CLEAR, test.name))

        logging.debug("linking input files")
        link_inputfiles_to_outputpath(test)
        logging.debug("getting and linking resources to output path")
        get_and_link_resources(test)

        os.chdir(test.outpath)

        for script, doit in [(test.run, args.run), (test.check, args.check)]:

            if not doit:
                continue

            result = run_in_test(test, script_file=script)
            if result["returncode"] != 0:
                print("{}\t[X]{}".format(RED, CLEAR))
                for f in ["stdout", "stderr"]:
                    out = ["\t{}» ({}){}  {}".format(RED, f, CLEAR, l)
                            for l in result[f].split("\n")]
                    print("\n".join(["\t...."] + tail(out, args.tail)))
            else:
                print("{}\t[ok]{} ({})".format(GREEN, CLEAR, script))

            for content in ["stdout", "stderr"]:
                fname = script + "." + content
                logging.debug("writing %s", fname)
                with open(fname, "w+") as f:
                    f.write(result[content])


        os.chdir(cwd)


if __name__ == "__main__":
    main()
