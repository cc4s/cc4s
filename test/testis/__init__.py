import argparse
import pathlib
import hashlib
import os
import sys
import yaml
import logging
from collections import namedtuple
import subprocess as sp
import shlex
import urllib.request
import tqdm


__version__ = "0.0.1"
__author__ = "Alejandro Gallo"
__email__ = "aamsgallo@gmail.com"
__license__ = "GPLv3"
STORE_FOLDER = pathlib.Path("~/.testis-store").expanduser()


TestCase = namedtuple("TestCase", "tags name path description runner resources")
Resource = namedtuple("Resource", "url out hash")


def get_tests_in_dir(folder):
    return [pathlib.Path(p)
            for p, _, _ in os.walk(folder)
            if (pathlib.Path(p) / "info.yaml").exists()]


def folder_to_test(folder):
    info = folder / "info.yaml"

    _data = read_yaml(info)

    return TestCase( name=_data.get("name") or folder.name
                   , tags=_data.get("tags", "").split(" ")
                   , path=pathlib.Path(os.path.abspath(folder))
                   , description=_data.get("description", "")
                   , runner=_data.get("runner", "./run-test")
                   , resources=parse_resources(_data.get("resources", []))
                   )


def parse_resources(raw_resources):
    assert isinstance(raw_resources, list)
    return [ Resource( url=r["url"]
                     , out=r["out"]
                     , hash=hash_string(r["url"])
                     ) for r in raw_resources ]


def hash_string(r):
    m = hashlib.md5()
    m.update(r.encode())
    return m.hexdigest()


def get_resource(r):
    out = STORE_FOLDER / r.hash / r.out
    if not out.exists():
        out.parent.mkdir(exist_ok=True)
        logging.info("downloading %s ⇒ %s", r.out, out)
        data = urllib.request.urlopen(r.url).read()
        with open(out, 'wb+') as f:
            f.write(data)


def link_resource(r, basedir):
    i = (STORE_FOLDER / r.hash / r.out).resolve()
    o = basedir / r.out
    if o.exists():
        return
    logging.info("linking %s ⇒ %s", o, i)
    o.symlink_to(i)


def read_yaml(filepath):
    with open(filepath) as f:
        return yaml.load(f, Loader=yaml.Loader) or {}


def run_test(test):
    cmd = shlex.split(test.runner)
    logging.debug("running: %s", cmd)

    for r in test.resources:
        get_resource(r)
        link_resource(r, test.path)

    with sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE) as p:
        p.wait()
        return dict( stdout=p.stdout.read().decode()
                   , stderr=p.stderr.read().decode()
                   , returncode=p.returncode
                   )


def call(cmd):
    assert isinstance(cmd, str)
    cmd = shlex.split(cmd.format(**os.environ))
    with sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE) as p:
        p.wait()
        ret = p.returncode
        if p.returncode != 0:
            out = p.stdout.read().decode()
            err = p.stdout.read().decode()
            raise Exception(out + "\n" + err)


def tail(xs, n):
    return xs[-min(len(xs), n): -1]


def main():
    parser = argparse.ArgumentParser("testis", "A minimal test runner")
    parser.add_argument("folders",
                        help="Test folder where to find tests",
                        nargs="+")
    parser.add_argument("-d", help="Debug mode", action="store_true")
    parser.add_argument("--list-tags", help="List tags", action="store_true")
    parser.add_argument(
        "--tail",
        help="How many lines to see from output when error occurs",
        type=int,
        default=10)
    args = parser.parse_args()

    STORE_FOLDER.mkdir(exist_ok=True)

    if args.d:
        logging.warning("verbose mode")
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    # get only folders
    args.folders = [pathlib.Path(f) for f in args.folders if os.path.isdir(f)]

    logging.info("Looking for tests in %s folders", len(args.folders))

    tests = []
    logging.info("Finding tests")
    for f in args.folders:
        logging.debug("in %s", f)
        tests.extend(map(folder_to_test, get_tests_in_dir(f)))

    logging.info("Found %s test folders", len(tests))

    if args.list_tags:
        for t in set(sum([t.tags for t in tests], [])):
            print(t)
        return

    cwd = pathlib.Path.cwd()
    logging.info("Running tests")
    for test in tests:
        os.chdir(test.path)
        print("∷", test.name)
        result = run_test(test)
        if result["returncode"] != 0:
            print("\x1b[31m\t[X]\x1b[0m")
            for f in ["stdout", "stderr"]:
                out = ["\t\x1b[31m» ({})\x1b[0m  {}".format(f, l)
                        for l in result[f].split("\n")]
                print("\n".join(["\t...."] + tail(out, args.tail)))
        else:
            print("\x1b[32m\t[ok]\x1b[0m")
        os.chdir(cwd)
