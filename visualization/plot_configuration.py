#!/usr/bin/python3

import argparse
import asyncio
import fileinput
import json
import sys

from ggb_backend import GeogebraPyppeteerBackend


async def execute_commands(commands, path):
    gb = GeogebraPyppeteerBackend({'headless' : True})
    await gb.run()

    for command in commands:
        await gb.evaluate_command(command)

    await gb.applet_screenshot(path)
    await gb.stop()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", default="example.png")
    parser.add_argument("--json", type=bool, default=False)
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    args = parser.parse_args()

    s = ''.join(args.infile.readlines())
    j = json.loads(s)
    print(j)

    loop = asyncio.get_event_loop()
    loop.run_until_complete(execute_commands(j["commands"], args.path))
    loop.close()

if __name__ == "__main__":
    main()


