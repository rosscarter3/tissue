from pyvista import PolyData #type: ignore
import pyvista as pv #type: ignore
import argparse
import init


def print_init(fn: str, nvars: int) -> None:
    m = pv.PolyData(fn)
    t = init.to_init(m=m, nvars=nvars)

    print(t.to_init())


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("fn")
    parser.add_argument("-n", "--nvars", type=int)

    return parser


def command_line_runner():
    parser = get_parser()
    args = parser.parse_args()

    print_init(args.fn, args.nvars)
    

if __name__ == "__main__":
    command_line_runner()
