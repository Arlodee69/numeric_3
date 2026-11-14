import argparse
import importlib


def main():
    parser = argparse.ArgumentParser(description="Лабораторные 3.1–3.5, вариант 2")
    parser.add_argument("lab", choices=["3.1", "3.2", "3.3", "3.4", "3.5", "all"])
    parser.add_argument("--no-plot", action="store_true", help="Только расчёты")
    parser.add_argument("--show", action="store_true", help="Открывать графики в окне")
    args = parser.parse_args()
    if args.no_plot and args.show:
        parser.error("--no-plot и --show нельзя использовать вместе")
    labs = [f"3.{i}" for i in range(1, 6)] if args.lab == "all" else [args.lab]
    for lab in labs:
        name = "lab_" + lab.replace(".", "_")
        options = []
        if lab in ("3.1", "3.2", "3.3"):
            if args.no_plot:
                options.append("--no-plot")
            if args.show:
                options.append("--show")
        importlib.import_module(name).main(options)


if __name__ == "__main__":
    main()
