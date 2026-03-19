import argparse

from iris.features import run_feature_generation
from iris.predict import run_prediction


def main() -> None:
    parser = argparse.ArgumentParser(
        description="IRIS: RNA-ligand docking feature generation and pose re-ranking"
    )
    subparsers = parser.add_subparsers(dest="command")

    features_parser = subparsers.add_parser(
        "features",
        help="Generate IRIS features for a single complex directory",
    )
    features_parser.add_argument(
        "input_dir",
        help="Path to the input complex directory",
    )
    features_parser.add_argument(
        "pymol_path",
        help="Path to the PyMOL executable",
    )
    features_parser.add_argument(
        "pymol_plugin_path",
        help="Path to the PyMOL plugin directory",
    )

    predict_parser = subparsers.add_parser(
        "predict",
        help="Run IRIS pose re-ranking on a single complex directory",
    )
    predict_parser.add_argument(
        "input_dir",
        help="Path to the input complex directory",
    )
    predict_parser.add_argument(
        "--model-path",
        default=None,
        help="Optional path to a custom trained IRIS model (.pkl)",
    )

    args = parser.parse_args()

    if args.command == "features":
        run_feature_generation(
            args.input_dir,
            args.pymol_path,
            args.pymol_plugin_path,
        )
        print("Features generated successfully!")
        return

    if args.command == "predict":
        output_path = run_prediction(args.input_dir, model_path=args.model_path)
        print(f"Processed and wrote sorted SDF file to {output_path}")
        return

    parser.print_help()


if __name__ == "__main__":
    main()
